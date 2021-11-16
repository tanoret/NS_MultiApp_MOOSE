//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//#include "Transient.h"
#include "CustomTransient.h"

// MOOSE includes
#include "Factory.h"
#include "SubProblem.h"
#include "TimeStepper.h"
#include "MooseApp.h"
#include "Conversion.h"
#include "FEProblem.h"
#include "NonlinearSystem.h"
#include "Control.h"
#include "TimePeriod.h"
#include "MooseMesh.h"
#include "TimeIntegrator.h"
#include "Console.h"
#include "INSFVPressureVariable.h"

#include "libmesh/implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

#include "AuxiliarySystem.h"
#include "NonlinearSystem.h"

// C++ Includes
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <petscsnes.h>

registerMooseObject("MooseApp", CustomTransient);

defineLegacyParams(CustomTransient);

InputParameters
CustomTransient::validParams()
{
  InputParameters params = Transient::validParams();
  params.addClassDescription("Executioner for time varying simulations.");

  params += FEProblemSolve::validParams();

  std::vector<Real> sync_times(1);
  sync_times[0] = -std::numeric_limits<Real>::max();

  /**
   * For backwards compatibility we'll allow users to set the TimeIntegration scheme inside of the
   * executioner block
   * as long as the TimeIntegrator does not have any additional parameters.
   */
  MooseEnum schemes("implicit-euler explicit-euler crank-nicolson bdf2 explicit-midpoint dirk "
                    "explicit-tvd-rk-2 newmark-beta",
                    "implicit-euler");

  params.addParam<Real>("start_time", 0.0, "The start time of the simulation");
  params.addParam<Real>("end_time", 1.0e30, "The end time of the simulation");
  params.addParam<Real>("dt", 1., "The timestep size between solves");
  params.addParam<Real>("dtmin", 2.0e-14, "The minimum timestep size in an adaptive run");
  params.addParam<Real>("dtmax", 1.0e30, "The maximum timestep size in an adaptive run");
  params.addParam<bool>(
      "reset_dt", false, "Use when restarting a calculation to force a change in dt.");
  params.addParam<unsigned int>("num_steps",
                                std::numeric_limits<unsigned int>::max(),
                                "The number of timesteps in a transient run");
  params.addParam<int>("n_startup_steps", 0, "The number of timesteps during startup");

  params.addDeprecatedParam<bool>("trans_ss_check",
                                  false,
                                  "Whether or not to check for steady state conditions",
                                  "Use steady_state_detection instead");
  params.addDeprecatedParam<Real>("ss_check_tol",
                                  1.0e-08,
                                  "Whenever the relative residual changes by less "
                                  "than this the solution will be considered to be "
                                  "at steady state.",
                                  "Use steady_state_tolerance instead");
  params.addDeprecatedParam<Real>(
      "ss_tmin",
      0.0,
      "Minimum amount of time to run before checking for steady state conditions.",
      "Use steady_state_start_time instead");

  params.addParam<bool>(
      "steady_state_detection", false, "Whether or not to check for steady state conditions");
  params.addParam<Real>("steady_state_tolerance",
                        1.0e-08,
                        "Whenever the relative residual changes by less "
                        "than this the solution will be considered to be "
                        "at steady state.");
  params.addParam<Real>(
      "steady_state_start_time",
      0.0,
      "Minimum amount of time to run before checking for steady state conditions.");
  params.addParam<bool>(
      "normalize_solution_diff_norm_by_dt",
      true,
      "Whether to divide the solution difference norm by dt. If taking 'small' "
      "time steps you probably want this to be true. If taking very 'large' timesteps in an "
      "attempt to *reach* a steady-state, you probably want this parameter to be false.");

  params.addParam<std::vector<std::string>>("time_periods", "The names of periods");
  params.addParam<std::vector<Real>>("time_period_starts", "The start times of time periods");
  params.addParam<std::vector<Real>>("time_period_ends", "The end times of time periods");
  params.addParam<bool>(
      "abort_on_solve_fail", false, "abort if solve not converged rather than cut timestep");
  params.addParam<bool>(
      "error_on_dtmin",
      true,
      "Throw error when timestep is less than dtmin instead of just aborting solve.");
  params.addParam<MooseEnum>("scheme", schemes, "Time integration scheme used.");
  params.addParam<Real>("timestep_tolerance",
                        2.0e-14,
                        "the tolerance setting for final timestep size and sync times");

  params.addParam<bool>("use_multiapp_dt",
                        false,
                        "If true then the dt for the simulation will be "
                        "chosen by the MultiApps.  If false (the "
                        "default) then the minimum over the master dt "
                        "and the MultiApps is used");

  params.addParamNamesToGroup(
      "steady_state_detection steady_state_tolerance steady_state_start_time",
      "Steady State Detection");

  params.addParamNamesToGroup("start_time dtmin dtmax n_startup_steps trans_ss_check ss_check_tol "
                              "ss_tmin abort_on_solve_fail timestep_tolerance use_multiapp_dt",
                              "Advanced");

  params.addParamNamesToGroup("time_periods time_period_starts time_period_ends", "Time Periods");

  params.addParam<bool>("verbose_print", false, "If true then print the matrix of coefs and rhs.");

  params.addParam<bool>("momentum_predictor_bool", false, "If true then do SIMPLE transfers.");

  return params;
}

CustomTransient::CustomTransient(const InputParameters & parameters)
  : Transient(parameters),
    _problem(_fe_problem),
    _feproblem_solve(*this),
    _nl(_fe_problem.getNonlinearSystemBase()),
    _time_scheme(getParam<MooseEnum>("scheme").getEnum<Moose::TimeIntegratorType>()),
    _t_step(_problem.timeStep()),
    _time(_problem.time()),
    _time_old(_problem.timeOld()),
    _dt(_problem.dt()),
    _dt_old(_problem.dtOld()),
    _unconstrained_dt_custom(declareRecoverableData<Real>("unconstrained_dt_custom", -1)),
    _at_sync_point_custom(declareRecoverableData<bool>("at_sync_point_custom", false)),
    _last_solve_converged_custom(declareRecoverableData<bool>("last_solve_converged_custom", true)),
    _xfem_repeat_step(false),
    _end_time(getParam<Real>("end_time")),
    _dtmin(getParam<Real>("dtmin")),
    _dtmax(getParam<Real>("dtmax")),
    _num_steps(getParam<unsigned int>("num_steps")),
    _n_startup_steps(getParam<int>("n_startup_steps")),
    _steady_state_detection(getParam<bool>("steady_state_detection")),
    _steady_state_tolerance(getParam<Real>("steady_state_tolerance")),
    _steady_state_start_time(getParam<Real>("steady_state_start_time")),
    _sync_times(_app.getOutputWarehouse().getSyncTimes()),
    _abort(getParam<bool>("abort_on_solve_fail")),
    _error_on_dtmin(getParam<bool>("error_on_dtmin")),
    _time_interval_custom(declareRecoverableData<bool>("time_interval_custom", false)),
    _start_time(getParam<Real>("start_time")),
    _timestep_tolerance(getParam<Real>("timestep_tolerance")),
    _target_time_custom(declareRecoverableData<Real>("target_time_custom", -std::numeric_limits<Real>::max())),
    _use_multiapp_dt(getParam<bool>("use_multiapp_dt")),
    _solution_change_norm_custom(declareRecoverableData<Real>("solution_change_norm_custom", 0.0)),
    _sln_diff(_nl.addVector("sln_diff", false, PARALLEL)),
    _normalize_solution_diff_norm_by_dt(getParam<bool>("normalize_solution_diff_norm_by_dt")),
    _verbose_print(getParam<bool>("verbose_print")),
    _momentum_predictor_bool(getParam<bool>("momentum_predictor_bool"))
{
  _fixed_point_solve->setInnerSolve(_feproblem_solve);

  // Handle deprecated parameters
  if (!parameters.isParamSetByAddParam("trans_ss_check"))
    _steady_state_detection = getParam<bool>("trans_ss_check");

  if (!parameters.isParamSetByAddParam("ss_check_tol"))
    _steady_state_tolerance = getParam<Real>("ss_check_tol");

  if (!parameters.isParamSetByAddParam("ss_tmin"))
    _steady_state_start_time = getParam<Real>("ss_tmin");

  _t_step = 0;
  _dt = 0;
  _next_interval_output_time = 0.0;

  // Either a start_time has been forced on us, or we want to tell the App about what our start time
  // is (in case anyone else is interested.
  if (_app.hasStartTime())
    _start_time = _app.getStartTime();
  else if (parameters.isParamSetByUser("start_time"))
    _app.setStartTime(_start_time);

  _time = _time_old = _start_time;
  _problem.transient(true);

  setupTimeIntegrator();

  if (_app.halfTransient()) // Cut timesteps and end_time in half...
  {
    _end_time /= 2.0;
    _num_steps /= 2.0;

    if (_num_steps == 0) // Always do one step in the first half
      _num_steps = 1;
  }

}

void
CustomTransient::init()
{
  if (!_time_stepper.get())
  {
    InputParameters pars = _app.getFactory().getValidParams("ConstantDT");
    pars.set<SubProblem *>("_subproblem") = &_problem;
    pars.set<Transient *>("_executioner") = static_cast<Transient *>(this);

    /**
     * We have a default "dt" set in the Transient parameters but it's possible for users to set
     * other parameters explicitly that could provide a better calculated "dt". Rather than provide
     * difficult to understand behavior using the default "dt" in this case, we'll calculate "dt"
     * properly.
     */
    if (!_pars.isParamSetByAddParam("end_time") && !_pars.isParamSetByAddParam("num_steps") &&
        _pars.isParamSetByAddParam("dt"))
      pars.set<Real>("dt") = (getParam<Real>("end_time") - getParam<Real>("start_time")) /
                             static_cast<Real>(getParam<unsigned int>("num_steps"));
    else
      pars.set<Real>("dt") = getParam<Real>("dt");

    pars.set<bool>("reset_dt") = getParam<bool>("reset_dt");
    _time_stepper = _app.getFactory().create<TimeStepper>("ConstantDT", "TimeStepper", pars);
  }

  _problem.execute(EXEC_PRE_MULTIAPP_SETUP);
  _problem.initialSetup();

  /**
   * If this is a restart run, the user may want to override the start time, which we already set in
   * the constructor. "_time" however will have been "restored" from the restart file. We need to
   * honor the original request of the developer now that the restore has been completed. This must
   * occur before we init the time stepper (since that prints out the start time). The multiapp case
   * is also bit complicated. If we didn't set a start time, the app won't have it yet, so we just
   * restart the old time from the current time.
   */
  if (_app.isRestarting())
  {
    if (_app.hasStartTime())
      _time = _time_old = _app.getStartTime();
    else
      _time_old = _time;
  }

  _time_stepper->init();

  if (_app.isRecovering()) // Recover case
  {
    if (_t_step == 0)
      mooseError("Internal error in Transient executioner: _t_step is equal to 0 while recovering "
                 "in init().");

    _dt_old = _dt;
  }
}

void
CustomTransient::preExecute()
{
  _time_stepper->preExecute();

  if (!_app.isRecovering())
  {
    _t_step = 0;
    _dt = 0;
    _next_interval_output_time = 0.0;
    if (!_app.isRestarting())
      _time = _time_old = _start_time;

    _problem.outputStep(EXEC_INITIAL);

    computeDT();
    _dt = getDT();
    if (_dt == 0)
      mooseError("Time stepper computed zero time step size on initial which is not allowed.\n"
                 "1. If you are using an existing time stepper, double check the values in your "
                 "input file or report an error.\n"
                 "2. If you are developing a new time stepper, make sure that initial time step "
                 "size in your code is computed correctly.");
    _nl.getTimeIntegrator()->init();
  }
}

void
CustomTransient::preStep()
{
  _time_stepper->preStep();
}

void
CustomTransient::postStep()
{
  std::cout << "Entering post step momentum predictor." << std::endl;

  _time_stepper->postStep();

  // Little PetSc obbejcts

  // Petsc primitive data types
  // PetscInt n = 20;
  // PetscPrintf(PETSC_COMM_WORLD, "n = %d \n", n);
  // PetscScalar v = -3.5, w = 1e9;
  // PetscPrintf(PETSC_COMM_WORLD, "v = %f - w = %e \n", v, w);
  // PetscReal x = 1e-9;
  // PetscPrintf(PETSC_COMM_WORLD, "x = %e \n", x);

  // Petsc parallel handling
  // PetscErrorCode rank, size;
  // MPI_Comm_size(PETSC_COMM_WORLD, &size);
  // MPI_Comm_size(PETSC_COMM_WORLD, &rank);
  // PetscPrintf(PETSC_COMM_WORLD,"Number of processors = %d, rank = %d\n",size,rank);

  // Vector handling
  // MPI_Comm comm = MPI_COMM_WORLD;
  // Vec xvec, yvec;
  // PetscInt vec_rank = 10;
  // VecCreate(comm, &xvec);
  // VecSetSizes(xvec, vec_rank+1, vec_rank+1);
  // VecSetType(xvec, VECMPI);
  // VecSet(xvec, 1.0);
  // std::cout << "xvec: " << std::endl; VecView(xvec, PETSC_VIEWER_STDOUT_WORLD);
  // for(PetscInt i=0; i <= vec_rank; i++)
  // {
  //   v = i;
  //   VecSetValues(xvec, 1, &i, &v, INSERT_VALUES);
  // }
  // VecAssemblyBegin(xvec);
  // VecAssemblyEnd(xvec);
  // std::cout << "xvec mod: " << std::endl; VecView(xvec, PETSC_VIEWER_STDOUT_WORLD);
  // VecDuplicate(xvec, &yvec);
  // PetscInt indexes_to_get[] = {2,7};
  // PetscScalar storing_values[2];
  // for(PetscInt i = 0; i < 2; i++)
  // {
  //   storing_values[i] = VecGetValues(yvec, 1, &indexes_to_get[i], &storing_values[i]);
  //   PetscPrintf(comm, "%d \n", storing_values[i]);
  //   VecSetValues(xvec, 1, &indexes_to_get[i], &storing_values[i], INSERT_VALUES);
  // }
  // VecAssemblyBegin(xvec);
  // VecAssemblyEnd(xvec);
  // std::cout << "xvec mod 2: " << std::endl; VecView(xvec, PETSC_VIEWER_STDOUT_WORLD);
  // VecDestroy(&xvec);
  // VecDestroy(&yvec);

  // Getting handlers
  PetscInt mesh_dimension  = feProblem().mesh().dimension();
  PetscInt active_elements = static_cast<PetscInt>(feProblem().mesh().getMesh().n_active_elem());
  if(_verbose_print)
  {
    std::cout << "Mesh dimensions: " << mesh_dimension << std::endl;
    std::cout << "Mesh DOFS per dimension: " << active_elements << std::endl;
  }

  // Set Viewer printing format.
  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_DENSE);
  System & sys = _nl.system();
  NonlinearImplicitSystem & isys = dynamic_cast<NonlinearImplicitSystem&>(sys);
  SparseMatrix<Number> * mat = isys.matrix;
  PetscMatrix<Number> * pmat = dynamic_cast<PetscMatrix<Number> *>(mat);
  if(_verbose_print)
  {
    std::cout << "Matrix of coefs: " << std::endl;
    MatView(pmat->mat(), PETSC_VIEWER_STDOUT_WORLD);
  }

  // Vec _Ainv_x, _Hu_x, _rhs_x;
  // Vec _Ainv_y, _Hu_y, _rhs_y;
  // Vec _Ainv_z, _Hu_z, _rhs_z;

  // Vectors to store Ainv and H*u
  // MPI_Comm comm = MPI_COMM_WORLD;
  // if(mesh_dimension > 0)
  // {
  //   // Create vectors
  //   VecCreate(comm, &_Ainv_x);
  //   VecCreate(comm, &_Hu_x);
  //   VecCreate(comm, &_rhs_x);
  //
  //   // Set Sizes
  //   VecSetSizes(_Ainv_x, active_elements, PETSC_DECIDE);
  //   VecSetSizes(_Hu_x, active_elements, PETSC_DECIDE);
  //   VecSetSizes(_rhs_x, active_elements, PETSC_DECIDE);
  //
  //   // Set types
  //   VecSetFromOptions(_Ainv_x);
  //   VecSetFromOptions(_Hu_x);
  //   VecSetFromOptions(_rhs_x);
  //
  //   // Assign zero to all entries
  //   VecZeroEntries(_Ainv_x);
  //   VecZeroEntries(_Hu_x);
  //   VecZeroEntries(_rhs_x);
  //
  //   // Print if verbose
  //   if(_verbose_print)
  //   {
  //     std::cout << "Ainv_x: " << std::endl;
  //     VecView(_Ainv_x, PETSC_VIEWER_STDOUT_WORLD);
  //     std::cout << "Hu_x: " << std::endl;
  //     VecView(_Hu_x, PETSC_VIEWER_STDOUT_WORLD);
  //     std::cout << "rhs_x: " << std::endl;
  //     VecView(_rhs_x, PETSC_VIEWER_STDOUT_WORLD);
  //   }
  // }
  //
  // if(mesh_dimension > 1)
  // {
  //   // Create and populate vectors by copy
  //   VecDuplicate(_Ainv_x, &_Ainv_y);
  //   VecDuplicate(_Hu_x, &_Hu_y);
  //   VecDuplicate(_rhs_x, &_rhs_y);
  //
  //   // Print if verbose
  //   if(_verbose_print)
  //   {
  //     std::cout << "Ainv_y: " << std::endl;
  //     VecView(_Ainv_y, PETSC_VIEWER_STDOUT_WORLD);
  //     std::cout << "Hu_y: " << std::endl;
  //     VecView(_Hu_y, PETSC_VIEWER_STDOUT_WORLD);
  //     std::cout << "rhs_y: " << std::endl;
  //     VecView(_rhs_y, PETSC_VIEWER_STDOUT_WORLD);
  //   }
  // }
  //
  // if(mesh_dimension > 2)
  // {
  //   // Create and populate vectors by copy
  //   VecDuplicate(_Ainv_x, &_Ainv_z);
  //   VecDuplicate(_Hu_x, &_Hu_z);
  //   VecDuplicate(_rhs_x, &_rhs_z);
  //
  //   // Print if verbose
  //   if(_verbose_print)
  //   {
  //     std::cout << "Ainv_z: " << std::endl;
  //     VecView(_Ainv_z, PETSC_VIEWER_STDOUT_WORLD);
  //     std::cout << "Hu_z: " << std::endl;
  //     VecView(_Hu_z, PETSC_VIEWER_STDOUT_WORLD);
  //     std::cout << "rhs_z: " << std::endl;
  //     VecView(_rhs_z, PETSC_VIEWER_STDOUT_WORLD);
  //   }
  // }

  if(! _momentum_predictor_bool)
  {

    Vec _rhs_loc;
    MPI_Comm comm = MPI_COMM_WORLD;

    VecCreate(comm, &_rhs_loc);
    VecSetSizes(_rhs_loc, active_elements, PETSC_DECIDE);
    VecSetFromOptions(_rhs_loc);
    VecZeroEntries(_rhs_loc);

    std::unique_ptr<NumericVector<Number>> loc_zero_sol = isys.rhs->zero_clone();
    std::unique_ptr<NumericVector<Number>> loc_zero_rhs = isys.rhs->zero_clone();

    feProblem().computeResidualSys(isys, *loc_zero_sol.get(), *loc_zero_rhs.get());
    PetscVector<Number> * loc_prhs = dynamic_cast<PetscVector<Number> *>(loc_zero_sol.get());
    VecCopy(loc_prhs->vec(), _rhs_loc);

    VecScale(_rhs_loc, -1.0);

    if(_verbose_print)
    {
      std::cout << "RHS: " << std::endl;
      VecView(_rhs_loc, PETSC_VIEWER_STDOUT_WORLD);
    }

    VecDestroy(&_rhs_loc);
    //feProblem().getAuxiliarySystem().solution().close();

    std::cout << "Tutto bene." << std::endl;

  } else

  {
    // Declaration
    Vec _Ainv, _Hu, _rhs;
    Vec vec_dummy;
    Mat MC;
    MPI_Comm comm = MPI_COMM_WORLD;

    // Create vectors
    VecCreate(comm, &_Ainv);
    VecCreate(comm, &_Hu);
    VecCreate(comm, &_rhs);

    // Set Sizes
    VecSetSizes(_Ainv, active_elements*mesh_dimension, PETSC_DECIDE);
    VecSetSizes(_Hu, active_elements*mesh_dimension, PETSC_DECIDE);
    VecSetSizes(_rhs, active_elements*mesh_dimension, PETSC_DECIDE);

    // Set types
    VecSetFromOptions(_Ainv);
    VecSetFromOptions(_Hu);
    VecSetFromOptions(_rhs);

    // Assign zero to all entries
    VecZeroEntries(_Ainv);
    VecZeroEntries(_Hu);
    VecZeroEntries(_rhs);

    // Print if verbose
    // if(_verbose_print)
    // {
    //   std::cout << "Ainv: " << std::endl;
    //   VecView(_Ainv, PETSC_VIEWER_STDOUT_WORLD);
    //   std::cout << "Hu: " << std::endl;
    //   VecView(_Hu, PETSC_VIEWER_STDOUT_WORLD);
    //   std::cout << "rhs: " << std::endl;
    //   VecView(_rhs, PETSC_VIEWER_STDOUT_WORLD);
    // }


    // ** Asigning Ainv
    // Get Diagonal of A
    //PetscScalar loc_value;
    NumericVector<Number> * loc_residual = isys.rhs;
    PetscVector<Number> * ploc_residual = dynamic_cast<PetscVector<Number> *>(loc_residual);
    VecDuplicate(ploc_residual->vec(), &_Ainv);
    MatGetDiagonal(pmat->mat(), _Ainv);
    VecDuplicate(ploc_residual->vec(), &vec_dummy);
    VecSet(vec_dummy, 1.0);
    VecPointwiseDivide(_Ainv, vec_dummy, _Ainv);
    if(_verbose_print)
    {
      std::cout << "Ainv: " << std::endl;
      VecView(_Ainv, PETSC_VIEWER_STDOUT_WORLD);
    }

    // if(_verbose_print)
    // {
    //   std::cout << "Holder: " << std::endl;
    //   VecView(holder, PETSC_VIEWER_STDOUT_WORLD);
    // }
    //
    // PetscInt loc_dim = 0, loc_index = 0;
    // if(mesh_dimension > loc_dim)
    // {
    //   VecAssemblyBegin(_Ainv_x);
    //   for (PetscInt i = 0;
    //        i <= active_elements;
    //        i++)
    //   {
    //     loc_index = i*mesh_dimension;
    //     VecGetValues(holder, 1, &loc_index, &loc_value);
    //     VecSetValues(_Ainv_x, 1, &i, &loc_value, INSERT_VALUES);
    //   }
    //   VecAssemblyEnd(_Ainv_x);
    //   if(_verbose_print)
    //   {
    //     std::cout << "_Ainv_x: " << std::endl;
    //     VecView(_Ainv_x, PETSC_VIEWER_STDOUT_WORLD);
    //   }
    // }
    // loc_dim = 1;
    // if(mesh_dimension > loc_dim)
    // {
    //   VecAssemblyBegin(_Ainv_y);
    //   for (PetscInt i = 0;
    //        i <= active_elements;
    //        i++)
    //   {
    //     loc_index = i*mesh_dimension+1;
    //     VecGetValues(holder, 1, &loc_index, &loc_value);
    //     VecSetValues(_Ainv_y, 1, &i, &loc_value, INSERT_VALUES);
    //   }
    //   VecAssemblyEnd(_Ainv_y);
    //   if(_verbose_print)
    //   {
    //     std::cout << "_Ainv_y: " << std::endl;
    //     VecView(_Ainv_y, PETSC_VIEWER_STDOUT_WORLD);
    //   }
    // }
    // loc_dim = 2;
    // if(mesh_dimension > loc_dim)
    // {
    //   for (PetscInt i = 0;
    //        i <= active_elements;
    //        i++)
    //   {
    //     loc_index = i*mesh_dimension+2;
    //     VecGetValues(holder, 1, &loc_index, &loc_value);
    //     VecSetValues(_Ainv_z, 1, &i, &loc_value, INSERT_VALUES);
    //   }
    //   VecAssemblyEnd(_Ainv_z);
    //   if(_verbose_print)
    //   {
    //     std::cout << "_Ainv_z: " << std::endl;
    //     VecView(_Ainv_z, PETSC_VIEWER_STDOUT_WORLD);
    //   }
    // }

    // Creating HU
    MatDuplicate(pmat->mat(), MAT_COPY_VALUES, &MC);
    VecZeroEntries(vec_dummy);
    MatDiagonalSet(MC, vec_dummy, INSERT_VALUES);
    if(_verbose_print)
    {
      std::cout << "H matrix: " << std::endl;
      MatView(MC, PETSC_VIEWER_STDOUT_WORLD);
    }
    NumericVector<Number> * loc_solution = isys.solution.get();
    PetscVector<Number> * ploc_solution = dynamic_cast<PetscVector<Number> *>(loc_solution);
    // if(_verbose_print)
    // {
    //   std::cout << "RHS: " << std::endl;
    //   VecView(ploc_solution->vec(), PETSC_VIEWER_STDOUT_WORLD);
    // }
    VecDuplicate(vec_dummy, &_Hu);
    MatMult(MC, ploc_solution->vec(), _Hu);
    //VecPointwiseMult(_Hu, _Hu, _Ainv);
    //VecScale(_Hu, -1.0);
    if (_verbose_print)
    {
      std::cout << "_Hu: " << std::endl;
      VecView(_Hu, PETSC_VIEWER_STDOUT_WORLD);
    }

    // loc_dim = 0;
    // if(mesh_dimension > loc_dim)
    // {
    //   VecAssemblyBegin(_Hu_x);
    //   for (PetscInt i = 0;
    //        i <= active_elements;
    //        i++)
    //   {
    //     loc_index = i*mesh_dimension;
    //     VecGetValues(_Hu_global, 1, &loc_index, &loc_value);
    //     VecSetValues(_Hu_x, 1, &i, &loc_value, INSERT_VALUES);
    //   }
    //   VecAssemblyEnd(_Hu_x);
    //   if(_verbose_print)
    //   {
    //     std::cout << "_Hu_x: " << std::endl;
    //     VecView(_Hu_x, PETSC_VIEWER_STDOUT_WORLD);
    //   }
    // }
    // loc_dim = 1;
    // if(mesh_dimension > loc_dim)
    // {
    //   VecAssemblyBegin(_Hu_y);
    //   for (PetscInt i = 0;
    //        i <= active_elements;
    //        i++)
    //   {
    //     loc_index = i*mesh_dimension+1;
    //     VecGetValues(_Hu_global, 1, &loc_index, &loc_value);
    //     VecSetValues(_Hu_y, 1, &i, &loc_value, INSERT_VALUES);
    //   }
    //   VecAssemblyEnd(_Hu_y);
    //   if(_verbose_print)
    //   {
    //     std::cout << "_Hu_y: " << std::endl;
    //     VecView(_Hu_y, PETSC_VIEWER_STDOUT_WORLD);
    //   }
    // }
    // loc_dim = 2;
    // if(mesh_dimension > loc_dim)
    // {
    //   VecAssemblyBegin(_Hu_z);
    //   for (PetscInt i = 0;
    //        i <= active_elements;
    //        i++)
    //   {
    //     loc_index = i*mesh_dimension+2;
    //     VecGetValues(_Hu_global, 1, &loc_index, &loc_value);
    //     VecSetValues(_Hu_z, 1, &i, &loc_value, INSERT_VALUES);
    //   }
    //   VecAssemblyEnd(_Hu_z);
    //   if(_verbose_print)
    //   {
    //     std::cout << "_Hu_z: " << std::endl;
    //     VecView(_Hu_z, PETSC_VIEWER_STDOUT_WORLD);
    //   }
    // }

    // Getting RHS
    std::unique_ptr<NumericVector<Number>> zero_sol = isys.rhs->zero_clone();
    std::unique_ptr<NumericVector<Number>> zero_rhs = isys.rhs->zero_clone();
    feProblem().computeResidualSys(isys, *zero_sol.get(), *zero_rhs.get());
    PetscVector<Number> * prhs = dynamic_cast<PetscVector<Number> *>(zero_rhs.get());
    VecCopy(prhs->vec(), _rhs);
    VecScale(_rhs, -1.0);
    if(_verbose_print)
    {
      std::cout << "RHS: " << std::endl;
      VecView(_rhs, PETSC_VIEWER_STDOUT_WORLD);
    }

    // loc_dim = 0;
    // if(mesh_dimension > loc_dim)
    // {
    //   VecAssemblyBegin(_rhs_x);
    //   for (PetscInt i = 0;
    //        i <= active_elements;
    //        i++)
    //   {
    //     loc_index = i*mesh_dimension;
    //     VecGetValues(prhs->vec(), 1, &loc_index, &loc_value);
    //     VecSetValues(_rhs_x, 1, &i, &loc_value, INSERT_VALUES);
    //   }
    //   VecAssemblyEnd(_rhs_x);
    //   if(_verbose_print)
    //   {
    //     std::cout << "_rhs_x: " << std::endl;
    //     VecView(_rhs_x, PETSC_VIEWER_STDOUT_WORLD);
    //   }
    // }
    // loc_dim = 1;
    // if(mesh_dimension > loc_dim)
    // {
    //   VecAssemblyBegin(_rhs_y);
    //   for (PetscInt i = 0;
    //        i <= active_elements;
    //        i++)
    //   {
    //     loc_index = i*mesh_dimension+1;
    //     VecGetValues(prhs->vec(), 1, &loc_index, &loc_value);
    //     VecSetValues(_rhs_y, 1, &i, &loc_value, INSERT_VALUES);
    //   }
    //   VecAssemblyEnd(_rhs_y);
    //   if(_verbose_print)
    //   {
    //     std::cout << "_rhs_y: " << std::endl;
    //     VecView(_rhs_y, PETSC_VIEWER_STDOUT_WORLD);
    //   }
    // }
    // loc_dim = 2;
    // if(mesh_dimension > loc_dim)
    // {
    //   VecAssemblyBegin(_rhs_z);
    //   for (PetscInt i = 0;
    //        i <= active_elements;
    //        i++)
    //   {
    //     loc_index = i*mesh_dimension+2;
    //     VecGetValues(prhs->vec(), 1, &loc_index, &loc_value);
    //     VecSetValues(_rhs_z, 1, &i, &loc_value, INSERT_VALUES);
    //   }
    //   VecAssemblyEnd(_rhs_z);
    //   if(_verbose_print)
    //   {
    //     std::cout << "_rhs_z: " << std::endl;
    //     VecView(_rhs_z, PETSC_VIEWER_STDOUT_WORLD);
    //   }
    // }

    // Inserting variables into the auxiliary system
    AuxiliarySystem & aux_sys = feProblem().getAuxiliarySystem();

    // unsigned int pressure_number = aux_sys.system().variable_number("pressure");
    // auto pressure = dynamic_cast<const INSFVPressureVariable *>(aux_sys.system().variable(pressure_number));
    unsigned int loc_dim = 0;
    if(mesh_dimension > loc_dim)
    {
      unsigned int u_num = _nl.system().variable_number("u");
      unsigned int Ainv_x_num = aux_sys.system().variable_number("Ainv_x");
      unsigned int Hu_x_num = aux_sys.system().variable_number("Hu_x");
      unsigned int rhs_x_num = aux_sys.system().variable_number("RHS_x");
      auto begin = feProblem().mesh().getMesh().active_local_elements_begin();
      auto end = feProblem().mesh().getMesh().active_local_elements_end();
      PetscScalar val;

      for (auto it = begin; it != end; ++it)
      {
        const Elem * elem = *it;
        dof_id_type u_dof = elem->dof_number(_nl.number(), u_num, 0);
        dof_id_type Ainv_x_dof = elem->dof_number(aux_sys.number(), Ainv_x_num, 0);
        dof_id_type Hu_x_dof = elem->dof_number(aux_sys.number(), Hu_x_num, 0);
        dof_id_type rhs_x_dof = elem->dof_number(aux_sys.number(), rhs_x_num, 0);

        if (_verbose_print)
          std::cout << "xdofs: " << Ainv_x_dof << " " << Hu_x_dof << " " << rhs_x_dof << " " << u_dof << std::endl;

        PetscInt petsc_dof = (PetscInt) u_dof;

        VecGetValues(_Ainv, 1, &petsc_dof, &val);
        if (_verbose_print) std::cout << "Val Ainv_x: " << val << std::endl;
        aux_sys.solution().set(Ainv_x_dof, val);
        VecGetValues(_Hu, 1, &petsc_dof, &val);
        if (_verbose_print) std::cout << "Val Hu_x: " << val << std::endl;
        aux_sys.solution().set(Hu_x_dof, val);
        VecGetValues(_rhs, 1, &petsc_dof, &val);
        if (_verbose_print) std::cout << "Val RHS_x: " << val << std::endl;
        aux_sys.solution().set(rhs_x_dof, val);
      }
    }

    loc_dim = 1;
    if(mesh_dimension > loc_dim)
    {
      unsigned int v_num = _nl.system().variable_number("v");
      unsigned int Ainv_y_num = aux_sys.system().variable_number("Ainv_y");
      unsigned int Hu_y_num = aux_sys.system().variable_number("Hu_y");
      unsigned int rhs_y_num = aux_sys.system().variable_number("RHS_y");
      auto begin = feProblem().mesh().getMesh().active_local_elements_begin();
      auto end = feProblem().mesh().getMesh().active_local_elements_end();
      PetscScalar val;

      for (auto it = begin; it != end; ++it)
      {
        const Elem * elem = *it;
        dof_id_type v_dof = elem->dof_number(_nl.number(), v_num, 0);
        dof_id_type Ainv_y_dof = elem->dof_number(aux_sys.number(), Ainv_y_num, 0);
        dof_id_type Hu_y_dof = elem->dof_number(aux_sys.number(), Hu_y_num, 0);
        dof_id_type rhs_y_dof = elem->dof_number(aux_sys.number(), rhs_y_num, 0);

        if (_verbose_print)
          std::cout << "ydofs: " << Ainv_y_dof << " " << Hu_y_dof << " " << rhs_y_dof << " " << v_dof << std::endl;

        PetscInt petsc_dof = (PetscInt) v_dof;

        VecGetValues(_Ainv, 1, &petsc_dof, &val);
        if (_verbose_print) std::cout << "Val Ainv_y: " << val << std::endl;
        aux_sys.solution().set(Ainv_y_dof, val);
        VecGetValues(_Hu, 1, &petsc_dof, &val);
        if (_verbose_print) std::cout << "Val Hu_y: " << val << std::endl;
        aux_sys.solution().set(Hu_y_dof, val);
        VecGetValues(_rhs, 1, &petsc_dof, &val);
        if (_verbose_print) std::cout << "Val RHS_y: " << val << std::endl;
        aux_sys.solution().set(rhs_y_dof, val);
      }
    }

    loc_dim = 2;
    if(mesh_dimension > loc_dim)
    {
      unsigned int w_num = _nl.system().variable_number("w");
      unsigned int Ainv_z_num = aux_sys.system().variable_number("Ainv_z");
      unsigned int Hu_z_num = aux_sys.system().variable_number("Hu_z");
      unsigned int rhs_z_num = aux_sys.system().variable_number("RHS_z");
      auto begin = feProblem().mesh().getMesh().active_local_elements_begin();
      auto end = feProblem().mesh().getMesh().active_local_elements_end();
      PetscScalar val;

      for (auto it = begin; it != end; ++it)
      {
        const Elem * elem = *it;
        dof_id_type w_dof = elem->dof_number(_nl.number(), w_num, 0);
        dof_id_type Ainv_z_dof = elem->dof_number(aux_sys.number(), Ainv_z_num, 0);
        dof_id_type Hu_z_dof = elem->dof_number(aux_sys.number(), Hu_z_num, 0);
        dof_id_type rhs_z_dof = elem->dof_number(aux_sys.number(), rhs_z_num, 0);

        if (_verbose_print)
          std::cout << "zdofs: " << Ainv_z_dof << " " << Hu_z_dof << " " << rhs_z_dof << " " << w_dof << std::endl;

        PetscInt petsc_dof = (PetscInt) w_dof;

        VecGetValues(_Ainv, 1, &petsc_dof, &val);
        if (_verbose_print) std::cout << "Val Ainv_z: " << val << std::endl;
        aux_sys.solution().set(Ainv_z_dof, val);
        VecGetValues(_Hu, 1, &petsc_dof, &val);
        if (_verbose_print) std::cout << "Val Hu_z: " << val << std::endl;
        aux_sys.solution().set(Hu_z_dof, val);
        VecGetValues(_rhs, 1, &petsc_dof, &val);
        if (_verbose_print) std::cout << "Val RHS_z: " << val << std::endl;
        aux_sys.solution().set(rhs_z_dof, val);
      }
    }

    VecDestroy(&_Ainv);
    VecDestroy(&_Hu);
    VecDestroy(&_rhs);
    VecDestroy(&vec_dummy);
    MatDestroy(&MC);
    aux_sys.solution().close();
  }


    // insert Ainv into Aux
    // get Ainv_num

    // unsigned int v_num = _nl.system().variable_number("v");
    // unsigned int w_num = _nl.system().variable_number("w");
    // unsigned int Ainv_y_num = aux_sys.system().variable_number("Ainv_y");
    // unsigned int Ainv_z_num = aux_sys.system().variable_number("Ainv_z");
    // unsigned int Hu_y_num = aux_sys.system().variable_number("Hu_y");
    // unsigned int Hu_z_num = aux_sys.system().variable_number("Hu_z");
    //
    // //dof_id_type n_local_fluid_elems = 0;
    // auto begin = feProblem().mesh().getMesh().active_local_elements_begin();
    // auto end = feProblem().mesh().getMesh().active_local_elements_end();
    // for (auto it = begin; it != end; ++it)
    // {
    //   const Elem * elem = *it;
    //   dof_id_type vel_dof = elem->dof_number(_nl.number(), u_num, 0);
    //   dof_id_type Ainv_x_dof = elem->dof_number(aux_sys.number(), Ainv_num, 0);
    //
    //   // I am setting vel_dof to Ainv_dof because it's the same in this
    //   // stupid test but in general, just use vel_dof as computed above
    //   // the test I picked had a nodal var as nonlinear var
    //   vel_dof = Ainv_dof;
    //
    //   PetscScalar val;
    //   PetscInt petsc_dof = (PetscInt) vel_dof;
    //   VecGetValues(Ainv, 1, &petsc_dof, &val);
    //   aux_sys.solution().set(Ainv_dof, val);
    //   std::cout << Ainv_dof << " " << vel_dof << std::endl;
    // }
    // aux_sys.solution().close();

  // Getting the correct dimension
  // NumericVector<Number> * rhs = isys.rhs;
  // PetscVector<Number> * prhs = dynamic_cast<PetscVector<Number> *>(rhs);
  // if(_verbose_print)
  //   std::cout << "Residual: " << std::endl; VecView(prhs->vec(), PETSC_VIEWER_STDOUT_WORLD);


  // Creating zero vectors for allocation of Ainv and H*u
  // NumericVector<Number> * residual = isys.rhs;
  // PetscVector<Number> * presidual = dynamic_cast<PetscVector<Number> *>(residual);

  // NumericVector<Number> * rhs = isys.rhs;
  // PetscVector<Number> * prhs = dynamic_cast<PetscVector<Number> *>(rhs);
  // if(_verbose_print)
  //   std::cout << "Residual: " << std::endl; VecView(prhs->vec(), PETSC_VIEWER_STDOUT_WORLD);

  // NumericVector<Number> * sol = isys.solution.get();
  // PetscVector<Number> * psol = dynamic_cast<PetscVector<Number> *>(sol);
  // if(_verbose_print)
  //   std::cout << "RHS: " << std::endl; VecView(psol->vec(), PETSC_VIEWER_STDOUT_WORLD);
  //
  // // get the actual rhs
  // // std::unique_ptr<NumericVector<Number>> zero_sol = rhs->zero_clone();
  // // std::unique_ptr<NumericVector<Number>> zero_rhs = rhs->zero_clone();
  // // feProblem().computeResidualSys(isys, *zero_sol.get(), *zero_rhs.get());
  // //zero_sol->close();
  //
  // PetscVector<Number> * pzero_rhs = dynamic_cast<PetscVector<Number> *>(zero_rhs.get());
  // zero_rhs->close();
  // if (!pzero_rhs)
  //   mooseError("casting failed");
  //
  // Mat M;
  // MatDuplicate(pmat->mat(), MAT_COPY_VALUES, &M);
  // Vec Ainv;
  // VecDuplicate(prhs->vec(), &Ainv);
  // Vec AinvHu;
  // VecDuplicate(prhs->vec(), &AinvHu);
  // Vec zerovec;
  // VecDuplicate(prhs->vec(), &zerovec);
  // VecSet(zerovec, 0.0);
  // Vec onevec;
  // VecDuplicate(prhs->vec(), &onevec);
  // VecSet(onevec, 1.0);
  // Vec monevec;
  // VecDuplicate(prhs->vec(), &monevec);
  // VecSet(monevec, -1.0);
  // MatGetDiagonal(M, Ainv);
  // VecPointwiseDivide(Ainv, onevec, Ainv);
  //
  // // compute H, note we multiply Ainv to it already
  // MatDiagonalSet(M, zerovec, INSERT_VALUES);
  // MatMult(M, psol->vec(), AinvHu);
  // VecPointwiseMult(AinvHu, monevec, AinvHu);
  // // NOTE: we get the negative rhs
  // VecAXPY(AinvHu, -1.0, pzero_rhs->vec());
  // VecPointwiseMult(AinvHu, Ainv, AinvHu);
  // //VecView(pzero_rhs->vec(), PETSC_VIEWER_STDOUT_WORLD);
  //
  // // get aux system
  // AuxiliarySystem & aux_sys = feProblem().getAuxiliarySystem();
  //
  // // insert Ainv into Aux
  // // get Ainv_num
  // unsigned int Ainv_num = aux_sys.system().variable_number("Ainv");
  // unsigned int vel_num = _nl.system().variable_number("u");
  //
  // //dof_id_type n_local_fluid_elems = 0;
  // auto begin = feProblem().mesh().getMesh().active_local_elements_begin();
  // auto end = feProblem().mesh().getMesh().active_local_elements_end();
  // for (auto it = begin; it != end; ++it)
  // {
  //   const Elem * elem = *it;
  //   dof_id_type Ainv_dof = elem->dof_number(aux_sys.number(), Ainv_num, 0);
  //   dof_id_type vel_dof = elem->dof_number(_nl.number(), vel_num, 0);
  //
  //   // I am setting vel_dof to Ainv_dof because it's the same in this
  //   // stupid test but in general, just use vel_dof as computed above
  //   // the test I picked had a nodal var as nonlinear var
  //   vel_dof = Ainv_dof;
  //
  //   PetscScalar val;
  //   PetscInt petsc_dof = (PetscInt) vel_dof;
  //   VecGetValues(Ainv, 1, &petsc_dof, &val);
  //   aux_sys.solution().set(Ainv_dof, val);
  //   std::cout << Ainv_dof << " " << vel_dof << std::endl;
  // }
  // aux_sys.solution().close();
}

void
CustomTransient::execute()
{
  preExecute();

  // Start time loop...
  while (keepGoing())
  {
    incrementStepOrReject();
    preStep();
    computeDT();
    takeStep();
    endStep();
    postStep();
  }

  if (lastSolveConverged())
  {
    _t_step++;

    /*
     * Call the multi-app executioners endStep and
     * postStep methods when doing Picard or when not automatically advancing sub-applications for
     * some other reason. We do not perform these calls for loose-coupling/auto-advancement
     * problems because Transient::endStep and Transient::postStep get called from
     * TransientMultiApp::solveStep in that case.
     */
    if (!_fixed_point_solve->autoAdvance())
    {
      _problem.finishMultiAppStep(EXEC_TIMESTEP_BEGIN, /*recurse_through_multiapp_levels=*/true);
      _problem.finishMultiAppStep(EXEC_TIMESTEP_END, /*recurse_through_multiapp_levels=*/true);
    }
  }

  if (!_app.halfTransient())
  {
    TIME_SECTION("final", 1, "Executing Final Objects");
    _problem.execMultiApps(EXEC_FINAL);
    _problem.finalizeMultiApps();
    _problem.execute(EXEC_FINAL);
    _problem.outputStep(EXEC_FINAL);
  }

  // This method is to finalize anything else we want to do on the problem side.
  _problem.postExecute();

  // This method can be overridden for user defined activities in the Executioner.
  postExecute();
}

void
CustomTransient::computeDT()
{
  _time_stepper->computeStep(); // This is actually when DT gets computed
}

void
CustomTransient::incrementStepOrReject()
{
  if (lastSolveConverged())
  {
    if (!_xfem_repeat_step)
    {
#ifdef LIBMESH_ENABLE_AMR
      if (_t_step != 0)
        _problem.adaptMesh();
#endif

      _time_old = _time;
      _t_step++;

      _problem.advanceState();

      if (_t_step == 1)
        return;

      /*
       * Call the multi-app executioners endStep and
       * postStep methods when doing Picard or when not automatically advancing sub-applications for
       * some other reason. We do not perform these calls for loose-coupling/auto-advancement
       * problems because Transient::endStep and Transient::postStep get called from
       * TransientMultiApp::solveStep in that case.
       */
      if (!_fixed_point_solve->autoAdvance())
      {
        _problem.finishMultiAppStep(EXEC_TIMESTEP_BEGIN);
        _problem.finishMultiAppStep(EXEC_TIMESTEP_END);
      }

      /*
       * Ensure that we increment the sub-application time steps so that
       * when dt selection is made in the master application, we are using
       * the correct time step information
       */
      _problem.incrementMultiAppTStep(EXEC_TIMESTEP_BEGIN);
      _problem.incrementMultiAppTStep(EXEC_TIMESTEP_END);
    }
  }
  else
  {
    _problem.restoreMultiApps(EXEC_TIMESTEP_BEGIN, true);
    _problem.restoreMultiApps(EXEC_TIMESTEP_END, true);
    _time_stepper->rejectStep();
    _time = _time_old;
  }
}

void
CustomTransient::takeStep(Real input_dt)
{
  _dt_old = _dt;

  if (input_dt == -1.0)
    _dt = computeConstrainedDT();
  else
    _dt = input_dt;

  _time_stepper->preSolve();

  // Increment time
  _time = _time_old + _dt;

  _problem.timestepSetup();

  _problem.onTimestepBegin();

  _time_stepper->step();
  _xfem_repeat_step = _fixed_point_solve->XFEMRepeatStep();

  _last_solve_converged_custom = _time_stepper->converged();

  if (!lastSolveConverged())
  {
    _console << "Aborting as solve did not converge" << std::endl;
    return;
  }

  if (!(_problem.haveXFEM() && _fixed_point_solve->XFEMRepeatStep()))
  {
    if (lastSolveConverged())
      _time_stepper->acceptStep();
    else
      _time_stepper->rejectStep();
  }

  _time = _time_old;

  _time_stepper->postSolve();

  _solution_change_norm_custom =
      relativeSolutionDifferenceNorm() / (_normalize_solution_diff_norm_by_dt ? _dt : Real(1));

  return;
}

void
CustomTransient::endStep(Real input_time)
{
  if (input_time == -1.0)
    _time = _time_old + _dt;
  else
    _time = input_time;

  if (lastSolveConverged())
  {
    if (_xfem_repeat_step)
      _time = _time_old;
    else
    {
      _nl.getTimeIntegrator()->postStep();

      // Compute the Error Indicators and Markers
      _problem.computeIndicators();
      _problem.computeMarkers();

      // Perform the output of the current time step
      _problem.outputStep(EXEC_TIMESTEP_END);

      // output
      if (_time_interval_custom && (_time + _timestep_tolerance >= _next_interval_output_time))
        _next_interval_output_time += _time_interval_custom_output_interval;
    }
  }
}

Real
CustomTransient::computeConstrainedDT()
{
  //  // If start up steps are needed
  //  if (_t_step == 1 && _n_startup_steps > 1)
  //    _dt = _input_dt/(double)(_n_startup_steps);
  //  else if (_t_step == 1+_n_startup_steps && _n_startup_steps > 1)
  //    _dt = _input_dt;

  Real dt_cur = _dt;
  std::ostringstream diag;

  // After startup steps, compute new dt
  if (_t_step > _n_startup_steps)
    dt_cur = getDT();

  else
  {
    diag << "Timestep < n_startup_steps, using old dt: " << std::setw(9) << std::setprecision(6)
         << std::setfill('0') << std::showpoint << std::left << _dt << " tstep: " << _t_step
         << " n_startup_steps: " << _n_startup_steps << std::endl;
  }
  _unconstrained_dt_custom = dt_cur;

  if (_verbose)
    _console << diag.str();

  diag.str("");
  diag.clear();

  // Allow the time stepper to limit the time step
  _at_sync_point_custom = _time_stepper->constrainStep(dt_cur);

  // Don't let time go beyond next time interval output if specified
  if ((_time_interval_custom) && (_time + dt_cur + _timestep_tolerance >= _next_interval_output_time))
  {
    dt_cur = _next_interval_output_time - _time;
    _at_sync_point_custom = true;

    diag << "Limiting dt for time interval output at time: " << std::setw(9) << std::setprecision(6)
         << std::setfill('0') << std::showpoint << std::left << _next_interval_output_time
         << " dt: " << std::setw(9) << std::setprecision(6) << std::setfill('0') << std::showpoint
         << std::left << dt_cur << std::endl;
  }

  // If a target time is set and the current dt would exceed it, limit dt to match the target
  if (_target_time_custom > -std::numeric_limits<Real>::max() + _timestep_tolerance &&
      _time + dt_cur + _timestep_tolerance >= _target_time_custom)
  {
    dt_cur = _target_time_custom - _time;
    _at_sync_point_custom = true;

    diag << "Limiting dt for target time: " << std::setw(9) << std::setprecision(6)
         << std::setfill('0') << std::showpoint << std::left << _next_interval_output_time
         << " dt: " << std::setw(9) << std::setprecision(6) << std::setfill('0') << std::showpoint
         << std::left << dt_cur << std::endl;
  }

  // Constrain by what the multi apps are doing
  Real multi_app_dt = _problem.computeMultiAppsDT(EXEC_TIMESTEP_BEGIN);
  if (_use_multiapp_dt || multi_app_dt < dt_cur)
  {
    dt_cur = multi_app_dt;
    _at_sync_point_custom = false;
    diag << "Limiting dt for MultiApps: " << std::setw(9) << std::setprecision(6)
         << std::setfill('0') << std::showpoint << std::left << dt_cur << std::endl;
  }
  multi_app_dt = _problem.computeMultiAppsDT(EXEC_TIMESTEP_END);
  if (multi_app_dt < dt_cur)
  {
    dt_cur = multi_app_dt;
    _at_sync_point_custom = false;
    diag << "Limiting dt for MultiApps: " << std::setw(9) << std::setprecision(6)
         << std::setfill('0') << std::showpoint << std::left << dt_cur << std::endl;
  }

  if (_verbose)
    _console << diag.str();

  return dt_cur;
}

Real
CustomTransient::getDT()
{
  return _time_stepper->getCurrentDT();
}

bool
CustomTransient::keepGoing()
{
  bool keep_going = !_problem.isSolveTerminationRequested();

  // Check for stop condition based upon steady-state check flag:
  if (lastSolveConverged())
  {
    if (!_xfem_repeat_step)
    {
      if (_steady_state_detection == true && _time > _steady_state_start_time)
      {
        // Check solution difference relative norm against steady-state tolerance
        if (_solution_change_norm_custom < _steady_state_tolerance)
        {
          _console << "Steady-State Solution Achieved at time: " << _time << std::endl;
          // Output last solve if not output previously by forcing it
          keep_going = false;
        }
        else // Keep going
        {
          // Print steady-state relative error norm
          _console << "Steady-State Relative Differential Norm: " << _solution_change_norm_custom
                   << std::endl;
        }
      }

      // Check for stop condition based upon number of simulation steps and/or solution end time:
      if (static_cast<unsigned int>(_t_step) >= _num_steps)
        keep_going = false;

      if ((_time >= _end_time) || (fabs(_time - _end_time) <= _timestep_tolerance))
        keep_going = false;
    }
  }
  else if (_abort)
  {
    _console << "Aborting as solve did not converge and input selected to abort" << std::endl;
    keep_going = false;
  }
  else if (!_error_on_dtmin && _dt <= _dtmin)
  {
    _console << "Aborting as timestep already at or below dtmin" << std::endl;
    keep_going = false;
  }

  return keep_going;
}

void
CustomTransient::estimateTimeError()
{
}

bool
CustomTransient::lastSolveConverged() const
{
  return _last_solve_converged_custom;
}

void
CustomTransient::postExecute()
{
  _time_stepper->postExecute();
}

void
CustomTransient::setTargetTime(Real target_time_custom)
{
  _target_time_custom = target_time_custom;
}

Real
CustomTransient::getSolutionChangeNorm()
{
  return _solution_change_norm_custom;
}

void
CustomTransient::setupTimeIntegrator()
{
  if (_pars.isParamSetByUser("scheme") && _problem.hasTimeIntegrator())
    mooseError("You cannot specify time_scheme in the Executioner and independently add a "
               "TimeIntegrator to the system at the same time");

  if (!_problem.hasTimeIntegrator())
  {
    // backwards compatibility
    std::string ti_str;
    using namespace Moose;

    switch (_time_scheme)
    {
      case TI_IMPLICIT_EULER:
        ti_str = "ImplicitEuler";
        break;
      case TI_EXPLICIT_EULER:
        ti_str = "ExplicitEuler";
        break;
      case TI_CRANK_NICOLSON:
        ti_str = "CrankNicolson";
        break;
      case TI_BDF2:
        ti_str = "BDF2";
        break;
      case TI_EXPLICIT_MIDPOINT:
        ti_str = "ExplicitMidpoint";
        break;
      case TI_LSTABLE_DIRK2:
        ti_str = "LStableDirk2";
        break;
      case TI_EXPLICIT_TVD_RK_2:
        ti_str = "ExplicitTVDRK2";
        break;
      case TI_NEWMARK_BETA:
        ti_str = "NewmarkBeta";
        break;
      default:
        mooseError("Unknown scheme: ", _time_scheme);
        break;
    }

    InputParameters params = _app.getFactory().getValidParams(ti_str);
    _problem.addTimeIntegrator(ti_str, ti_str, params);
  }
}

std::string
CustomTransient::getTimeStepperName()
{
  if (_time_stepper)
  {
    TimeStepper & ts = *_time_stepper;
    return demangle(typeid(ts).name());
  }
  else
    return std::string();
}

Real
CustomTransient::relativeSolutionDifferenceNorm()
{
  const NumericVector<Number> & current_solution = *_nl.currentSolution();
  const NumericVector<Number> & old_solution = _nl.solutionOld();

  _sln_diff = current_solution;
  _sln_diff -= old_solution;

  return (_sln_diff.l2_norm() / current_solution.l2_norm());
}
