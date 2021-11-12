//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

//#include "Executioner.h"
#include "Transient.h"

// System includes
#include <string>
#include <fstream>

// Forward Declarations
class CustomTransient;
class TimeStepper;
class FEProblemBase;

template <>
InputParameters validParams<CustomTransient>();

/**
 * CustomTransient executioners usually loop through a number of timesteps... calling solve()
 * for each timestep.
 */
class CustomTransient : public Transient
{
public:
  /**
   * Constructor
   *
   * @param parameters The parameters object holding data for the class to use.
   * @return Whether or not the solve was successful.
   */
  static InputParameters validParams();

  CustomTransient(const InputParameters & parameters);

  virtual void init() override;

  virtual void execute() override;

  /**
   * Do whatever is necessary to advance one step.
   */
  virtual void takeStep(Real input_dt = -1.0) override;

  /**
   * @return The fully constrained dt for this timestep
   */
  virtual Real computeConstrainedDT() override;
  virtual void estimateTimeError() override;

  /**
   * @return The the computed dt to use for this timestep.
   */
  virtual Real getDT() override;

  /**
   * CustomTransient loop will continue as long as this keeps returning true.
   */
  virtual bool keepGoing() override;

  /**
   * Whether or not the last solve converged.
   */
  virtual bool lastSolveConverged() const override;

  virtual void preExecute() override;

  virtual void postExecute() override;

  virtual void computeDT() override;

  virtual void preStep() override;

  virtual void postStep() override;

  /**
   * This is where the solve step is actually incremented.
   */
  virtual void incrementStepOrReject() override;

  virtual void endStep(Real input_time = -1.0) override;

  /**
   * Can be used to set the next "target time" which is a time to nail perfectly.
   * Useful for driving MultiApps.
   */
  virtual void setTargetTime(Real target_time_custom) override;

  /**
   * Get the current time.
   */
  virtual Real getTime() override { return _time; };

  /**
   * Get the current target time
   * @return target time
   */
  virtual Real getTargetTime() override { return _target_time_custom; }

  /**
   * Set the current time.
   */
  virtual void setTime(Real t) override { _time = t; };

  /**
   * Set the old time.
   */
  virtual void setTimeOld(Real t) override { _time_old = t; };

  /**
   * Get the Relative L2 norm of the change in the solution.
   */
  Real getSolutionChangeNorm();

  /**
   * Pointer to the TimeStepper
   * @return Pointer to the time stepper for this Executioner
   */
  TimeStepper * getTimeStepper() { return _time_stepper.get(); }

  /**
   * Set the timestepper to use.
   *
   * @param ts The TimeStepper to use
   */
  void setTimeStepper(std::shared_ptr<TimeStepper> ts) { _time_stepper = ts; }

  /**
   * Get the timestepper.
   */
  virtual std::string getTimeStepperName() override;

  /**
   * Get the time scheme used
   * @return MooseEnum with the time scheme
   */
  Moose::TimeIntegratorType getTimeScheme() { return _time_scheme; }

  /**
   * Get the set of sync times
   * @return The reference to the set of sync times
   */
  std::set<Real> & syncTimes() { return _sync_times; }

  /**
   * Get the maximum dt
   * @return The maximum dt
   */
  Real & dtMax() { return _dtmax; }

  /**
   * Get the minimal dt
   * @return The minimal dt
   */
  Real & dtMin() { return _dtmin; }

  /**
   * Return the start time
   * @return The start time
   */
  Real getStartTime() { return _start_time; }

  /**
   * Get the end time
   * @return The end time
   */
  Real & endTime() { return _end_time; }

  /**
   * Get the timestep tolerance
   * @return The timestep tolerance
   */
  Real & timestepTol() { return _timestep_tolerance; }

  /**
   * Set the timestep tolerance
   * @param tolerance timestep tolerance
   */
  virtual void setTimestepTolerance(const Real & tolerance) override { _timestep_tolerance = tolerance; }

  /**
   * Is the current step at a sync point (sync times, time interval, target time, etc)?
   * @return Bool indicataing whether we are at a sync point
   */
  bool atSyncPoint() { return _at_sync_point_custom; }

  /**
   * Get the unconstrained dt
   * @return Value of dt before constraints were applied
   */
  Real unconstrainedDT() { return _unconstrained_dt_custom; }

  void parentOutputPositionChanged() override { _fe_problem.parentOutputPositionChanged(); }

  /**
   * The relative L2 norm of the difference between solution and old solution vector.
   */
  virtual Real relativeSolutionDifferenceNorm() override;

  /**
   * Set the number of time steps
   * @param num_steps number of time steps
   */
  virtual void forceNumSteps(const unsigned int num_steps) override { _num_steps = num_steps; }

  /// Return the solve object wrapped by time stepper
  virtual SolveObject * timeStepSolveObject() override { return _fixed_point_solve.get(); }

protected:
  /// Here for backward compatibility
  FEProblemBase & _problem;

  /// inner-most solve object to perform Newton solve with PETSc on every time step
  FEProblemSolve _feproblem_solve;

  /// Reference to nonlinear system base for faster access
  NonlinearSystemBase & _nl;

  Moose::TimeIntegratorType _time_scheme;
  std::shared_ptr<TimeStepper> _time_stepper;

  /// Current timestep.
  int & _t_step;
  /// Current time
  Real & _time;
  /// Previous time
  Real & _time_old;
  /// Current delta t... or timestep size.
  Real & _dt;
  Real & _dt_old;

  Real & _unconstrained_dt_custom;
  bool & _at_sync_point_custom;

  /// Whether or not the last solve converged
  bool & _last_solve_converged_custom;

  /// Whether step should be repeated due to xfem modifying the mesh
  bool _xfem_repeat_step;

  Real _end_time;
  Real _dtmin;
  Real _dtmax;
  unsigned int _num_steps;
  int _n_startup_steps;

  /**
   * Steady state detection variables:
   */
  bool _steady_state_detection;
  Real _steady_state_tolerance;
  Real _steady_state_start_time;

  std::set<Real> & _sync_times;

  bool _abort;
  /// This parameter controls how the system will deal with _dt <= _dtmin
  /// If true, the time stepper is expected to throw an error
  /// If false, the executioner will continue through EXEC_FINAL
  const bool _error_on_dtmin;

  ///if to use time interval output
  bool & _time_interval_custom;
  Real _next_interval_output_time;
  Real _time_interval_custom_output_interval;

  Real _start_time;
  Real _timestep_tolerance;
  Real & _target_time_custom;
  bool _use_multiapp_dt;

  Real & _solution_change_norm_custom;

  /// The difference of current and old solutions
  NumericVector<Number> & _sln_diff;

  void setupTimeIntegrator();

  /// Whether to divide the solution difference norm by dt. If taking 'small' time steps this member
  /// should probably be true. If taking very 'large' timesteps in an attempt to reach a
  /// steady-state, this member should probably be be false.
  const bool _normalize_solution_diff_norm_by_dt;
  const bool & _verbose_print;
  const bool & _momentum_predictor_bool;
};
