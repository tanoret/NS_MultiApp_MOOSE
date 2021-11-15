//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVFluxKernel.h"

// Forward Declarations

/**
 * This class computes the pressure Poisson solve which is part of
 * the "split" scheme used for solving the incompressible Navier-Stokes
 * equations.
 */
class FVNavStokesPressurePredictor_p: public FVFluxKernel
{
public:
  static InputParameters validParams();

  FVNavStokesPressurePredictor_p(const InputParameters & parameters);

  virtual ~FVNavStokesPressurePredictor_p() {}

protected:

  // Override QP residual
  virtual ADReal computeQpResidual() override;

  // Thermo-physical properties
  // const Moose::Functor<ADReal> & _rho;
  // const Moose::Functor<ADReal> & _mu;

  /// pressure variable
  //const INSFVPressureVariable * const _p_var;
  /// pressure variable old
  // const INSFVPressureVariable * const _p_old;
  // /// x-velocity
  // const INSFVVelocityVariable * const _u_star;
  // /// y-velocity
  // const INSFVVelocityVariable * const _v_star;
  // /// z-velocity
  // const INSFVVelocityVariable * const _w_star;

  /// Transfer Variables
  const MooseVariableFVReal * const _Ainv_x;
  const MooseVariableFVReal * const _Ainv_y;
  const MooseVariableFVReal * const _Ainv_z;

  const MooseVariableFVReal * const _Hu_x;
  const MooseVariableFVReal * const _Hu_y;
  const MooseVariableFVReal * const _Hu_z;

  // const Moose::Functor<ADReal> & _Hu_x;
  // const Moose::Functor<ADReal> & _Hu_y;
  // const Moose::Functor<ADReal> & _Hu_z;

  // const MooseVariableFVReal * const _rhs_x;
  // const MooseVariableFVReal * const _rhs_y;
  // const MooseVariableFVReal * const _rhs_z;

  // Variable numberings
  // unsigned _u_vel_star_var_number;
  // unsigned _v_vel_star_var_number;
  // unsigned _w_vel_star_var_number;

  //// Access to the current element
  const Elem * const & _current_elem;
  unsigned int counter;

};
