//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

// Forward Declarations

/**
 * This class computes the "Chorin" Predictor equation in fully-discrete
 * (both time and space) form.
 */
class NavStokesPredictor_p : public Kernel
{
public:
  static InputParameters validParams();

  NavStokesPredictor_p(const InputParameters & parameters);

  virtual ~NavStokesPredictor_p() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Old Velocity
  const VariableValue & _u_vel_old;
  const VariableValue & _v_vel_old;
  const VariableValue & _w_vel_old;

  // Star Velocity
  const VariableValue & _u_vel_star;
  const VariableValue & _v_vel_star;
  const VariableValue & _w_vel_star;

  // Star Velocity Gradients
  const VariableGradient & _grad_u_vel_star;
  const VariableGradient & _grad_v_vel_star;
  const VariableGradient & _grad_w_vel_star;

  // Pressure Gradient old
  const VariableGradient & _grad_p_old;

  // Star velocity numbers
  unsigned _u_vel_star_var_number;
  unsigned _v_vel_star_var_number;
  unsigned _w_vel_star_var_number;

  // Parameters
  unsigned _component;

  // Material properties
  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _rho;
};
