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
 * This class computes the pressure Poisson solve which is part of
 * the "split" scheme used for solving the incompressible Navier-Stokes
 * equations.
 */
class NavStokesPressurePoisson_p: public Kernel
{
public:
  static InputParameters validParams();

  NavStokesPressurePoisson_p(const InputParameters & parameters);

  virtual ~NavStokesPressurePoisson_p() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Gradients of the "star" velocity
  const VariableGradient & _grad_u_star;
  const VariableGradient & _grad_v_star;
  const VariableGradient & _grad_w_star;

  // Pressure Gradient old
  const VariableGradient & _grad_p_old;

  // Variable numberings
  unsigned _u_vel_star_var_number;
  unsigned _v_vel_star_var_number;
  unsigned _w_vel_star_var_number;

  // Material properties
  const MaterialProperty<Real> & _rho;
};
