// This file is part of the MOOSE framework
// https://www.mooseframework.org
//
// All rights reserved, see COPYRIGHT for full restrictions
// https://github.com/idaholab/moose/blob/master/COPYRIGHT
//
// Licensed under LGPL 2.1, please see LICENSE for details
// https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * Auxiliary kernel responsible for computing the Darcy velocity given
 * several fluid properties and the pressure gradient.
 */
class Corrector : public AuxKernel
{
public:
  static InputParameters validParams();

  Corrector(const InputParameters & parameters);

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue() override;

  // "Star" velocity components
  const VariableValue & _u_vel_star;
  const VariableValue & _v_vel_star;
  const VariableValue & _w_vel_star;

  // Pressure gradients
  const VariableGradient & _grad_p;
  // Pressure Gradient old
  const VariableGradient & _grad_p_old;

  // Parameters
  unsigned _component;

  // Material properties
  const MaterialProperty<Real> & _rho;
};
