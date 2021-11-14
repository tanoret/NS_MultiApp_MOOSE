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

#include "INSFVVelocityVariable.h"
#include "INSFVPressureVariable.h"

/**
 * Auxiliary kernel responsible for computing the Darcy velocity given
 * several fluid properties and the pressure gradient.
 */
class FVCorrector : public AuxKernel
{
public:
  static InputParameters validParams();

  FVCorrector(const InputParameters & parameters);

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue() override;

  // Pressures
  const INSFVPressureVariable * const _p_var;
  const INSFVPressureVariable * const _p_old;

  // Velocity components
  //const INSFVVelocityVariable * const _vel_star;

  /// Transfer Variables
  const MooseVariableFVReal * const _Ainv;
  const MooseVariableFVReal * const _Hhat;

  //// Access to the current element
  //const Elem * const & _current_elem;

  /// Access to current direction
  const unsigned int _index;
  const Real _pressure_relaxation;
  const VariableValue & _u_old;

};
