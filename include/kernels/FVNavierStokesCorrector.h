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

#include "INSFVVelocityVariable.h"
#include "INSFVPressureVariable.h"

class FVNavierStokesCorrector : public FVFluxKernel
{
public:
  static InputParameters validParams();

  FVNavierStokesCorrector(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;
  // void computeResidual() override;
  // void computeJacobian() override;
  // void computeOffDiagJacobian() override;

  // Pressures
  const INSFVPressureVariable * const _p_var;
  const INSFVPressureVariable * const _p_old;

  /// Transfer Variables
  const MooseVariableFVReal * const _Ainv;
  const MooseVariableFVReal * const _Hhat;

  /// Access to current direction
  const unsigned int _index;
  const Real _advection_relaxation;

};
