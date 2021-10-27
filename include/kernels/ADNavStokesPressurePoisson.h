//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernel.h"

// Forward Declarations

/**
 * This class computes the pressure Poisson solve which is part of
 * the "split" scheme used for solving the incompressible Navier-Stokes
 * equations.
 */
class ADNavStokesPressurePoisson : public ADKernel
{
public:
  static InputParameters validParams();

  ADNavStokesPressurePoisson(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  /// gradient of velocity
  const ADVectorVariableGradient & _grad_velocity;

  const ADMaterialProperty<Real> & _rho;

};
