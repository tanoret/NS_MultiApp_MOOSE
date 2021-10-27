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
 * This class computes the "Chorin" Corrector equation in fully-discrete
 * (both time and space) form.
 */
class ADNavStokesCorrector : public ADVectorKernel
{
public:
  static InputParameters validParams();

  ADNavStokesCorrector(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  //velocity
  const ADVectorVariableValue & _velocity;

  // Pressure gradients
  const ADVariableGradient & _grad_p;

  // Material properties
  const ADMaterialProperty<Real> & _rho;
};
