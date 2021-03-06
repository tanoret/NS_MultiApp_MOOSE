//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVElementalKernel.h"

class FVFunctorDivergence : public FVElementalKernel
{
public:
  static InputParameters validParams();

  FVFunctorDivergence(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  const Real _sign;
  const Moose::Functor<ADReal> & _x;
  const Moose::Functor<ADReal> & _y;
  const Moose::Functor<ADReal> & _z;
};
