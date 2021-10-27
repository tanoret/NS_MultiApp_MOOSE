//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// Forward Declarations
class NormalizationAuxOld;

template <>
InputParameters validParams<NormalizationAuxOld>();

/**
 * This auxiliary kernel normalizes a variable based on a postprocessor.
 * Typically this postprocessor is a norm of the variable to be normalized.
 * The option to shift the value is also provided.
 */
class NormalizationAuxOld : public AuxKernel
{
public:
  static InputParameters validParams();

  NormalizationAuxOld(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _src;
  const Real * const _pp_on_source;
  const Real * const _shift;
  Real _normal_factor;
};
