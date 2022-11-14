//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVCopyKernel.h"
#include "MooseMesh.h"

registerMooseObject("AirfoilAppApp", FVCopyKernel);

InputParameters
FVCopyKernel::validParams()
{
  auto params = AuxKernel::validParams();

  params.addClassDescription("Object to copy all kernel.");
  params.addRequiredCoupledVar("FVVar", "FVVar from momenutm predictor.");

  return params;
}

FVCopyKernel::FVCopyKernel(const InputParameters & parameters)
  :
    // Base class
    AuxKernel(parameters),

    // Get coupled predictor variables
    _FVVar(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("FVVar", 0)))
{
}

Real
FVCopyKernel::computeValue()
{
  /// Diffusion term
  using namespace Moose::FV;

  return _FVVar->getElemValue(_current_elem).value();
}
