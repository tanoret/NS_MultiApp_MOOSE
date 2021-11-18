//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "pressureRelaxation.h"
#include "MooseMesh.h"

registerMooseObject("AirfoilAppApp", pressureRelaxation);

InputParameters
pressureRelaxation::validParams()
{
  auto params = AuxKernel::validParams();

  params.addClassDescription("Object pressure predictor equation.");
  // Coupled variables
  params.addRequiredCoupledVar("pressure", "The pressure variable.");
  params.addRequiredCoupledVar("pressure_old", "The pressure variable.");

  params.addRangeCheckedParam<Real>("pressure_relaxation", 1.0,
                                    "(0 <= pressure_relaxation) & (pressure_relaxation <= 1)",
                                    "Pressure field relaxation must be between 0 and 1 both included.");

  return params;
}

pressureRelaxation::pressureRelaxation(const InputParameters & parameters)
  :
    // Base class
    AuxKernel(parameters),

    // Get coupled variables
    _p_var(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure", 0))),
    _p_old(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure_old", 0))),

    // Get pressure relaxation factor
    _pressure_relaxation(getParam<Real>("pressure_relaxation"))
{
}

Real
pressureRelaxation::computeValue()
{
  auto new_pressure = _p_old->getElemValue(_current_elem) +
                      _pressure_relaxation * _p_var->getElemValue(_current_elem);

  return new_pressure.value();
}
