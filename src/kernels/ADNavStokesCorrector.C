//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADNavStokesCorrector.h"
#include "MooseMesh.h"

registerMooseObject("AirfoilAppApp", ADNavStokesCorrector);

InputParameters
ADNavStokesCorrector::validParams()
{
  InputParameters params = ADKernel::validParams();

  params.addClassDescription("This class computes the 'Chorin' Corrector equation in "
                             "fully-discrete (both time and space) form.");
  // Coupled variables
  params.addRequiredCoupledVar("u", "star velocity vector");
  params.addRequiredCoupledVar("p", "pressure");

  params.addParam<MaterialPropertyName>("rho_name", "rho", "density name");

  return params;
}

ADNavStokesCorrector::ADNavStokesCorrector(const InputParameters & parameters)
  : ADVectorKernel(parameters),

    // Current velocities
    _velocity(adCoupledVectorValue("u")),

    // Pressure gradient
    _grad_p(adCoupledGradient("p")),

    // Material properties
    _rho(getADMaterialProperty<Real>("rho_name"))
{
}

ADReal
ADNavStokesCorrector::computeQpResidual()
{
  // The symmetric part
  ADReal symmetric_part = (_u[_qp] - _velocity[_qp]) * _test[_i][_qp];

  // The pressure part, don't forget to multiply by dt!
  ADReal pressure_part = (_dt / _rho[_qp]) * _grad_p[_qp] * _test[_i][_qp];

  return symmetric_part + pressure_part;
}
