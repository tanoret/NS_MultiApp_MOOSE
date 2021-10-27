//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADNavStokesPressurePoisson.h"
#include "MooseMesh.h"

registerMooseObject("AirfoilAppApp", ADNavStokesPressurePoisson);

InputParameters
ADNavStokesPressurePoisson::validParams()
{
  InputParameters params = ADKernel::validParams();

  params.addClassDescription("This class computes the pressure Poisson solve which is part of the "
                             "'split' scheme used for solving the incompressible Navier-Stokes "
                             "equations.");
  // Coupled variables
  params.addRequiredCoupledVar("u", "star velocity vector");

  // Optional parameters
  params.addParam<MaterialPropertyName>("rho_name", "rho", "density_name");

  return params;
}

ADNavStokesPressurePoisson::ADNavStokesPressurePoisson(const InputParameters & parameters)
  : ADKernel(parameters),

    // Gradients
    _grad_velocity(adCoupledVectorGradient("u")),
    // Material properties
    _rho(getADMaterialProperty<Real>("rho_name"))
{
}

ADReal
ADNavStokesPressurePoisson::computeQpResidual()
{
  // Return the result
  return _grad_u[_qp] * _grad_test[_i][_qp] +
  (_rho[_qp] / _dt) * _grad_velocity[_qp].tr() * _test[_i][_qp];
}
