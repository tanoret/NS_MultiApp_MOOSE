//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Corrector.h"

registerMooseObject("AirfoilAppApp", Corrector);

InputParameters
Corrector::validParams()
{
  InputParameters params = AuxKernel::validParams();

  // Add a "coupling paramater" to get a variable from the input file.
  params.addRequiredCoupledVar("u_star", "star x-velocity");
  params.addCoupledVar("v_star", "star y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w_star", "star z-velocity"); // only required in 3D

  params.addRequiredCoupledVar("p", "pressure");

  // Required parameters
  params.addRequiredParam<unsigned>(
      "component",
      "0,1,2 depending on if we are solving the x,y,z component of the Corrector equation");

  // Optional parameters
  params.addParam<MaterialPropertyName>("rho_name", "rho", "density name");

  return params;
}

Corrector::Corrector(const InputParameters & parameters)
  : AuxKernel(parameters),

  _u_vel_star(coupledValue("u_star")),
  _v_vel_star(_mesh.dimension() >= 2 ? coupledValue("v_star") : _zero),
  _w_vel_star(_mesh.dimension() == 3 ? coupledValue("w_star") : _zero),

  // Pressure at previous time num_steps
  _grad_p(coupledGradient("p")),
  _grad_p_old(coupledGradientOld("p")),

  _component(getParam<unsigned>("component")),

  // Material properties
  _rho(getMaterialProperty<Real>("rho_name"))

{
}

Real
Corrector::computeValue()
{
  RealVectorValue U_star(_u_vel_star[_qp], _v_vel_star[_qp], _w_vel_star[_qp]);

  return U_star(_component) + (_dt / _rho[_qp]) * (_grad_p_old[_qp](_component) -  _grad_p[_qp](_component));
}
