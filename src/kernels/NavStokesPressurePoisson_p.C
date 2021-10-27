//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NavStokesPressurePoisson_p.h"
#include "MooseMesh.h"

registerMooseObject("AirfoilAppApp", NavStokesPressurePoisson_p);

InputParameters
NavStokesPressurePoisson_p::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addClassDescription("This class computes the pressure Poisson solve which is part of the "
                             "'split' scheme used for solving the incompressible Navier-Stokes "
                             "equations.");
  // Coupled variables
  params.addRequiredCoupledVar("u_star", "star x-velocity");
  params.addCoupledVar("v_star", "star y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w_star", "star z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("p", "pressure");

  // Optional parameters
  params.addParam<MaterialPropertyName>("rho_name", "rho", "density_name");

  return params;
}

NavStokesPressurePoisson_p::NavStokesPressurePoisson_p(const InputParameters & parameters)
  : Kernel(parameters),

    // Gradients
    _grad_u_star(coupledGradient("u_star")),
    _grad_v_star(_mesh.dimension() >= 2 ? coupledGradient("v_star") : _grad_zero),
    _grad_w_star(_mesh.dimension() == 3 ? coupledGradient("w_star") : _grad_zero),

    _grad_p_old(coupledGradientOld("p")),

    // Variable numberings
    _u_vel_star_var_number(coupled("u_star")),
    _v_vel_star_var_number(_mesh.dimension() >= 2 ? coupled("v_star") : libMesh::invalid_uint),
    _w_vel_star_var_number(_mesh.dimension() == 3 ? coupled("w_star") : libMesh::invalid_uint),

    // Material properties
    _rho(getMaterialProperty<Real>("rho_name"))
{
}

Real
NavStokesPressurePoisson_p::computeQpResidual()
{
  // Laplacian part
  Real laplacian_part = _grad_u[_qp] * _grad_test[_i][_qp];

  Real laplacian_part_old =  - _grad_p_old[_qp] * _grad_test[_i][_qp];

  // Divergence part, don't forget to *divide* by _dt
  Real div_part = (_rho[_qp] / _dt) *
                  (_grad_u_star[_qp](0) + _grad_v_star[_qp](1) + _grad_w_star[_qp](2)) *
                  _test[_i][_qp];

  // Return the result
  return laplacian_part + laplacian_part_old + div_part;
}

Real
NavStokesPressurePoisson_p::computeQpJacobian()
{
  return _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

Real
NavStokesPressurePoisson_p::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_star_var_number)
    return (_rho[_qp] / _dt) * _grad_phi[_j][_qp](0) * _test[_i][_qp];

  else if (jvar == _v_vel_star_var_number)
    return (_rho[_qp] / _dt) * _grad_phi[_j][_qp](1) * _test[_i][_qp];

  else if (jvar == _w_vel_star_var_number)
    return (_rho[_qp] / _dt) * _grad_phi[_j][_qp](2) * _test[_i][_qp];

  else
    return 0;
}
