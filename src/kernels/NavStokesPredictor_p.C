//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NavStokesPredictor_p.h"
#include "MooseMesh.h"

registerMooseObject("AirfoilAppApp", NavStokesPredictor_p);

InputParameters
NavStokesPredictor_p::validParams()
{
  InputParameters params = Kernel::validParams();

  params.addClassDescription("This class computes the Predictor equation in "
                             "fully-discrete (both time and space) form.");
  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", "z-velocity"); // only required in 3D

  // Make star also be required, even though we might not use it?
  params.addRequiredCoupledVar("u_star", "star x-velocity");
  params.addCoupledVar("v_star", "star y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w_star", "star z-velocity"); // only required in 3D

  params.addRequiredCoupledVar("p", "pressure");

  // Required parameters
  params.addRequiredRangeCheckedParam<unsigned>(
      "component",
      "component>=0 & component<=2",
      "0,1,2 depending on if we are solving the x,y,z component of the Predictor equation");

  // Optional parameters
  params.addParam<MaterialPropertyName>("mu_name", "mu", "The name of the dynamic viscosity");
  params.addParam<MaterialPropertyName>("rho_name", "rho", "The name of the density");

  return params;
}

NavStokesPredictor_p::NavStokesPredictor_p(const InputParameters & parameters)
  : Kernel(parameters),

      // Old velocities
      _u_vel_old(coupledValue("u")),
      _v_vel_old(_mesh.dimension() >= 2 ? coupledValue("v") : _zero),
      _w_vel_old(_mesh.dimension() == 3 ? coupledValue("w") : _zero),

      // Star velocities
      _u_vel_star(coupledValue("u_star")),
      _v_vel_star(_mesh.dimension() >= 2 ? coupledValue("v_star") : _zero),
      _w_vel_star(_mesh.dimension() == 3 ? coupledValue("w_star") : _zero),

      // Star Velocity Gradients
      _grad_u_vel_star(coupledGradient("u_star")),
      _grad_v_vel_star(_mesh.dimension() >= 2 ? coupledGradient("v_star") : _grad_zero),
      _grad_w_vel_star(_mesh.dimension() == 3 ? coupledGradient("w_star") : _grad_zero),

      // Pressure Gradient previous timestep
      _grad_p_old(coupledGradient("p")),

      // Star velocity numberings
      _u_vel_star_var_number(coupled("u_star")),
      _v_vel_star_var_number(_mesh.dimension() >= 2 ? coupled("v_star") : libMesh::invalid_uint),
      _w_vel_star_var_number(_mesh.dimension() == 3 ? coupled("w_star") : libMesh::invalid_uint),

      // Required parameters
      _component(getParam<unsigned>("component")),

      // Material properties
      _mu(getMaterialProperty<Real>("mu_name")),
      _rho(getMaterialProperty<Real>("rho_name"))
{
}

Real
NavStokesPredictor_p::computeQpResidual()
{
  // Vector object for test function
  RealVectorValue test;
  test(_component) = _test[_i][_qp];

  // Tensor object for test function gradient
  RealTensorValue grad_test;
  for (unsigned k = 0; k < 3; ++k)
    grad_test(_component, k) = _grad_test[_i][_qp](k);

  // Note: _u is the component'th entry of "u_star" in Chorin's method.
  RealVectorValue U(_u_vel_star[_qp], _v_vel_star[_qp], _w_vel_star[_qp]);
  RealTensorValue grad_U(_grad_u_vel_star[_qp], _grad_v_vel_star[_qp], _grad_w_vel_star[_qp]);
  RealVectorValue U_old(_u_vel_old[_qp], _v_vel_old[_qp], _w_vel_old[_qp]);

  Real time_derivative = (U(_component) - U_old(_component)) * _test[_i][_qp];

  // Convective part.  Remember to multiply by _dt!
  Real convective_part =  (grad_U * U_old) * test;

  // Pressure Gradients
  Real pressure_gradient = (1 / _rho[_qp]) * _grad_p_old[_qp](_component) * _test[_i][_qp];

  // Viscous part - we are using the Laplacian form here for simplicity.
  // Remember to multiply by _dt!
  Real viscous_part =  (_mu[_qp] / _rho[_qp]) * grad_U.contract(grad_test);

  return time_derivative + (convective_part + viscous_part + pressure_gradient) * _dt;
}

Real
NavStokesPredictor_p::computeQpJacobian()
{
  // The mass matrix part is always there.
  Real mass_part = _phi[_j][_qp] * _test[_i][_qp];

  // The on-diagonal Jacobian contribution depends on whether the predictor uses the
  // 'new' or 'star' velocity.
  Real other_part = 0.;

  RealVectorValue U_old(_u_vel_old[_qp], _v_vel_old[_qp], _w_vel_old[_qp]);

  Real convective_part =
      _dt * ((U_old * _grad_phi[_j][_qp])) * _test[_i][_qp];

  Real viscous_part =
      _dt * ((_mu[_qp] / _rho[_qp]) * (_grad_phi[_j][_qp] * _grad_test[_i][_qp]));

  other_part = convective_part + viscous_part;

  return mass_part + other_part;
}

Real
NavStokesPredictor_p::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_star_var_number)
  {
    return _dt * _phi[_j][_qp] * _grad_u[_qp](0) * _test[_i][_qp];
  }

  else if (jvar == _v_vel_star_var_number)
  {
    return _dt * _phi[_j][_qp] * _grad_u[_qp](1) * _test[_i][_qp];
  }

  else if (jvar == _w_vel_star_var_number)
  {
    return _dt * _phi[_j][_qp] * _grad_u[_qp](2) * _test[_i][_qp];
  }

  else
    return 0;

}
