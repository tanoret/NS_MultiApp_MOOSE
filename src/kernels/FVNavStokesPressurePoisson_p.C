//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVNavStokesPressurePoisson_p.h"
#include "MooseMesh.h"

registerMooseObject("AirfoilAppApp", FVNavStokesPressurePoisson_p);

InputParameters
FVNavStokesPressurePoisson_p::validParams()
{

  auto params = FVFluxKernel::validParams();

  params.addRequiredParam<MooseFunctorName>("mu", "The viscosity functor material property");
  params.addRequiredParam<MaterialPropertyName>("rho", "Density functor material property");

  params.addClassDescription("Object pressure Poisson eqaution");
  // Coupled variables
  params.addRequiredCoupledVar("u_star", "star x-velocity");
  params.addCoupledVar("v_star", "star y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w_star", "star z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("pressure_old", "pressure old");

  // Set velocity interpolation method for the RSH term
  MooseEnum velocity_interp_method("average rc", "rc");

  params.addParam<MooseEnum>(
      "velocity_interp_method",
      velocity_interp_method,
      "The interpolation to use for the velocity. Options are "
      "'average' and 'rc' which stands for Rhie-Chow. The default is Rhie-Chow.");

  return params;
}

FVNavStokesPressurePoisson_p::FVNavStokesPressurePoisson_p(const InputParameters & parameters)
  :
    // Base class
    FVFluxKernel(parameters),

    // Material properties
    _rho(getFunctor<ADReal>("rho")),
    _mu(getFunctor<ADReal>("mu")),

    // Get coupled variables
    _p_old(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure_old", 0))),
    _u_star(dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("u_star", 0))),
    _v_star(_mesh.dimension() >= 2 ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("v_star", 0)) : nullptr),
    _w_star(_mesh.dimension() == 3 ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("w_star", 0)) : nullptr),

    // Variable numberings
    _u_vel_star_var_number(coupled("u_star")),
    _v_vel_star_var_number(_mesh.dimension() >= 2 ? coupled("v_star") : libMesh::invalid_uint),
    _w_vel_star_var_number(_mesh.dimension() == 3 ? coupled("w_star") : libMesh::invalid_uint),

    // Get current element
    _current_elem(_assembly.elem())

{
}


ADReal
FVNavStokesPressurePoisson_p::computeQpResidual()
{

  /// Diffusion term
  using namespace Moose::FV;

  const auto elem_face = elemFromFace();
  const auto neighbor_face = neighborFromFace();
  const auto mu_elem = _mu(elem_face);
  const auto mu_neighbor = _mu(neighbor_face);

  // Compute the diffusion driven by the velocity gradient
  // Interpolate viscosity divided by porosity on the face
  ADReal mu_face;
  interpolate(Moose::FV::InterpMethod::Average,
              mu_face,
              mu_elem,
              mu_neighbor,
              *_face_info,
              true);

  // Compute face superficial velocity gradient
  auto dudn = gradUDotNormal();

  // First term of residual
  ADReal residual = mu_face * dudn;

  /// Divergence term
  const auto _grad_u_x = _u_star->adGradSln(_current_elem)(0);
  auto u_div = _grad_u_x;

  if (_mesh.dimension() >= 2)
  {
    const auto _grad_v_x = _v_star->adGradSln(_current_elem)(1);
    u_div = u_div + _grad_v_x;
  }
  if (_mesh.dimension() >= 2)
  {
    const auto _grad_w_z = _w_star->adGradSln(_current_elem)(2);
    u_div = u_div + _grad_w_z;
  }

  if(_is_transient)
    residual -= (_rho(_current_elem)/_dt) * u_div;

  return residual;
}
