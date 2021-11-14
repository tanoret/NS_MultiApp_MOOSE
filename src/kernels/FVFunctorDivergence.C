//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVFunctorDivergence.h"
#include "Function.h"
#include "FVUtils.h"

registerMooseObject("AirfoilAppApp", FVFunctorDivergence);

InputParameters
FVFunctorDivergence::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription(
      "Computes the divergence of a vector field composed of component functors");
  params.addRequiredParam<Real>("sign", "The sign of the divergence term");
  params.addRequiredParam<MooseFunctorName>("x_functor", "The x-component functor.");
  params.addParam<MooseFunctorName>("y_functor", 0, "The y-component functor.");
  params.addParam<MooseFunctorName>("z_functor", 0, "The z-component functor.");
  return params;
}

FVFunctorDivergence::FVFunctorDivergence(const InputParameters & parameters)
  : FVElementalKernel(parameters),
    _sign(getParam<Real>("sign")),
    _x(getFunctor<ADReal>("x_functor")),
    _y(getFunctor<ADReal>("y_functor")),
    _z(getFunctor<ADReal>("z_functor"))
    //_current_elem(_assembly.elem())
{
}

ADReal
FVFunctorDivergence::computeQpResidual()
{
  ADReal divergence = 0;

  auto action_functor = [this, &divergence](const Elem & libmesh_dbg_var(functor_elem),
                                            const Elem * neighbor,
                                            const FaceInfo * const fi,
                                            const Point & surface_vector,
                                            Real,
                                            bool) {
    mooseAssert(&functor_elem == _current_elem, "these should be the same");
    const auto elem_sub_id = _current_elem->subdomain_id();
    const auto neighbor_sub_id = neighbor ? neighbor->subdomain_id() : Moose::INVALID_BLOCK_ID;
    const auto cd_face = std::make_tuple(fi,
                                         Moose::FV::LimiterType::CentralDifference,
                                         true,
                                         std::make_pair(elem_sub_id, neighbor_sub_id));
    ADRealVectorValue face_vec(_x(cd_face), _y(cd_face), _z(cd_face));
    std::cout << _current_elem->vertex_average() << fi->faceCentroid() << face_vec * surface_vector << std::endl;
    divergence += face_vec * surface_vector;
  };

  Moose::FV::loopOverElemFaceInfo(*_current_elem, _mesh, _subproblem, action_functor);
  return _sign * divergence; /// _assembly.elemVolume();
}
