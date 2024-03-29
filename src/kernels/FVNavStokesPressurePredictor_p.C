//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVNavStokesPressurePredictor_p.h"
#include "MooseMesh.h"

registerMooseObject("AirfoilAppApp", FVNavStokesPressurePredictor_p);

InputParameters
FVNavStokesPressurePredictor_p::validParams()
{

  auto params = FVFluxKernel::validParams();

  params.addClassDescription("Object pressure predictor equation.");
  // Coupled variables
  // params.addRequiredCoupledVar("pressure", "The pressure variable.");
  // params.addRequiredCoupledVar("pressure_old", "The pressure old variable.");
  // params.addRequiredCoupledVar("u_star", "The velocity in the x direction.");
  // params.addCoupledVar("v_star", "The velocity in the y direction.");
  // params.addCoupledVar("w_star", "The velocity in the z direction.");
  params.addRequiredCoupledVar("Ainv_x", "Ainv_x from momenutm predictor.");
  params.addCoupledVar("Ainv_y", "Ainv_y from momenutm predictor.");
  params.addCoupledVar("Ainv_z", "Ainv_z from momenutm predictor.");
  params.addRequiredParam<MooseFunctorName>("Hu_x", "The x-component functor.");
  params.addParam<MooseFunctorName>("Hu_y", 0, "The y-component functor.");
  params.addParam<MooseFunctorName>("Hu_z", 0, "The z-component functor.");
  // params.addRequiredCoupledVar("rhs_x", "rhs_x from momenutm predictor.");
  // params.addCoupledVar("rhs_y", "rhs_y from momenutm predictor.");
  // params.addCoupledVar("rhs_z", "rhs_z from momenutm predictor.");

  // Set velocity interpolation method for the RSH term
  // MooseEnum velocity_interp_method("average rc", "rc");
  //
  // params.addParam<MooseEnum>(
  //     "velocity_interp_method",
  //     velocity_interp_method,
  //     "The interpolation to use for the velocity. Options are "
  //     "'average' and 'rc' which stands for Rhie-Chow. The default is Rhie-Chow.");

  return params;
}

FVNavStokesPressurePredictor_p::FVNavStokesPressurePredictor_p(const InputParameters & parameters)
  : // Base class
    FVFluxKernel(parameters),

    // Material properties
    // _rho(getFunctor<ADReal>("rho")),
    // _mu(getFunctor<ADReal>("mu")),

    // Get coupled variables
    // _p_var(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure", 0))),
    // _p_old(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure_old", 0))),
    // _u_star(dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("u_star", 0))),
    // _v_star(_mesh.dimension() >= 2 ? dynamic_cast<const INSFVVelocityVariable
    // *>(getFieldVar("v_star", 0)) : nullptr), _w_star(_mesh.dimension() == 3 ? dynamic_cast<const
    // INSFVVelocityVariable *>(getFieldVar("w_star", 0)) : nullptr),

    // Get coupled predictor variables
    _Ainv_x(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Ainv_x", 0))),
    _Ainv_y(_mesh.dimension() >= 2
                ? dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Ainv_y", 0))
                : nullptr),
    _Ainv_z(_mesh.dimension() >= 3
                ? dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Ainv_z", 0))
                : nullptr),

    _Hu_x(getFunctor<ADReal>("Hu_x")),
    _Hu_y(getFunctor<ADReal>("Hu_y")),
    _Hu_z(getFunctor<ADReal>("Hu_z")),

    // _rhs_x(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("rhs_x", 0))),
    // _rhs_y(_mesh.dimension() >= 2 ? dynamic_cast<const MooseVariableFVReal
    // *>(getFieldVar("rhs_y", 0)) : nullptr), _rhs_z(_mesh.dimension() >= 3 ? dynamic_cast<const
    // MooseVariableFVReal *>(getFieldVar("rhs_z", 0)) : nullptr),

    // Variable numberings
    // _u_vel_star_var_number(coupled("u_star")),
    // _v_vel_star_var_number(_mesh.dimension() >= 2 ? coupled("v_star") : libMesh::invalid_uint),
    // _w_vel_star_var_number(_mesh.dimension() == 3 ? coupled("w_star") : libMesh::invalid_uint),

    // Get current element
    _current_elem(_assembly.elem()),
    counter(0)

{
  std::cout << "Constructor OK" << std::endl;
}

ADReal
FVNavStokesPressurePredictor_p::computeQpResidual()
{

  /// Diffusion term
  using namespace Moose::FV;

  const Elem * const elem = &_face_info->elem();
  const Elem * const neighbor = _face_info->neighborPtr();

  ADRealVectorValue elem_Ainv(_Ainv_x->getElemValue(elem) * _face_info->elemVolume());
  if (_Ainv_y)
    elem_Ainv(1) = _Ainv_y->getElemValue(elem) * _face_info->elemVolume();
  if (_Ainv_z)
    elem_Ainv(2) = _Ainv_z->getElemValue(elem) * _face_info->elemVolume();

  mooseAssert((neighbor == &_face_info->elem()) || (neighbor == _face_info->neighborPtr()),
              "Surely the neighbor has to match one of the face information's elements, right?");

  ADRealVectorValue neighbor_Ainv(_Ainv_x->getNeighborValue(neighbor, *_face_info, elem_Ainv(0)) * _face_info->neighborVolume());
  if (_Ainv_y)
    neighbor_Ainv(1) = _Ainv_y->getNeighborValue(neighbor, *_face_info, elem_Ainv(1)) * _face_info->neighborVolume();
  if (_Ainv_z)
    neighbor_Ainv(2) = _Ainv_z->getNeighborValue(neighbor, *_face_info, elem_Ainv(2)) * _face_info->neighborVolume();

  ADRealVectorValue interp_Ainv_face;
  if (onBoundary(*_face_info))
  {
    interp_Ainv_face(0) = _Ainv_x->getBoundaryFaceValue(*_face_info) * _assembly.elemVolume();
    if (_Ainv_y)
      interp_Ainv_face(1) = _Ainv_y->getBoundaryFaceValue(*_face_info) * _assembly.elemVolume();
    if (_Ainv_z)
      interp_Ainv_face(2) = _Ainv_z->getBoundaryFaceValue(*_face_info) * _assembly.elemVolume();
  }
  else
    Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
                           interp_Ainv_face,
                           elem_Ainv,
                           neighbor_Ainv,
                           *_face_info,
                           true);

  ADRealVectorValue Ainv_gradp(interp_Ainv_face(0) * _var.adGradSln(*_face_info)(0));
  if (_Ainv_y)
    Ainv_gradp(1) = interp_Ainv_face(1) * _var.adGradSln(*_face_info)(1);
  if (_Ainv_z)
    Ainv_gradp(2) = interp_Ainv_face(2) * _var.adGradSln(*_face_info)(2);

  ADReal residual = -Ainv_gradp * _face_info->normal(); //* interp_Ainv_face;

  /// Divergence term
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
    ADRealVectorValue face_vec(_Hu_x(cd_face), _Hu_y(cd_face), _Hu_z(cd_face));
    divergence += face_vec * surface_vector;
  };

  Moose::FV::loopOverElemFaceInfo(*_current_elem, _mesh, _subproblem, action_functor);
  residual += divergence; // / _assembly.elemVolume();

  return residual;
}
