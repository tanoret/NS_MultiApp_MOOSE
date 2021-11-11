//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVNavStokesHhatSource.h"
#include "MooseMesh.h"

registerMooseObject("AirfoilAppApp", FVNavStokesHhatSource);

InputParameters
FVNavStokesHhatSource::validParams()
{

  auto params = FVFluxKernel::validParams();

  params.addClassDescription("Object pressure predictor equation.");
  // Coupled variables
  params.addRequiredCoupledVar("Hu_x", "Hu_x from momenutm predictor.");
  params.addCoupledVar("Hu_y", "Hu_y from momenutm predictor.");
  params.addCoupledVar("Hu_z", "Hu_z from momenutm predictor.");

  return params;
}

FVNavStokesHhatSource::FVNavStokesHhatSource(const InputParameters & parameters)
  :
    // Base class
    FVFluxKernel(parameters),

    _Hu_x(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Hu_x", 0))),
    _Hu_y(_mesh.dimension() >= 2 ? dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Hu_y", 0)) : nullptr),
    _Hu_z(_mesh.dimension() >= 3 ? dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Hu_z", 0)) : nullptr),

    // Get current element
    _current_elem(_assembly.elem())

{
  std::cout << "Constructor OK" << std::endl;

}


ADReal
FVNavStokesHhatSource::computeQpResidual()
{

  /// Diffusion term
  using namespace Moose::FV;

  /// Divergence term
  const Elem * const elem = &_face_info->elem();
  const Elem * const neighbor = _face_info->neighborPtr();

  // Interpolating H to the faces
  ADRealVectorValue elem_Hu(_Hu_x->getElemValue(elem));
  if (_Hu_y)
    elem_Hu(1) = _Hu_y->getElemValue(elem);
  if (_Hu_z)
    elem_Hu(2) = _Hu_z->getElemValue(elem);

  ADRealVectorValue neighbor_Hu(_Hu_x->getNeighborValue(neighbor, *_face_info, elem_Hu(0)));
  if (_Hu_y)
    neighbor_Hu(1) = _Hu_y->getNeighborValue(neighbor, *_face_info, elem_Hu(1));
  if (_Hu_z)
    neighbor_Hu(2) = _Hu_z->getNeighborValue(neighbor, *_face_info, elem_Hu(2));

  ADRealVectorValue interp_Hu_face;
  Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
                         interp_Hu_face,
                         elem_Hu,
                         neighbor_Hu,
                         *_face_info,
                         true);

  auto residual = interp_Hu_face * _face_info->normal();

  return residual;
}
