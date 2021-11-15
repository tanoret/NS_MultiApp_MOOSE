//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVNavierStokesCorrector.h"
#include "MooseMesh.h"

registerMooseObject("AirfoilAppApp", FVNavierStokesCorrector);

InputParameters
FVNavierStokesCorrector::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params.addClassDescription(
      "Computes corrector for integration with MultiApp");

  // Coupled variables
  // params.addRequiredCoupledVar("pressure", "The pressure variable.");
  // params.addRequiredCoupledVar("pressure_old", "The pressure variable.");
  // params.addRequiredCoupledVar("Ainv", "Ainv from momenutm predictor.");
  // params.addRequiredCoupledVar("Hhat", "Hu from momenutm predictor.");
  params.addRequiredParam<MooseFunctorName>("pressure", "TThe pressure variable functor.");
  params.addRequiredParam<MooseFunctorName>("pressure_old", "The pressure old variable functor.");
  params.addRequiredParam<MooseFunctorName>("Ainv", "Ainv from momenutm predictor functor.");
  params.addRequiredParam<MooseFunctorName>("Hhat", "Hhat from momenutm predictor functor.");

  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>(
      "momentum_component",
      momentum_component,
      "The component of the momentum equation that this kernel applies to.");

  params.addRangeCheckedParam<Real>("advection_relaxation", 1.0,
                                    "(0 <= advection_relaxation) & (advection_relaxation <= 1)",
                                    "Pressure field relaxation must be between 0 and 1 both included.");

  return params;
}

FVNavierStokesCorrector::FVNavierStokesCorrector(const InputParameters & parameters)
  : FVFluxKernel(parameters),

  // Get coupled variables
  // _p_var(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure", 0))),
  // _p_old(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure_old", 0))),
  _p_var(getFunctor<ADReal>("pressure")),
  _p_old(getFunctor<ADReal>("pressure_old")),

  // Get coupled predictor variables
  // _Ainv(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Ainv", 0))),
  // _Hhat(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Hhat", 0))),
  _Ainv(getFunctor<ADReal>("Ainv")),
  _Hhat(getFunctor<ADReal>("Hhat")),

  // Get current direction
  _index(getParam<MooseEnum>("momentum_component")),

  // Get pressure relaxation factor
  _advection_relaxation(getParam<Real>("advection_relaxation"))
{
}

ADReal
FVNavierStokesCorrector::computeQpResidual()
{
  /// Diffusion term
  using namespace Moose::FV;

  const Elem * const elem = &_face_info->elem();
  const Elem * const neighbor = _face_info->neighborPtr();

  const auto elem_sub_id = elem->subdomain_id();
  const auto neighbor_sub_id = neighbor ? neighbor->subdomain_id() : Moose::INVALID_BLOCK_ID;
  const auto cd_face = std::make_tuple(_face_info,
                                       Moose::FV::LimiterType::CentralDifference,
                                       true,
                                       std::make_pair(elem_sub_id, neighbor_sub_id));

  ADReal _Ainv_face(_Ainv(cd_face));
  ADReal _p_grad_face(_p_var.gradient(cd_face)(_index));
  ADReal _H_hat_face(_Hhat(cd_face));

  ADReal _new_vel(_H_hat_face - _Ainv_face * _p_grad_face);

  const auto var_interface = _var(cd_face);

  const auto var_old_elem = _var.getElementalValueOld(&_face_info->elem(), 0);
  const auto var_old_neighbor = onBoundary(*_face_info) ? var_old_elem : _var.getElementalValueOld(_face_info->neighborPtr(), 0);
  ADReal old_var_interface;
  Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
                         old_var_interface,
                         var_old_elem,
                         var_old_neighbor,
                         *_face_info,
                         true);

  auto residual = _advection_relaxation * _new_vel
                  + (1 - _advection_relaxation) * old_var_interface
                  - var_interface;

  return residual;


  // ADReal _p_grad_face();
  //
  //
  // ADReal elem_Ainv(_Ainv->getElemValue(elem));
  // ADReal neighbor_Ainv(_Ainv->getNeighborValue(neighbor, *_face_info, elem_Ainv));
  //
  // ADReal interp_Ainv_face;
  // Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
  //                        interp_Ainv_face,
  //                        elem_Ainv,
  //                        neighbor_Ainv,
  //                        *_face_info,
  //                        true);
  //
  // ADReal pressure_contribution = interp_Ainv_face * _p_var->adGradSln(*_face_info)(_index);
  //
  // ADReal elem_Hhat(_Hhat->getElemValue(elem));
  // ADReal neighbor_Hhat(_Hhat->getNeighborValue(neighbor, *_face_info, elem_Hhat));
  //
  // ADReal interp_Hhat_face;
  // Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
  //                        interp_Hhat_face,
  //                        elem_Hhat,
  //                        neighbor_Hhat,
  //                        *_face_info,
  //                        true);
  //
  // auto Hhat_contribution = interp_Hhat_face;
  //
  // // auto new_vel_face = Hhat_contribution - pressure_contribution;
  //
  // const auto face = std::make_tuple(_face_info, Moose::FV::LimiterType::CentralDifference, true, faceArgSubdomains());
  // const auto var_interface = _var(face);
  //
  // const auto var_old_elem = _var.getElementalValueOld(&_face_info->elem(), 0);
  // const auto var_old_neighbor = onBoundary(*_face_info) ? var_old_elem : _var.getElementalValueOld(_face_info->neighborPtr(), 0);
  // ADReal old_var_interface;
  // Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
  //                        old_var_interface,
  //                        var_old_elem,
  //                        var_old_neighbor,
  //                        *_face_info,
  //                        true);
  //
  // auto residual = _advection_relaxation * (Hhat_contribution - pressure_contribution)
  //                 + (1 - _advection_relaxation) * old_var_interface
  //                 - var_interface;
  //
  // return residual;
}
