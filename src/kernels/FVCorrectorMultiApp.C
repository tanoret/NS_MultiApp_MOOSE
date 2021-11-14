//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVCorrectorMultiApp.h"
#include "MooseMesh.h"

registerMooseObject("AirfoilAppApp", FVCorrectorMultiApp);

InputParameters
FVCorrectorMultiApp::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription(
      "Computes corrector for integration with MultiApp");

  // Coupled variables
  params.addRequiredCoupledVar("pressure", "The pressure variable.");
  params.addRequiredCoupledVar("pressure_old", "The pressure variable.");
  params.addRequiredCoupledVar("Ainv", "Ainv from momenutm predictor.");
  params.addRequiredCoupledVar("Hhat", "Hu from momenutm predictor.");

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

FVCorrectorMultiApp::FVCorrectorMultiApp(const InputParameters & parameters)
  : FVElementalKernel(parameters),

  // Get coupled variables
  _p_var(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure", 0))),
  _p_old(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure_old", 0))),

  // Get coupled predictor variables
  _Ainv(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Ainv", 0))),
  _Hhat(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Hhat", 0))),

  // Get current direction
  _index(getParam<MooseEnum>("momentum_component")),

  // Get pressure relaxation factor
  _advection_relaxation(getParam<Real>("advection_relaxation"))
{
}

ADReal
FVCorrectorMultiApp::computeQpResidual()
{
  /// Diffusion term
  using namespace Moose::FV;

  auto new_pressure_grad = _p_var->adGradSln(_current_elem)(_index);

  std::cout << "P: " << _p_var->getElemValue(_current_elem) << " - Grad: " << new_pressure_grad << std::endl;

  // Computing RHS term
  auto _p_term = _Ainv->getElemValue(_current_elem) * new_pressure_grad * _assembly.elemVolume();

  // Assign new expression because we will be returning its value
  auto _new_vel = _Hhat->getElemValue(_current_elem) - _p_term;

  // std::cout << "Pressure: " << _p_var->getElemValue(_current_elem) << std::endl;
  // std::cout << "Pressure Grad: " << _p_var->adGradSln(_current_elem)(_index) << std::endl;
  // std::cout << "Ainv Pressure Grad: " << _p_term << std::endl;
  // std::cout << "Hhat: " << _Hhat->getElemValue(_current_elem) << std::endl;
  // std::cout << "New Vel: " << _new_vel << std::endl;
  // std::cout << "Relax: " << _advection_relaxation << std::endl;

  // std::cout << "Hhat: " << _Hhat->getElemValue(_current_elem).value()
  //           << "Pressure: " << _p_var->getElemValue(_current_elem).value()
  //           << "Pressure grad: " << new_pressure_grad.value()
  //           << " Ainv Press grad: " << _p_term.value()
  //           << " New vel: " << _new_vel.value() << "\n" << std::endl;

  // Current solution = _new_vel

  auto residual = _advection_relaxation * _new_vel
                 + (1 - _advection_relaxation) * _var.getElementalValueOld(_current_elem)
                 - _u[_qp];

  //std::cout << "Res. der.: " << residual.derivatives() << std::endl;
  // std::cout << "Velocity adv: "
  // << _advection_relaxation * _new_vel + (1 - _advection_relaxation) * _var.getElementalValueOld(_current_elem)
  // << std::endl;

  return residual;
}

void
FVCorrectorMultiApp::computeResidual()
{
  prepareVectorTag(_assembly, _var.number());
  _local_re(0) += MetaPhysicL::raw_value(computeQpResidual());
  //std::cout << "local re: " << _local_re(0) << std::endl;
  accumulateTaggedLocalResidual();
}
