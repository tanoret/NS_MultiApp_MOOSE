//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVHhat.h"
#include "MooseMesh.h"

#include "FEProblem.h"
#include "SubProblem.h"
#include "AuxiliarySystem.h"
#include "MooseTypes.h"
#include "Assembly.h"

#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/quadrature.h"
#include "libmesh/boundary_info.h"

registerMooseObject("AirfoilAppApp", FVHhat);

InputParameters
FVHhat::validParams()
{
  auto params = AuxKernel::validParams();

  params.addClassDescription("Object pressure predictor equation.");
  // Coupled variables
  params.addRequiredCoupledVar("pressure", "The pressure old variable.");
  params.addRequiredCoupledVar("Ainv", "Ainv from momenutm predictor.");
  params.addRequiredCoupledVar("Hu", "Hu from momenutm predictor.");
  params.addRequiredCoupledVar("rhs", "rhs from momenutm predictor.");

  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>(
      "momentum_component",
      momentum_component,
      "The component of the momentum equation that this kernel applies to.");

  return params;
}

FVHhat::FVHhat(const InputParameters & parameters)
  :
    // Base class
    AuxKernel(parameters),

    // Get coupled variables
    _p_mom_predictor(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure", 0))),

    // Get coupled predictor variables
    _Ainv(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Ainv", 0))),
    _Hu(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Hu", 0))),
    _rhs(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("rhs", 0))),

    // Get current element
    //_current_elem(_assembly.elem()),

    // Get current direction
    _index(getParam<MooseEnum>("momentum_component"))
{
}

Real
FVHhat::computeValue()
{
  /// Diffusion term
  using namespace Moose::FV;

  // auto _rhs_corr = _rhs->getElemValue(_current_elem) +
  //                  _Ainv->getElemValue(_current_elem) * _p_mom_predictor->adGradSln(_current_elem)(_index) * _assembly.elemVolume();

  // Constructing Hhat
  // auto _Hhat = - 1.0 * _Ainv->getElemValue(_current_elem) * _Hu->getElemValue(_current_elem)
  //              + _Ainv->getElemValue(_current_elem) * _rhs_corr;

 auto _Hhat = _rhs->getElemValue(_current_elem) - _Hu->getElemValue(_current_elem);

  // std::cout << "RHS: " << _rhs->getElemValue(_current_elem).value() << std::endl;
  // std::cout << "pgrad: " << (_p_mom_predictor->adGradSln(_current_elem)(_index) * _assembly.elemVolume()).value() << std::endl;
  // std::cout << "AinHu: " << (- 1.0 * _Ainv->getElemValue(_current_elem) * _Hu->getElemValue(_current_elem)).value() << std::endl;
  // std::cout << "rhs corr: " << (_Ainv->getElemValue(_current_elem) * _rhs_corr).value() << std::endl;
  // std::cout << "H hat: " << _Hhat.value() << std::endl;

  return _Hhat.value();
}
