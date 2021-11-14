//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVCorrector.h"
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

registerMooseObject("AirfoilAppApp", FVCorrector);

InputParameters
FVCorrector::validParams()
{
  auto params = AuxKernel::validParams();

  params.addClassDescription("Object pressure predictor equation.");
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

  params.addRangeCheckedParam<Real>("pressure_relaxation", 1.0,
                                    "(0 <= pressure_relaxation) & (pressure_relaxation <= 1)",
                                    "Pressure field relaxation must be between 0 and 1 both included.");

  return params;
}

FVCorrector::FVCorrector(const InputParameters & parameters)
  :
    // Base class
    AuxKernel(parameters),

    // Get coupled variables
    _p_var(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure", 0))),
    _p_old(dynamic_cast<const INSFVPressureVariable *>(getFieldVar("pressure_old", 0))),
    //_vel_star(dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("vel_star", 0))),

    // Get coupled predictor variables
    _Ainv(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Ainv", 0))),
    _Hhat(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Hhat", 0))),

    // Get current element
    //_current_elem(_assembly.elem()),

    // Get current direction
    _index(getParam<MooseEnum>("momentum_component")),

    // Get pressure relaxation factor
    _pressure_relaxation(getParam<Real>("pressure_relaxation")),
    _u_old(uOld())
{
}

Real
FVCorrector::computeValue()
{
  /// Diffusion term
  using namespace Moose::FV;

  auto new_pressure_grad = _p_var->adGradSln(_current_elem)(_index);

  std::cout << "Grad: " << new_pressure_grad << std::endl;

  // Computing RHS term
  auto _p_term = _Ainv->getElemValue(_current_elem) * new_pressure_grad * _assembly.elemVolume();

  // Assign new expression because we will be returning its value
  auto _new_vel = _Hhat->getElemValue(_current_elem) - _p_term;

  auto _old_vel = _u_old[_qp];

  std::cout << "Hhat: " << _Hhat->getElemValue(_current_elem).value()
            << " Press grad: " << _p_term.value()
            << " New vel: " << _new_vel.value() << std::endl;

  return (_pressure_relaxation) * _new_vel.value() +
         (1. - _pressure_relaxation) * _old_vel;
}
