//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVNavierStokesTimeRelax.h"

#include "SystemBase.h"

registerADMooseObject("MooseApp", FVNavierStokesTimeRelax);

InputParameters
FVNavierStokesTimeRelax::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription(
      "Residual contribution from time derivative of a variable for the finite volume method.");
  params.set<MultiMooseEnum>("vector_tags") = "time";
  params.set<MultiMooseEnum>("matrix_tags") = "system time";

  params.addRangeCheckedParam<Real>("velocity_relaxation", 1.0,
                                    "(0 <= velocity_relaxation) & (velocity_relaxation <= 1)",
                                    "Pressure field relaxation must be between 0 and 1 both included.");

  params.addDeprecatedParam<bool>("add_time_derivative",
                                  false,
                                  "Whether or not to check for steady state conditions",
                                  "Use steady_state_detection instead");

  params.addRequiredCoupledVar("Ainv", "Ainv from previous momenutm predictor for solver relaxation.");
  params.addRequiredCoupledVar("Ainv_old", "Ainv_old from previous momenutm predictor for solver relaxation.");

  return params;
}

FVNavierStokesTimeRelax::FVNavierStokesTimeRelax(const InputParameters & parameters)
  : FVElementalKernel(parameters),
    _add_time_derivative(getParam<bool>("add_time_derivative")),
    _velocity_relaxation(getParam<Real>("velocity_relaxation")),
    _Ainv(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Ainv", 0))),
    _Ainv_old(dynamic_cast<const MooseVariableFVReal *>(getFieldVar("Ainv_old", 0))),
    _u_dot(_var.adUDot())
{
  std::cout << "Constructor OK" << std::endl;
}

ADReal
FVNavierStokesTimeRelax::computeQpResidual()
{
  ADReal residual(0.0);
  Real _loc_vel_relaxation(0.0);
  Real _diagonal_entry(0.0);
  Real _diagonal_entry_old(0.0);

  if (_add_time_derivative)
    residual += _u_dot[_qp];
  //
  // Implementing lagged relaxation
  // Should be fine for steady state and we don't really care for transient
  if (_velocity_relaxation < 1.0)
  {
    _loc_vel_relaxation = (1.0 - _velocity_relaxation) / _velocity_relaxation;
    auto _inverse_diagonal_entry = _Ainv->getElemValue(_current_elem); //.value();
    auto _inverse_diagonal_entry_old = _Ainv_old->getElementalValueOld(_current_elem);
    auto _elem_trace_min = std::pow(3.0, 1); // /std::sqrt(_assembly.elemVolume());
    auto _elem_trace_max = std::pow(3.0, 3); // /std::sqrt(_assembly.elemVolume());

    if (_inverse_diagonal_entry == 0)
    {
      _diagonal_entry = _elem_trace_min;
    }
    else
    {
      _diagonal_entry = std::max(std::min(1.0 / std::abs(_inverse_diagonal_entry.value()),
                                          _elem_trace_max),
                                          _elem_trace_min);
    }
    if (_inverse_diagonal_entry_old == 0)
    {
      _diagonal_entry_old = _elem_trace_min;
    }
    else
    {
      _diagonal_entry = std::max(std::min(1.0 / std::abs(_inverse_diagonal_entry_old),
                                          _elem_trace_max),
                                          _elem_trace_min);
    }

    // std::cout << "Diag: " << _loc_vel_relaxation * _diagonal_entry << std::endl;
    // std::cout << "Diag Old: " << _loc_vel_relaxation * _diagonal_entry_old << std::endl;
    // std::cout << "Diff: " << _loc_vel_relaxation * (_diagonal_entry - _diagonal_entry_old) << std::endl;

    //std::cout << _diagonal_entry << std::endl;
    residual += _loc_vel_relaxation * _diagonal_entry
                * (_u[_qp] - _var.getElementalValueOld(_current_elem));
  }

  return residual;

}
