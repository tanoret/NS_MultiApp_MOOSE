//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVFluxKernel.h"

// Forward Declarations

/**
 * This class computes the pressure Poisson solve which is part of
 * the "split" scheme used for solving the incompressible Navier-Stokes
 * equations.
 */
class FVNavStokesHhatSource: public FVFluxKernel
{
public:
  static InputParameters validParams();

  FVNavStokesHhatSource(const InputParameters & parameters);

  virtual ~FVNavStokesHhatSource() {}

protected:

  // Override QP residual
  virtual ADReal computeQpResidual() override;

  /// Transfer Variables
  const MooseVariableFVReal * const _Hu_x;
  const MooseVariableFVReal * const _Hu_y;
  const MooseVariableFVReal * const _Hu_z;

  //// Access to the current element
  const Elem * const & _current_elem;

};
