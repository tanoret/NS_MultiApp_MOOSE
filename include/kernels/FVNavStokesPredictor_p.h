//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#pragma once

#include "FVMatAdvection.h"
#include "FVElementalKernel.h"
#include "FVTimeKernel.h"
#include "FVFluxKernel.h"

#include "TheWarehouse.h"
#include "SubProblem.h"
#include "MooseApp.h"
#include "INSFVAttributes.h"

#include "FVElementalKernel.h"

#include <vector>
#include <set>

class INSFVVelocityVariable;
class INSFVPressureVariable;


class FVNavierStokesPredictor_p : public FVMatAdvection
                                  //public FVTimeKernel
                                  //public FVFluxKernel,
                                  //public FVElementalKernel
{
public:
  static InputParameters validParams();
  FVNavierStokesPredictor_p(const InputParameters & params);
  void initialSetup() override;

protected:

  // Problem dimension
  const unsigned int _dim;

  // Thermophysical properties
  const Moose::Functor<ADReal> & _rho; //Density
  const Moose::Functor<ADReal> & _mu; //Dynamic viscosity

  // Turbulent coefficients
  // const Moose::Functor<ADReal> & _mixing_len;

  /// index x|y|z
  const unsigned int _index;

  // Over-write virtual Qp residual mehtod
  virtual ADReal computeQpResidual() override;
  void residualSetup() override final { clearRCCoeffs(); }
  void jacobianSetup() override final { clearRCCoeffs(); }

  // Interpolation overload for the velocity in the advection kernel
  virtual void interpolate(Moose::FV::InterpMethod m, ADRealVectorValue & interp_v);

  //Rhie-Chow related
  const VectorValue<ADReal> & rcCoeff(const Elem & elem) const;
  virtual VectorValue<ADReal> coeffCalculator(const Elem & elem) const;
  void clearRCCoeffs();
  bool skipForBoundary(const FaceInfo & fi) const override;

  /// pressure variable
  const INSFVPressureVariable * const _p_var;
  /// x-velocity
  const INSFVVelocityVariable * const _u_var;
  /// y-velocity
  const INSFVVelocityVariable * const _v_var;
  /// z-velocity
  const INSFVVelocityVariable * const _w_var;
  // pressure kernel
  // const MooseVariableFVReal * const _p_mom_var;
  /// variable for time derivative
  // const ADVariableValue & _var_dot;

  /// The interpolation method to use for the velocity
  Moose::FV::InterpMethod _velocity_interp_method;

  /// Boundary IDs
  std::set<BoundaryID> _no_slip_wall_boundaries;
  std::set<BoundaryID> _slip_wall_boundaries;
  std::set<BoundaryID> _flow_boundaries;
  std::set<BoundaryID> _fully_developed_flow_boundaries; //sub-set of flow boundary IDs
  std::set<BoundaryID> _symmetry_boundaries;
  std::set<BoundaryID> _all_boundaries; // All the BoundaryIDs covered by our different types of INSFVBCs

  /// A map from elements to the 'a' coefficients used in the Rhie-Chow interpolation
  static std::unordered_map<const MooseApp *, std::vector<std::unordered_map<const Elem *, VectorValue<ADReal>>>> _rc_a_coeffs;

  //// Access to the current element
  const Elem * const & _current_elem;

private:
  void setupFlowBoundaries(BoundaryID bnd_id); //Query for \p INSFVBCs::INSFVFlowBC on \p bc_id and add if query successful

  template <typename T>
  void setupBoundaries(const BoundaryID bnd_id, INSFVBCs bc_type, std::set<BoundaryID> & bnd_ids); //Query for \p INSFVBCs on \p bc_id and add if query successful
};
