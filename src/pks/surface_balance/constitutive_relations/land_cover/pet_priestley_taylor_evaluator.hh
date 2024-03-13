/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

//! Evaluates potential evapotranpiration (PET) using Priestley & Taylor formulation.
/*!

This implementation is based on models provided in the PRMS-IV, Version 4, see
pages 90-93, Equations 1-57 to 1-60

Requires the following dependencies:

.. _pet-priestley-taylor-evaluator-spec:
.. admonition:: pet-priestley-taylor-evaluator-spec:

   * `"include limiter`" ``[bool]`` **false** If true, multiply potential ET by
     a limiter to get an actual ET.
   * `"limiter number of dofs`" ``[int]`` **1** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides how many.
   * `"limiter dof`" ``[int]`` **0** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides which one to use.
   * `"include 1 - limiter`" ``[bool]`` **false** If true, multiply potential
     ET by 1 - a limiter (e.g. a limiter that partitions between two pools) to
     get actual ET.
   * `"1 - limiter number of dofs`" ``[int]`` **1** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides how many.
   * `"1 - limiter dof`" ``[int]`` **0** Area fractions are often used
     as limiters, and these have multiple dofs.  This provides which one to use.
   * `"sublimate snow`" ``[bool]`` **false** If true, use latent heat of
      vaporization of snow, not water.

   KEYS:

   - `"air temperature`" **DOMAIN-air_temperature** Air temp, in [K]
   - `"surface temperature`" **DOMAIN-temperature** Ground or leaf temp, in [K].  Note this may be the
      same as air temperature.
   - `"elevation`" **DOMAIN-elevation** Elevation [m]
   - `"net radiation`" **DOMAIN-net_radiation** [W m^-2] Net radiation balance, positive to the ground.
   - `"limiter`" [-] See `"include limiter`" above.
   - `"1 - limiter`" [-] See `"include 1 - limiter`" above.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {




class PETPriestleyTaylorEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit PETPriestleyTaylorEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  PETPriestleyTaylorEvaluator(const PETPriestleyTaylorEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new PETPriestleyTaylorEvaluator(*this));
  }

  std::string getType() const override { return eval_type; }

protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  static const std::string eval_type;

 protected:
  Key domain_;
  Key evap_type_;
  Key air_temp_key_;
  Key surf_temp_key_;
  Key elev_key_;
  Key rad_key_;
  Key limiter_key_;
  Key one_minus_limiter_key_;

  double pt_alpha_;
  bool limiter_, one_minus_limiter_;
  int limiter_nvecs_, one_minus_limiter_nvecs_;
  int limiter_dof_, one_minus_limiter_dof_;
  bool compatible_;

  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, PETPriestleyTaylorEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
