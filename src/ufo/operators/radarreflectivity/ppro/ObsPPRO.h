/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_RADARREFLECTIVITY_PPRO_OBSPPRO_H_
#define UFO_OPERATORS_RADARREFLECTIVITY_PPRO_OBSPPRO_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/operators/radarreflectivity/ppro/ObsPPRO.interface.h"
#include "ufo/operators/radarreflectivity/MicrophysicsOptions.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

  // PPRO-specific operator type enumeration
  enum class PolarimetricOperatorOption {
    ZHANG21, TCWA2
  };

  struct PolarimetricOperatorOptionParameterTraitsHelper {
    typedef PolarimetricOperatorOption EnumType;
    static constexpr char enumTypeName[] = "PolarimetricOperatorOption";
    static constexpr util::NamedEnumerator<PolarimetricOperatorOption> namedValues[] = {
      { PolarimetricOperatorOption::ZHANG21, "Zhang21" },
      { PolarimetricOperatorOption::TCWA2, "TCWA2" }
    };
  };
} // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::PolarimetricOperatorOption> :
    public EnumParameterTraits<ufo::PolarimetricOperatorOptionParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

// -----------------------------------------------------------------------------
class ObsPPROParameters : public ObsOperatorParametersBase  {
  OOPS_CONCRETE_PARAMETERS(ObsPPROParameters, ObsOperatorParametersBase )
 public:
  /// Vertical Coordinate
  oops::Parameter<std::string> VertCoord
    {"VertCoord", 
     "vertical coordinate used by the model",
     "geometric_height",  // this should be consistent with var_prs defined in ufo_vars_mod
     this};
    
  oops::Parameter<MicrophysicsOption> microOption
    {"microphysics option",
     "microphysics option by name (e.g. Thompson)",
     MicrophysicsOption::THOMPSON,
     this};

  oops::Parameter<PolarimetricOperatorOption> pproOption
    {"polarimetric operator",
     "polarimetric operator by name (Zhang21 or TCWA2)",
     PolarimetricOperatorOption::ZHANG21,
     this};

  oops::Parameter<bool> use_variational
    {"use variational method",
     "use variational TL/AD method (P-PRO)",
     false,
     this};
  
  /*
   * The following list of hydrometeor species mixing ratios and number concentrations
   * need to be consistent with list of variables in ufo_variables_mod, which has some
   * confusing names that need to be resolved in the future - all to match CCPP
   * convention.  In the meantime, FV3 and MPAS use different variables.  The code
   * will default to FV3 names (rain_water, snow_water, and graupel), but all input
   * species can be overridden with these optional YAML parameters, which for MPAS
   * is mass_content_of_rain_in_atmosphere_layer for example.
   */

  oops::Parameter<std::string> var_rain_mixing_ratio
    {"var_rain_mixing_ratio",
     "Name of model rain mixing ratio variable",
     "rain_water",  // this should be consistent with var_qr
     this};

  oops::Parameter<std::string> var_snow_mixing_ratio
    {"var_snow_mixing_ratio",
     "Name of model snow mixing ratio variable",
     "snow_water",  // this should be consistent with var_qs
     this};

  oops::Parameter<std::string> var_graupel_mixing_ratio
    {"var_graupel_mixing_ratio",
     "Name of model graupel mixing ratio variable",
     "graupel",  // this should be consistent with var_qg
     this};

  oops::Parameter<std::string> var_hail_mixing_ratio
    {"var_hail_mixing_ratio",
     "Name of model hail mixing ratio variable",
     "hail",  // this should be consistent with var_qh
     this};

  oops::Parameter<std::string> var_rain_number_concentration
    {"var_rain_number_concentration",
     "Name of model rain number concentration variable",
     "rain_number_concentration",  // this should be consistent with var_nr
     this};

  oops::Parameter<std::string> var_snow_number_concentration
    {"var_snow_number_concentration",
     "Name of model snow number concentration variable",
     "snow_number_concentration",  // this should be consistent with var_ns
     this};

  oops::Parameter<std::string> var_graupel_number_concentration
    {"var_graupel_number_concentration",
     "Name of model graupel number concentration variable",
     "graupel_number_concentration",  // this should be consistent with var_ng
     this};

  oops::Parameter<std::string> var_hail_number_concentration
    {"var_hail_number_concentration",
     "Name of model hail number concentration variable",
     "graupel_hail_concentration",  // this should be consistent with var_nh
     this};

  oops::Parameter<std::string> var_graupel_vol_mixing_ratio
    {"var_graupel_vol_mixing_ratio",
     "Name of model graupel volume mixing ratio variable",
     "volume_mixing_ratio_of_graupel_in_air",  // this should be consistent with var_qvg
     this};

  oops::Parameter<std::string> var_hail_vol_mixing_ratio
    {"var_hail_vol_mixing_ratio",
     "Name of model hail volume mixing ratio variable",
     "volume_mixing_ratio_of_hail_in_air",  // this should be consistent with var_qvh
     this};
  
};

/// RadarReflectivity observation operator class
class ObsPPRO : public ObsOperatorBase,
                   private util::ObjectCounter<ObsPPRO> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  typedef ObsPPROParameters Parameters_;
  static const std::string classname() {return "ufo::ObsPPRO";}

  ObsPPRO(const ioda::ObsSpace &, const ObsPPROParameters &);
  virtual ~ObsPPRO();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                  const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperPPRO_;}
  const int & toFortran() const {return keyOperPPRO_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperPPRO_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARREFLECTIVITY_PPRO_OBSPPRO_H_
