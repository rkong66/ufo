/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_OBSPPRO_H_
#define UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_OBSPPRO_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/operators/radarpolarimetric/ppro/ObsPPRO.interface.h"
#include "ufo/operators/radarshared/MicrophysicsOptions.h"
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
 
  oops::Parameter<double> coeff_melt{
    "tuning coefficient for melting",                // YAML name
    "Coefficient scaling the geometric mean sqrt(qr*qx)/sqrt(nr*qx) in Liu et al. (2024) melting model", // description
    0.3,                                             // default
    this
  };

  oops::Parameter<std::string> melting_scheme{
    "melting scheme",                                // YAML name
    "Melting scheme option: 'liu24' (sqrt formula for all), 'zhang24' (original)", // description
    "liu24",                                         // default: liu24
    this
  };

  oops::Parameter<bool> enable_melting_transition{
    "enable melting transition",                     // YAML name
    "Enable smooth transition function for melting based on qx/qr ratio balance", // description
    false,                                           // default: off
    this
  };

  oops::Parameter<double> snow_ratio_low{
    "snow ratio low",
    "Lower bound of transition for snow (factor=0 below this)",
    0.05,
    this
  };

  oops::Parameter<double> snow_ratio_high{
    "snow ratio high",
    "Upper bound of transition for snow (factor=1 above this)",
    0.2,
    this
  };

  oops::Parameter<double> graupel_ratio_low{
    "graupel ratio low",
    "Lower bound of transition for graupel (factor=0 below this)",
    0.05,
    this
  };

  oops::Parameter<double> graupel_ratio_high{
    "graupel ratio high",
    "Upper bound of transition for graupel (factor=1 above this)",
    0.2,
    this
  };

  oops::Parameter<double> hail_ratio_low{
    "hail ratio low",
    "Lower bound of transition for hail (factor=0 below this)",
    0.05,
    this
  };

  oops::Parameter<double> hail_ratio_high{
    "hail ratio high",
    "Upper bound of transition for hail (factor=1 above this)",
    0.2,
    this
  };

  // Melting water content limit parameters
  oops::Parameter<bool> enable_melting_water_limit{
    "enable melting water limit",
    "Enable limiting melting water content (qmsr/qmgr/qmhr) to a fraction of qr. When false (default), no limit",
    false,
    this
  };

  oops::Parameter<double> melting_water_fraction{
    "melting water fraction",
    "Fraction of qr to limit melting water content (only used when enable_melting_water_limit is true)",
    0.3,
    this
  };

  // dmmax configuration: maximum mean diameter limits for hydrometeor species [mm]
  // Each species can be individually configured. If not specified, uses default.
  oops::Parameter<double> dmmax_rain{
    "dmmax rain",
    "Max mean diameter for rain [mm]",
    5.0,
    this
  };

  oops::Parameter<double> dmmax_pure_snow{
    "dmmax pure snow",
    "Max mean diameter for pure snow [mm]",
    10.0,
    this
  };

  oops::Parameter<double> dmmax_melting_snow{
    "dmmax melting snow",
    "Max mean diameter for melting snow [mm]",
    5.0,
    this
  };

  oops::Parameter<double> dmmax_pure_graupel{
    "dmmax pure graupel",
    "Max mean diameter for pure graupel [mm]",
    10.0,
    this
  };

  oops::Parameter<double> dmmax_melting_graupel{
    "dmmax melting graupel",
    "Max mean diameter for melting graupel [mm]",
    5.0,
    this
  };

  oops::Parameter<double> dmmax_pure_hail{
    "dmmax pure hail",
    "Max mean diameter for pure hail [mm]",
    10.0,
    this
  };

  oops::Parameter<double> dmmax_melting_hail{
    "dmmax melting hail",
    "Max mean diameter for melting hail [mm]",
    5.0,
    this
  };
 
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
#endif  // UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_OBSPPRO_H_
