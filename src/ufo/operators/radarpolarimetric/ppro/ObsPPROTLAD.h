/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_OBSPPROTLAD_H_
#define UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_OBSPPROTLAD_H_


#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/radarpolarimetric/ppro/ObsPPRO.h"
#include "ufo/operators/radarpolarimetric/ppro/ObsPPROTLAD.interface.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// PPRO observation operator
class ObsPPROTLAD : public LinearObsOperatorBase,
                          private util::ObjectCounter<ObsPPROTLAD> {
 public:
  typedef ObsPPROParameters Parameters_;

  static const std::string classname() {return "ufo::ObsPPROTLAD";}

  ObsPPROTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsPPROTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const QCFlags_t &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, const QCFlags_t &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperPPRO_;}
  const int & toFortran() const {return keyOperPPRO_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperPPRO_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARPOLARIMETRIC_PPRO_OBSPPROTLAD_H_
