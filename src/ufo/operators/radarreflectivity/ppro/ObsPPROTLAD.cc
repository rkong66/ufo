/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/radarreflectivity/ppro/ObsPPROTLAD.h"

#include <ostream>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsPPROTLAD>
                makerPPROTL_("PPRO");
// -----------------------------------------------------------------------------

ObsPPROTLAD::ObsPPROTLAD(const ioda::ObsSpace & odb,
                                           const Parameters_ & params)
  : LinearObsOperatorBase(odb), keyOperPPRO_(0), varin_()
{
  // dualpol operator requires the column indices of assimvariables (simulated variables), rather than the 
  // column indices of variables in ioda *.h5 data (obsvariables)
  ufo_PPRO_tlad_setup_f90(keyOperPPRO_, params.toConfiguration(), 
                                    odb.assimvariables(), varin_);

  oops::Log::trace() << "ObsPPROTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsPPROTLAD::~ObsPPROTLAD() {
  ufo_PPRO_tlad_delete_f90(keyOperPPRO_);
  oops::Log::trace() << "ObsPPROTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsPPROTLAD::setTrajectory(const GeoVaLs & geovals, ObsDiagnostics &,
                                const QCFlags_t & qc_flags_t) {
  oops::Log::trace() << "ObsPPROTLAD::setTrajectory entering" << std::endl;

  ufo_PPRO_tlad_settraj_f90(keyOperPPRO_, geovals.toFortran(),
                                           obsspace());

  oops::Log::trace() << "ObsPPROTLAD::setTrajectory exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsPPROTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                                const QCFlags_t & qc_flags_t) const {
  ufo_PPRO_simobs_tl_f90(keyOperPPRO_, geovals.toFortran(),
                                        obsspace(), ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsPPROTLAD::simulateObsTL exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsPPROTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                                const QCFlags_t & qc_flags_t) const {
  ufo_PPRO_simobs_ad_f90(keyOperPPRO_, geovals.toFortran(),
                                        obsspace(), ovec.nvars(), ovec.nlocs(), ovec.toFortran());

  oops::Log::trace() << "ObsPPROTLAD::simulateObsAD exiting" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsPPROTLAD::print(std::ostream & os) const {
  os << "ObsPPROTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
