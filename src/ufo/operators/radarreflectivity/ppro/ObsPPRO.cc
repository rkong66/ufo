/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/radarreflectivity/ppro/ObsPPRO.h"

#include <ostream>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsPPRO> makerPPRO_("PPRO");
// -----------------------------------------------------------------------------

ObsPPRO::ObsPPRO(const ioda::ObsSpace & odb,
                       const ObsPPROParameters & params)
  : ObsOperatorBase(odb), keyOperPPRO_(0), odb_(odb), varin_()
{
  // dualpol operator requires the column indices of assimvariables (simulated variables), rather than the 
  // column indices of variables in ioda *.h5 data (obsvariables)
  ufo_PPRO_setup_f90(keyOperPPRO_, params.toConfiguration(),
                                    odb.assimvariables(), varin_);

  oops::Log::trace() << "ObsPPRO created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsPPRO::~ObsPPRO() {
  ufo_PPRO_delete_f90(keyOperPPRO_);
  oops::Log::trace() << "ObsPPRO destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsPPRO::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                         ObsDiagnostics &,
                                         const QCFlags_t & qc_flags) const {
  ufo_PPRO_simobs_f90(keyOperPPRO_, gv.toFortran(), odb_, ovec.nvars(), ovec.nlocs(),
                         ovec.toFortran());
  oops::Log::trace() << "ObsPPRO: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsPPRO::print(std::ostream & os) const {
  os << "ObsPPRO::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
