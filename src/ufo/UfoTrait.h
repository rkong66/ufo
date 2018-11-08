/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UFOTRAIT_H_
#define UFO_UFOTRAIT_H_

#include <string>

#include "ioda/Locations.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "GeoVaLs.h"
#include "ObsBias.h"
#include "ObsBiasCovariance.h"
#include "ObsBiasIncrement.h"
#include "ObsOperator.h"
#include "LinearObsOperator.h"

namespace ufo {

struct UfoTrait {
  static std::string name() {return "UFO";}

  typedef ufo::GeoVaLs             GeoVaLs;
  typedef ioda::Locations          Locations;
  typedef ioda::ObsSpace           ObsSpace;
  typedef ioda::ObsVector          ObsVector;

  typedef ufo::ObsOperator         ObsOperator;
  typedef ufo::LinearObsOperator   LinearObsOperator;

  typedef ufo::ObsBias             ObsAuxControl;
  typedef ufo::ObsBiasIncrement    ObsAuxIncrement;
  typedef ufo::ObsBiasCovariance   ObsAuxCovariance;
};

}  // namespace ufo

#endif  // UFO_UFOTRAIT_H_
