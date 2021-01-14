/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_OBSFILTERS_H_
#define TEST_UFO_OBSFILTERS_H_

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/ObsFilters.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/interface/ObsTestsFixture.h"
#include "test/TestEnvironment.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsTraits.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit
{
  // Don't use the contracted output for these types: the current implementation works only
  // with integer types.
  // TODO(wsmigaj) Report this (especially for floats) as a bug in eckit?
  template <> struct VectorPrintSelector<float> { typedef VectorPrintSimple selector; };
  template <> struct VectorPrintSelector<util::DateTime> { typedef VectorPrintSimple selector; };
  template <> struct VectorPrintSelector<util::Duration> { typedef VectorPrintSimple selector; };
}  // namespace eckit

namespace ufo {
namespace test {

// -----------------------------------------------------------------------------

/// \brief Options used to configure comparison of a variable generated by a filter
/// against a reference variable loaded from the input IODA file.
class CompareVariablesParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CompareVariablesParameters, Parameters)

 public:
  /// Variable that should be contained against reference values.
  oops::RequiredParameter<ufo::Variable> test{"test", this};
  /// Variable containing the reference values.
  oops::RequiredParameter<ufo::Variable> reference{"reference", this};

  /// If set, the comparison will succeed if all corresponding elements of the test and reference
  /// variables differ by at most `absTol`. Otherwise the comparison will succeed only if all
  /// corresponding elements match exactly.
  oops::OptionalParameter<float> absTol{"absTol", this};
};

// -----------------------------------------------------------------------------

/// \brief Options used to configure a test running a sequence of filters on observations
/// from a single obs space.
///
/// Note: at least one of the options whose names end in `Benchmark` or the `compareVariables`
/// option needs to be set (otherwise the test won't test much).
class ObsTypeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsTypeParameters, Parameters)

 public:
  /// Options used to configure the observation space.
  oops::Parameter<eckit::LocalConfiguration> obsSpace{
    "obs space", eckit::LocalConfiguration(), this};

  /// Options used to configure observation filters.
  oops::Parameter<std::vector<oops::ObsFilterParametersWrapper<ObsTraits>>> obsFilters{
    "obs filters", {}, this};

  /// Options passed to the observation operator that will be applied during the test. If not set,
  /// no observation operator will be applied. To speed up tests of filters that depend on the
  /// values produced by the observation operator (model equivalents), these values can be
  /// precalculated and stored in the IODA file used to initialize the ObsSpace. In that case the
  /// `obs operator` keyword should be omitted and instead the `HofX` option should be set to the
  /// name of the group of ObsSpace variables containing the precalculated model equivalents.
  oops::OptionalParameter<eckit::LocalConfiguration> obsOperator{"obs operator", this};

  /// Group of variables storing precalculated model equivalents of observations. See the
  /// description of the `obs operator` option for more information.
  oops::OptionalParameter<std::string> hofx{"HofX", this};

  /// Options used to load GeoVaLs from a file. Required if any observation filters depend on
  /// GeoVaLs or of the `obs operator` option is set.
  oops::OptionalParameter<eckit::LocalConfiguration> geovals{"geovals", this};

  /// Options used to load observation diagnostics from a file. Required if any observation filters
  /// depend on observation diagnostics.
  oops::OptionalParameter<eckit::LocalConfiguration> obsDiagnostics{"obs diagnostics", this};

  /// Options used to configure the observation bias.
  oops::OptionalParameter<eckit::LocalConfiguration> obsBias{"obs bias", this};

  /// Indices of observations expected to pass quality control.
  ///
  /// The observations are numbered as in the input IODA file.
  oops::OptionalParameter<std::vector<size_t>> passedObservationsBenchmark{
    "passedObservationsBenchmark", this};
  /// Number of observations expected to pass quality control.
  oops::OptionalParameter<size_t> passedBenchmark{"passedBenchmark", this};

  /// Indices of observations expected to fail quality control.
  ///
  /// The observations are numbered as in the input IODA file.
  oops::OptionalParameter<std::vector<size_t>> failedObservationsBenchmark{
    "failedObservationsBenchmark", this};
  /// Number of observations expected to fail quality control.
  oops::OptionalParameter<size_t> failedBenchmark{"failedBenchmark", this};

  /// An integer corresponding to one of the constants in the QCflags namespace.
  oops::OptionalParameter<int> benchmarkFlag{"benchmarkFlag", this};

  /// Indices of observations expected to receive the `benchmarkFlag` flag.
  ///
  /// The observations are numbered as in the input IODA file.
  oops::OptionalParameter<std::vector<size_t>> flaggedObservationsBenchmark{
    "flaggedObservationsBenchmark", this};
  /// Number of observations expected to receive the `benchmarkFlag` flag.
  oops::OptionalParameter<size_t> flaggedBenchmark{"flaggedBenchmark", this};

  /// A list of options indicating variables whose final values should be compared against reference
  /// values loaded from the input IODA file.
  oops::OptionalParameter<std::vector<CompareVariablesParameters>> compareVariables{
    "compareVariables", this};
};

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the ObsFilters test.
class ObsFiltersParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsFiltersParameters, Parameters)

 public:
  /// Only observations taken at times lying in the (`window begin`, `window end`] interval
  /// will be included in observation spaces.
  oops::RequiredParameter<util::DateTime> windowBegin{"window begin", this};
  oops::RequiredParameter<util::DateTime> windowEnd{"window end", this};
  /// A list whose elements are used to configure tests running sequences of filters on
  /// observations from individual observation spaces.
  oops::Parameter<std::vector<ObsTypeParameters>> observations{"observations", {}, this};
};

// -----------------------------------------------------------------------------

//!
//! \brief Convert indices of observations held by this process to global observation indices.
//!
//! \param[inout] indices
//!   On input: local indices of observations held by this process. On output: corresponding
//!   global observation indices.
//! \param globalIdxFromLocalIdx
//!   A vector whose ith element is the global index of ith observation held by this process.
//!
void convertLocalObsIndicesToGlobal(std::vector<size_t> &indices,
                                    const std::vector<size_t> &globalIdxFromLocalIdx) {
  for (size_t &index : indices)
    index = globalIdxFromLocalIdx[index];
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations whose quality control flags satisfy the
//! predicate in at least one variable.
//!
//! \param qcFlags
//!   Vector of quality control flags for all observations.
//! \param comm
//!   The MPI communicator used by the ObsSpace.
//! \param globalIdxFromLocalIdx
//!   A vector whose ith element is the global index of ith observation held by this process.
//! \param predicate
//!   A function object taking an argument of type int and returning bool.
//!
template <typename Predicate>
std::vector<size_t> getObservationIndicesWhere(
    const ObsTraits::ObsDataVector<int> &qcFlags,
    const ioda::Distribution &obsDistribution,
    const std::vector<size_t> &globalIdxFromLocalIdx,
    const Predicate &predicate) {
  std::vector<size_t> indices;
  for (size_t locIndex = 0; locIndex < qcFlags.nlocs(); ++locIndex) {
    bool satisfied = false;
    for (size_t varIndex = 0; varIndex < qcFlags.nvars(); ++varIndex) {
      if (predicate(qcFlags[varIndex][locIndex])) {
        satisfied = true;
        break;
      }
    }
    if (satisfied) {
      indices.push_back(locIndex);
    }
  }

  convertLocalObsIndicesToGlobal(indices, globalIdxFromLocalIdx);
  obsDistribution.allGatherv(indices);
  std::sort(indices.begin(), indices.end());
  return indices;
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations that have passed quality control in
//! at least one variable.
//!
std::vector<size_t> getPassedObservationIndices(const ObsTraits::ObsDataVector<int> &qcFlags,
                                                const ioda::Distribution &obsDistribution,
                                                const std::vector<size_t> &globalIdxFromLocalIdx) {
  return getObservationIndicesWhere(qcFlags, obsDistribution, globalIdxFromLocalIdx,
                                    [](int qcFlag) { return qcFlag == 0; });
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations that have failed quality control in
//! at least one variable.
//!
std::vector<size_t> getFailedObservationIndices(const ObsTraits::ObsDataVector<int> &qcFlags,
                                                const ioda::Distribution &obsDistribution,
                                                const std::vector<size_t> &globalIdxFromLocalIdx) {
  return getObservationIndicesWhere(qcFlags, obsDistribution, globalIdxFromLocalIdx,
                                    [](int qcFlag) { return qcFlag != 0; });
}

// -----------------------------------------------------------------------------

//!
//! Return the indices of observations whose quality control flag is set to \p flag in
//! at least one variable.
//!
std::vector<size_t> getFlaggedObservationIndices(const ObsTraits::ObsDataVector<int> &qcFlags,
                                                 const ioda::Distribution &obsDistribution,
                                                 const std::vector<size_t> &globalIdxFromLocalIdx,
                                                 int flag) {
  return getObservationIndicesWhere(qcFlags, obsDistribution, globalIdxFromLocalIdx,
                                    [flag](int qcFlag) { return qcFlag == flag; });
}

// -----------------------------------------------------------------------------

//!
//! Return the number of elements of \p data with at least one nonzero component.
//!
size_t numNonzero(const ObsTraits::ObsDataVector<int> & data) {
  size_t result = 0;
  for (size_t locIndex = 0; locIndex < data.nlocs(); ++locIndex) {
    for (size_t varIndex = 0; varIndex < data.nvars(); ++varIndex) {
      if (data[varIndex][locIndex] != 0)
        ++result;
    }
  }
  return result;
}

// -----------------------------------------------------------------------------

//!
//! Return the number of elements of \p data with at least one component equal to \p value.
//!
size_t numEqualTo(const ObsTraits::ObsDataVector<int> & data, int value) {
  size_t result = 0;
  for (size_t locIndex = 0; locIndex < data.nlocs(); ++locIndex) {
    for (size_t varIndex = 0; varIndex < data.nvars(); ++varIndex) {
      if (data[varIndex][locIndex] == value)
        ++result;
    }
  }
  return result;
}

// -----------------------------------------------------------------------------

template <typename T>
void expectVariablesEqual(const ObsTraits::ObsSpace &obsspace,
                          const ufo::Variable &referenceVariable,
                          const ufo::Variable &testVariable)
{
  std::vector<T> reference(obsspace.nlocs());
  obsspace.get_db(referenceVariable.group(), referenceVariable.variable(), reference);
  std::vector<T> test(obsspace.nlocs());
  obsspace.get_db(testVariable.group(), testVariable.variable(), test);
  EXPECT_EQUAL(reference, test);
}

// -----------------------------------------------------------------------------

void expectVariablesApproximatelyEqual(const ObsTraits::ObsSpace &obsspace,
                                       const ufo::Variable &referenceVariable,
                                       const ufo::Variable &testVariable,
                                       float absTol)
{
  std::vector<float> reference(obsspace.nlocs());
  obsspace.get_db(referenceVariable.group(), referenceVariable.variable(), reference);
  std::vector<float> test(obsspace.nlocs());
  obsspace.get_db(testVariable.group(), testVariable.variable(), test);
  EXPECT(oops::are_all_close_absolute(reference, test, absTol));
}

// -----------------------------------------------------------------------------

void testFilters() {
  typedef ::test::ObsTestsFixture<ObsTraits> Test_;
  typedef oops::GeoVaLs<ufo::ObsTraits>           GeoVaLs_;
  typedef oops::ObsDiagnostics<ufo::ObsTraits>    ObsDiags_;
  typedef oops::ObsAuxControl<ufo::ObsTraits>     ObsAuxCtrl_;
  typedef oops::ObsFilters<ufo::ObsTraits>        ObsFilters_;
  typedef oops::ObsOperator<ufo::ObsTraits>       ObsOperator_;
  typedef oops::ObsVector<ufo::ObsTraits>         ObsVector_;
  typedef oops::ObsSpace<ufo::ObsTraits>          ObsSpace_;

  ObsFiltersParameters params;
  params.validateAndDeserialize(::test::TestEnvironment::config());

  for (std::size_t jj = 0; jj < Test_::obspace().size(); ++jj) {
/// identify parameters used for this group of observations
    const ObsTypeParameters &typeParams = params.observations.value()[jj];

/// init QC and error
    std::shared_ptr<oops::ObsDataVector<ufo::ObsTraits, float> > obserr
      (new oops::ObsDataVector<ufo::ObsTraits, float>(Test_::obspace()[jj],
               Test_::obspace()[jj].obsvariables(), "ObsError"));
    std::shared_ptr<oops::ObsDataVector<ufo::ObsTraits, int> >
      qcflags(new oops::ObsDataVector<ufo::ObsTraits, int>  (Test_::obspace()[jj],
               Test_::obspace()[jj].obsvariables()));

//  Create filters and run preProcess
    ObsFilters_ filters(Test_::obspace()[jj], typeParams.obsFilters, qcflags, obserr);
    filters.preProcess();

/// call priorFilter and postFilter if hofx is available
    oops::Variables geovars = filters.requiredVars();
    oops::Variables diagvars = filters.requiredHdiagnostics();
    if (typeParams.hofx.value() != boost::none) {
///   read GeoVaLs from file if required
      std::unique_ptr<const GeoVaLs_> gval;
      if (geovars.size() > 0) {
        if (typeParams.geovals.value() == boost::none)
          throw eckit::UserError(std::to_string(jj) + "th element of the 'observations' list "
                                                      "requires a 'geovals' section", Here());
        gval.reset(new GeoVaLs_(*typeParams.geovals.value(), Test_::obspace()[jj], geovars));
        filters.priorFilter(*gval);
      } else {
        oops::Log::info() << "Filters don't require geovals, priorFilter not called" << std::endl;
      }
///   read H(x) and ObsDiags from file
      oops::Log::info() << "HofX section specified, reading HofX from file" << std::endl;
      const std::string &hofxgroup = *typeParams.hofx.value();
      ObsVector_ hofx(Test_::obspace()[jj], hofxgroup);
      eckit::LocalConfiguration obsdiagconf;
      if (diagvars.size() > 0) {
        if (typeParams.obsDiagnostics.value() == boost::none)
          throw eckit::UserError(std::to_string(jj) + "th element of the 'observations' list "
                                                      "requires an 'obs diagnostics' section",
                                 Here());
        obsdiagconf = *typeParams.obsDiagnostics.value();
        oops::Log::info() << "Obs diagnostics section specified, reading obs diagnostics from file"
                          << std::endl;
      }
      const ObsDiags_ diags(obsdiagconf, Test_::obspace()[jj], diagvars);
      filters.postFilter(hofx, diags);
    } else if (typeParams.obsOperator.value() != boost::none) {
///   read GeoVaLs, compute H(x) and ObsDiags
      oops::Log::info() << "ObsOperator section specified, computing HofX" << std::endl;
      ObsOperator_ hop(Test_::obspace()[jj], *typeParams.obsOperator.value());
      // TODO(wsmigaj): ObsAuxCtrl currently expects to receive the top-level configuration
      // even though the implementations I've seen only reference elements of the "obs bias"
      // subconfiguration. If that's universally true, ObsAuxCtrl should be refactored so that
      // it can be passed just the contents of the "obs bias" section.
      const ObsAuxCtrl_ ybias(Test_::obspace()[jj], typeParams.toConfiguration());
      ObsVector_ hofx(Test_::obspace()[jj]);
      oops::Variables vars;
      vars += hop.requiredVars();
      vars += filters.requiredVars();
      if (typeParams.obsBias.value() != boost::none) vars += ybias.requiredVars();
      if (typeParams.geovals.value() == boost::none)
        throw eckit::UserError(std::to_string(jj) + "th element of the 'observations' list "
                                                    "requires a 'geovals' section", Here());
      const GeoVaLs_ gval(*typeParams.geovals.value(), Test_::obspace()[jj], vars);
      oops::Variables diagvars;
      diagvars += filters.requiredHdiagnostics();
      if (typeParams.obsBias.value() != boost::none) diagvars += ybias.requiredHdiagnostics();
      ObsDiags_ diags(Test_::obspace()[jj], hop.locations(), diagvars);
      filters.priorFilter(gval);
      hop.simulateObs(gval, hofx, ybias, diags);
      hofx.save("hofx");
      filters.postFilter(hofx, diags);
    } else if (geovars.size() > 0) {
///   Only call priorFilter
      if (typeParams.geovals.value() == boost::none)
        throw eckit::UserError(std::to_string(jj) + "th element of the 'observations' list "
                                                    "requires a 'geovals' section", Here());
      const GeoVaLs_ gval(*typeParams.geovals.value(), Test_::obspace()[jj], geovars);
      filters.priorFilter(gval);
      oops::Log::info() << "HofX or ObsOperator sections not provided for filters, " <<
                           "postFilter not called" << std::endl;
    } else {
///   no need to run priorFilter or postFilter
      oops::Log::info() << "GeoVaLs not required, HofX or ObsOperator sections not " <<
                           "provided for filters, only preProcess was called" << std::endl;
    }

    qcflags->save("EffectiveQC");
    const std::string errname = "EffectiveError";
    obserr->save(errname);

//  Compare with known results
    bool atLeastOneBenchmarkFound = false;
    const ObsTraits::ObsSpace &obsspace = Test_::obspace()[jj].obsspace();

    if (typeParams.passedObservationsBenchmark.value() != boost::none) {
      atLeastOneBenchmarkFound = true;
      const std::vector<size_t> &passedObsBenchmark =
          *typeParams.passedObservationsBenchmark.value();
      const std::vector<size_t> passedObs = getPassedObservationIndices(
            qcflags->obsdatavector(), obsspace.distribution(), obsspace.index());
      EXPECT_EQUAL(passedObs, passedObsBenchmark);
    }

    if (typeParams.passedBenchmark.value() != boost::none) {
      atLeastOneBenchmarkFound = true;
      const size_t passedBenchmark = *typeParams.passedBenchmark.value();
      size_t passed = numEqualTo(qcflags->obsdatavector(), ufo::QCflags::pass);
      obsspace.distribution().sum(passed);
      EXPECT_EQUAL(passed, passedBenchmark);
    }

    if (typeParams.failedObservationsBenchmark.value() != boost::none) {
      atLeastOneBenchmarkFound = true;
      const std::vector<size_t> &failedObsBenchmark =
          *typeParams.failedObservationsBenchmark.value();
      const std::vector<size_t> failedObs = getFailedObservationIndices(
            qcflags->obsdatavector(), obsspace.distribution(), obsspace.index());
      EXPECT_EQUAL(failedObs, failedObsBenchmark);
    }

    if (typeParams.failedBenchmark.value() != boost::none) {
      atLeastOneBenchmarkFound = true;
      const size_t failedBenchmark = *typeParams.failedBenchmark.value();
      size_t failed = numNonzero(qcflags->obsdatavector());
      obsspace.distribution().sum(failed);
      EXPECT_EQUAL(failed, failedBenchmark);
    }

    if (typeParams.benchmarkFlag.value() != boost::none) {
      const int flag = *typeParams.benchmarkFlag.value();

      if (typeParams.flaggedObservationsBenchmark.value() != boost::none) {
        atLeastOneBenchmarkFound = true;
        const std::vector<size_t> &flaggedObsBenchmark =
            *typeParams.flaggedObservationsBenchmark.value();
        const std::vector<size_t> flaggedObs =
            getFlaggedObservationIndices(qcflags->obsdatavector(), obsspace.distribution(),
                                         obsspace.index(), flag);
        EXPECT_EQUAL(flaggedObs, flaggedObsBenchmark);
      }

      if (typeParams.flaggedBenchmark.value() != boost::none) {
        atLeastOneBenchmarkFound = true;
        const size_t flaggedBenchmark = *typeParams.flaggedBenchmark.value();
        size_t flagged = numEqualTo(qcflags->obsdatavector(), flag);
        obsspace.distribution().sum(flagged);
        EXPECT_EQUAL(flagged, flaggedBenchmark);
      }
    }

    if (typeParams.compareVariables.value() != boost::none) {
      for (const CompareVariablesParameters &compareVariablesParams :
           *typeParams.compareVariables.value()) {
        atLeastOneBenchmarkFound = true;

        const ufo::Variable &referenceVariable = compareVariablesParams.reference;
        const ufo::Variable &testVariable = compareVariablesParams.test;

        switch (obsspace.dtype(referenceVariable.group(), referenceVariable.variable())) {
        case ioda::ObsDtype::Integer:
          expectVariablesEqual<int>(obsspace, referenceVariable, testVariable);
          break;
        case ioda::ObsDtype::String:
          expectVariablesEqual<std::string>(obsspace, referenceVariable, testVariable);
          break;
        case ioda::ObsDtype::DateTime:
          expectVariablesEqual<util::DateTime>(obsspace, referenceVariable, testVariable);
          break;
        case ioda::ObsDtype::Float:
          if (compareVariablesParams.absTol.value() == boost::none) {
            expectVariablesEqual<float>(obsspace, referenceVariable, testVariable);
          } else {
            const float tol = *compareVariablesParams.absTol.value();
            expectVariablesApproximatelyEqual(obsspace, referenceVariable, testVariable, tol);
          }
          break;
        case ioda::ObsDtype::None:
          ASSERT_MSG(false, "Reference variable not found in observation space");
        }
      }
    }

    EXPECT(atLeastOneBenchmarkFound);
  }
}

// -----------------------------------------------------------------------------

class ObsFilters : public oops::Test {
  typedef ::test::ObsTestsFixture<ObsTraits> Test_;
 public:
  ObsFilters() {}
  virtual ~ObsFilters() {}
 private:
  std::string testid() const override {return "test::ObsFilters";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("ufo/ObsFilters/testFilters")
      { testFilters(); });
  }

  void clear() const override {
    Test_::reset();
  }
};

// =============================================================================

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_OBSFILTERS_H_
