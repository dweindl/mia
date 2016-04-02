#ifndef CONFIG_H
#define CONFIG_H

#define STRINGIFY(str) # str
#define COMPOUND_GROUPING_FEATURE STRINGIFY(CMP_ID)

namespace mia {

// check settings.cpp


static const bool CMP_ID_USE_RI = true; /** use retention index */
static const int CMP_ID_RI_TOL = 100; /** tolerance for retention index difference */
static const double CMP_ID_SCORE_CUTOFF = 0.75; /** Cutoff for spectrum matching score */

static const double CMP_MATCHING_RI_TOL = 5; /** RI tolerance for peak matching different chromatograms */
//TODO: to config
static const double CMP_MATCHING_SCORE_CUTOFF = 0.85; /** Spec score for peak matching different chromatograms */

static const double GCMS_PURE_FACTOR = 0.5; /** ... */
static const double GCMS_IMPURE_FACTOR = 0.5; /** ... */
static const std::string CMP_ID_LIBRARY = ""; /** MD library file */

static const int LABELS_MAX_HITS = 1;
static const bool LID_FILTER_BY_CONF_INTERVAL = true; // NTFD: hardcoded
static const double LID_MIN_SIGNAL_TO_NOISE = 5; // NTFD: hardcoded
static const double LID_REQUIRED_SPEC_FREQ = 0.75;  // NTFD: hardcoded

/*std::list<std::pair<size_t,size_t> > filter;
filter.push_back(std::make_pair(0,100));
lid.setMassFilter(filter);*/

static const double LID_REQ_LABEL_AMOUNT = 0.05;
static const double LID_REQ_R2 = 0.95;
static const double LID_MIN_M0 = 0.0;
static const int LID_MIN_FRAG_NUM = 2;
// Removed in NTFD1.1 static const double LID_SENSITIVITY = 1/10000.0; // .0 !!
static const double LID_MAXIMAL_FRAG_DEV = 0.1;
static const double LID_CORRECTION_RATIO = 0.010934; // C tracer: correct for natural M+1 carbon abundance
static const int LID_MAX_MASS_ISOTOPOMER = 20;

static const double NW_GAP_PENALTY = 0.2; // Gap penalty for needleman wunsch scoring
static const bool NW_EXCLUDE_M0 = false;
static const double MID_DISTANCE_CUTOFF = 0.0; /** distance threshold for network edges */

static const bool NW_USE_LARGEST_COMMON_ION = false; /** Use largest *common* ion of group for network, instead individual largest */

}

#endif // CONFIG_H

