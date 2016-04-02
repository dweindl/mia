/* * MIA - Mass Isotopolome Analyzer
 * Copyright (C) 2013-15 Daniel Weindl <daniel@danielweindl.de>
 *
 * This file is part of MIA.
 *
 * MIA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * MIA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with MIA.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SETTINGS_H
#define SETTINGS_H

#include <set>
#include <iostream>
#include <vector>

namespace mia {

/**
 * @class Settings
 * @brief The Settings class holds all parameters for labeled compound detection,
 * identification, distance calculation and network generation.
 * @author Daniel Weindl
 */
class Settings
{
public:
    /**
     * @brief The SettingsGroup enum. Groups parameters into settings for gcms-lib, labid-lib,
     * compound identification and network generation.
     */
    enum SettingsGroup {
        GCMSSettings = 1,
        LidSettings = 2,
        CmpIdSettings = 4,
        NWSettings = 8,
        AllSettings = 0xFF
    };

    Settings();
    Settings(const Settings&);

    bool settingsChanged(const Settings& other, SettingsGroup which = AllSettings) const;
    std::string toString() const;
    std::string toXML() const;

    std::string experiment;         /**< Name of the dataset. */

    // Compound identification settings
    bool cmp_id_use_ri;         /**< Use retention index. */
    double cmp_id_ri_tol;       /**< Tolerance for retention index difference. */
    double cmp_id_score_cutoff; /**< Cutoff for spectrum matching score for identification. */
    double cmp_matching_ri_tol; /**< RI tolerance for peak matching different chromatograms */
    double cmp_matching_score_cutoff; /**< Spec score for peak matching different chromatograms */

    int labels_max_hits;        /**< Maximum number of labels to include in the name. -1: show all. */
    double gcms_pure_factor;    /**< See labid. */
    double gcms_impure_factor;  /**< See labid. */
    std::set<std::pair<float,float> > cmp_id_mass_filter; /** m/z range to exclude from library matching. */
    std::string cmp_id_library; /**< MD library file. */
    double tracer_atom_mass;    /**< Mass of the tracer isotope. To be used for filtering. */

    std::vector<std::string> unlabFiles; /**< Unlabeled files for NTFD. */
    std::vector<std::string> labFiles;  /**< Labeled files for NTFD. */

    // LabelIdentificator Settings
    bool lid_filter_by_conf_interval;   /**< NTFD: hardcoded */
    int lid_min_signal_to_noise;        /**< NTFD: hardcoded */
    double lid_required_spec_freq;      /**< NTFD: hardcoded */
    double lid_req_label_amount;        /**< Minimum sum of mass isotopomer abundances to keep ion in results list. */
    double lid_req_r2;                  /**< Filter by R^2 < lid_req_r2. */
    int lid_min_frag_num;               /**< Minimum number of labeled fragments a spectrum needs to hold. */
    double lid_sensitivity;             /**< Specificity for label detection. */
    double lid_maximal_frag_dev;        /**< Filter fragments where sum of absolute mass isotopomer abundances is > 1 + lid_maximal_frag_dev. */
    double lid_correction_ratio;        /**< M+1 correction. Use 0.010934 for C tracer. M+1/M+0. See Jennings et al. */
    int lid_max_mass_isotopomer;        /**< Filter ions with M+n | n > lid_max_mass_isotopomer. Defaults to 20. */
    double lid_min_m0;                  /**< Minimum M0 abundance (use little less than tracer percentage). Defaults to 0.45 (for 50:50 tracer:non-tracer). */
    /* TODO lid_mass_filter * m/z range to exclude from labeled ion detection. */

    // Network construction settings
    double nw_gap_penalty;              /**< Gap penalty for needleman wunsch scoring. */
    bool nw_exclude_m0;                 /**< Exclude M+0 abundance for distance calculation & scoring? */
    double mid_distance_cutoff;         /**< Distance threshold for network edges (lower/upper?). */
};

}
#endif // SETTINGS_H
