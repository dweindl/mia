//
// MIA - Mass Isotopolome Analyzer
// Copyright (C) 2013-15 Daniel Weindl <daniel@danielweindl.de>
//
// This file is part of MIA.
//
// MIA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// MIA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with MIA.  If not, see <http://www.gnu.org/licenses/>.
//

#include<sstream>
#include<list>

#include "settings.h"
#include "config.h"

namespace mia {

/**
 * @brief Default constructor setting default parameter values.
 */
Settings::Settings()
{
    cmp_id_use_ri = true;
    cmp_id_ri_tol = 100;
    cmp_id_score_cutoff = 0.75;
    cmp_id_library = "";

    cmp_matching_ri_tol = 10;
    cmp_matching_score_cutoff = 0.95;

    tracer_atom_mass = 13.0;
    gcms_pure_factor = 0.5;
    gcms_impure_factor = 0.5;

    lid_filter_by_conf_interval = true; // NTFD: hardcoded
    lid_min_signal_to_noise = 5; // NTFD: hardcoded
    lid_required_spec_freq = 0.75;  // NTFD: hardcoded

    std::list<std::pair<size_t,size_t> > filter;
    filter.push_back(std::make_pair(0,100)); // TODO used yet?

    lid_min_m0 = 0.45;
    lid_req_label_amount = 0.05;
    lid_req_r2 = 0.95;
    lid_min_frag_num = 2;
    lid_sensitivity = 1/10000.0; // .0 !!
    lid_maximal_frag_dev = 0.2;
    lid_correction_ratio = 0.010934; // C tracer
    lid_max_mass_isotopomer = 20;

    nw_gap_penalty = 0.2;
    nw_exclude_m0 = false;
    mid_distance_cutoff = 0.0;
    labels_max_hits = 2;
}

/**
 * @brief Copy constructor.
 * @param s Instance to copy from.
 */
Settings::Settings(const Settings &s)
{
    experiment = s.experiment;
    cmp_id_use_ri = s.cmp_id_use_ri;
    cmp_id_ri_tol = s.cmp_id_ri_tol;
    cmp_id_score_cutoff = s.cmp_id_score_cutoff;
    cmp_matching_ri_tol = s.cmp_matching_ri_tol;
    cmp_matching_score_cutoff = s.cmp_matching_score_cutoff;
    tracer_atom_mass = s.tracer_atom_mass;
    labels_max_hits = s.labels_max_hits;
    gcms_pure_factor = s.gcms_pure_factor;
    gcms_impure_factor = s.gcms_impure_factor;
    cmp_id_library = s.cmp_id_library;
    cmp_id_mass_filter = s.cmp_id_mass_filter;
    unlabFiles = s.unlabFiles;
    labFiles = s.labFiles;
    lid_filter_by_conf_interval = s.lid_filter_by_conf_interval;
    lid_min_signal_to_noise = s.lid_min_signal_to_noise;
    lid_required_spec_freq = s.lid_required_spec_freq;
    lid_req_label_amount = s.lid_req_label_amount;
    lid_req_r2 = s.lid_req_r2;
    lid_min_frag_num = s.lid_min_frag_num;
    lid_sensitivity = s.lid_sensitivity;
    lid_maximal_frag_dev = s.lid_maximal_frag_dev;
    lid_correction_ratio = s.lid_correction_ratio;
    lid_max_mass_isotopomer = s.lid_max_mass_isotopomer;
    lid_min_m0 = s.lid_min_m0;
    nw_gap_penalty = s.nw_gap_penalty;
    nw_exclude_m0 = s.nw_exclude_m0;
    mid_distance_cutoff = s.mid_distance_cutoff;
}

/**
 * @brief Check if certain parameter groups' settings are different from another instance.
 * @param other The other instance.
 * @param which Parameter group to be checked.
 * @return True if parameters differ.
 */
bool Settings::settingsChanged(const Settings &other, Settings::SettingsGroup which) const
{
    bool unchanged = true;

    if(which & CmpIdSettings) {
        unchanged = (unlabFiles == other.unlabFiles)
                && (labFiles == other.labFiles)
                && (cmp_id_use_ri == other.cmp_id_use_ri)
                && (cmp_id_ri_tol == other.cmp_id_ri_tol)
                && (cmp_id_score_cutoff == other.cmp_id_score_cutoff)
                && (labels_max_hits == other.labels_max_hits)
                && (cmp_id_library == other.cmp_id_library)
                && (cmp_matching_ri_tol == other.cmp_matching_ri_tol)
                && (cmp_matching_score_cutoff == other.cmp_matching_score_cutoff)
                && (cmp_id_mass_filter == other.cmp_id_mass_filter);
    }
    if(which & GCMSSettings) {
        unchanged = (gcms_pure_factor == other.gcms_pure_factor)
                && (gcms_impure_factor == other.gcms_impure_factor);
    }
    if(which & LidSettings) {
        unchanged = (lid_filter_by_conf_interval == other.lid_filter_by_conf_interval)
                && (lid_min_signal_to_noise == other.lid_min_signal_to_noise)
                && (lid_required_spec_freq == other.lid_required_spec_freq)
                && (lid_req_label_amount == other.lid_req_label_amount)
                && (lid_req_r2 == other.lid_req_r2)
                && (lid_min_frag_num == other.lid_min_frag_num)
                && (lid_sensitivity == other.lid_sensitivity)
                && (lid_maximal_frag_dev == other.lid_maximal_frag_dev)
                && (lid_correction_ratio == other.lid_correction_ratio)
                && (lid_min_m0 == other.lid_min_m0)
                && (tracer_atom_mass == other.tracer_atom_mass)
                && (lid_max_mass_isotopomer == other.lid_max_mass_isotopomer);
    }
    if(which & NWSettings) {
        unchanged = (nw_gap_penalty == other.nw_gap_penalty)
                && (nw_exclude_m0 == other.nw_exclude_m0)
                && (mid_distance_cutoff == other.mid_distance_cutoff);
    }
    if(which & AllSettings) {
        unchanged = (experiment == other.experiment);
    }
    return !unchanged;
}

/**
 * @brief Human-readable string of settings parameters.
 * @return Parameter string.
 */
std::string Settings::toString() const
{
    std::stringstream ss;
    ss << "experiment = " << experiment << std::endl;
    ss << "cmp_id_use_ri = " << cmp_id_use_ri << std::endl;
    ss << "cmp_id_ri_tol = " << cmp_id_ri_tol << std::endl;
    ss << "cmp_id_score_cutoff = " << cmp_id_score_cutoff << std::endl;
    ss << "tracer_atom_mass = " << tracer_atom_mass << std::endl;
    ss << "labels_max_hits = " << labels_max_hits << std::endl;
    ss << "gcms_pure_factor = " << gcms_pure_factor << std::endl;
    ss << "gcms_impure_factor = " << gcms_impure_factor << std::endl;
    ss << "cmp_id_library = " << cmp_id_library << std::endl;
    ss << "cmp_matching_ri_tol = " << cmp_matching_ri_tol << std::endl;
    ss << "cmp_matching_score_cutoff = " << cmp_matching_score_cutoff << std::endl;

    // TODO ss << "cmp_id_mass_filter = " << cmp_id_mass_filter << std::endl;
    for(std::vector<std::string>::const_iterator it = unlabFiles.begin(); it != unlabFiles.end(); ++it)
        ss << "unlabFiles = " << *it << std::endl;
    for(std::vector<std::string>::const_iterator it = labFiles.begin(); it != labFiles.end(); ++it)
        ss << "labFiles = " << *it << std::endl;
    ss << "lid_filter_by_conf_interval = " << lid_filter_by_conf_interval << std::endl;
    ss << "lid_min_signal_to_noise = " << lid_min_signal_to_noise << std::endl;
    ss << "lid_required_spec_freq = "<< lid_required_spec_freq << std::endl;
    ss << "lid_req_label_amount = " << lid_req_label_amount << std::endl;
    ss << "lid_req_r2 = " << lid_req_r2 << std::endl;
    ss << "lid_min_frag_num = " << lid_min_frag_num << std::endl;
    ss << "lid_sensitivity = " << lid_sensitivity << std::endl;
    ss << "lid_maximal_frag_dev = " << lid_maximal_frag_dev << std::endl;
    ss << "lid_correction_ratio = " << lid_correction_ratio << std::endl;
    ss << "lid_min_m0 = " << lid_min_m0 << std::endl;
    ss << "lid_max_mass_isotopomer = " << lid_max_mass_isotopomer << std::endl;
    ss << "nw_gap_penalty = " << nw_gap_penalty << std::endl;
    ss << "nw_exclude_m0 = " << nw_exclude_m0 << std::endl;
    ss << "mid_distance_cutoff = " << mid_distance_cutoff << std::endl;
    return ss.str();
}

std::string Settings::toXML() const
{
    std::stringstream ss;
    ss << "<experiment name=\"" << experiment << "\">" << std::endl;

    // settings
    ss << "<localsettings>"<<std::endl;

    ss << "<cmp_id_use_ri value=\"" << cmp_id_use_ri << "\"/>" << std::endl;
    ss << "<cmp_id_ri_tol value=\"" << cmp_id_ri_tol << "\"/>" << std::endl;
    ss << "<cmp_id_score_cutoff value=\"" << cmp_id_score_cutoff << "\"/>" << std::endl;
    ss << "<tracer_atom_mass value=\"" << tracer_atom_mass << "\"/>" << std::endl;
    ss << "<labels_max_hits value=\"" << labels_max_hits << "\"/>" << std::endl;
    ss << "<gcms_pure_factor value=\"" << gcms_pure_factor << "\"/>" << std::endl;
    ss << "<gcms_impure_factor value=\"" << gcms_impure_factor << "\"/>" << std::endl;
    ss << "<cmp_id_library value=\"" << cmp_id_library << "\"/>" << std::endl;
    ss << "<cmp_matching_ri_tol value=\"" << cmp_matching_ri_tol << "\"/>" << std::endl;
    ss << "<cmp_matching_score_cutoff value=\"" << cmp_matching_score_cutoff << "\"/>" << std::endl;
    // TODO cmp_id_mass_filter
    ss << "<lid_filter_by_conf_interval value=\"" << lid_filter_by_conf_interval << "\"/>" << std::endl;
    ss << "<lid_min_signal_to_noise value=\"" << lid_min_signal_to_noise << "\"/>" << std::endl;
    ss << "<lid_required_spec_freq value=\""<< lid_required_spec_freq << "\"/>" << std::endl;
    ss << "<lid_req_label_amount value=\"" << lid_req_label_amount << "\"/>" << std::endl;
    ss << "<lid_req_r2 value=\"" << lid_req_r2 << "\"/>" << std::endl;
    ss << "<lid_min_frag_num value=\"" << lid_min_frag_num << "\"/>" << std::endl;
    ss << "<lid_sensitivity value=\"" << lid_sensitivity << "\"/>" << std::endl;
    ss << "<lid_maximal_frag_dev value=\"" << lid_maximal_frag_dev<< "\"/>"  << std::endl;
    ss << "<lid_correction_ratio value=\"" << lid_correction_ratio << "\"/>" << std::endl;
    ss << "<lid_min_m0 value=\"" << lid_min_m0 << "\"/>" << std::endl;
    ss << "<lid_max_mass_isotopomer value=\"" << lid_max_mass_isotopomer << "\"/>" << std::endl;
    ss << "<nw_gap_penalty value=\"" << nw_gap_penalty << "\"/>" << std::endl;
    ss << "<nw_exclude_m0 value=\"" << nw_exclude_m0 << "\"/>" << std::endl;
    ss << "<mid_distance_cutoff value=\"" << mid_distance_cutoff << "\"/>" << std::endl;

    ss << "</localsettings>"<<std::endl;

    // Data files
    ss<<"<labeledfiles>"<<std::endl;
    for(std::vector<std::string>::const_iterator it = labFiles.begin(); it != labFiles.end(); ++it)
        ss << "<file name=\"" << *it << "\"/>" << std::endl;
    ss<<"</labeledfiles>"<<std::endl;

    ss<<"<unlabeledfiles>"<<std::endl;
    for(std::vector<std::string>::const_iterator it = unlabFiles.begin(); it != unlabFiles.end(); ++it)
        ss << "<file name=\"" << *it << "\"/>" << std::endl;
    ss<<"</unlabeledfiles>"<<std::endl;

    ss << "</experiment>" << std::endl;
    return ss.str();
}

}
