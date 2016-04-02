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

#include "misc.h"

#include <sstream>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <locale>

#include <QString>
#include <QVector>
#include <QStringList>
#include <QDataStream>

#include "librarysearch.h"
#include "libraryhit.h"
#include "gcmssettings.h"
#include "compound.h"

#include "nodecompound.h"
#include "labelingdataset.h"
#include "utilities.h"



namespace mia
{

/**
 * @brief Default constructor using default settings.
 */
LabelingDataset::LabelingDataset()
{
    settings = Settings();
    lid = 0;
    excludeLib = 0;
}

/**
 * @brief Constructor
 * @param _settings Settings to be used for compound detection, ...
 */
LabelingDataset::LabelingDataset(Settings _settings)
{
    settings = _settings;
    lid = 0;
    excludeLib = 0;
}

void LabelingDataset::setExcludeLib(gcms::LibrarySearch<int, float> *_excludeLib)
{
    excludeLib = _excludeLib;
}

bool LabelingDataset::isExcludedMetabolite(const gcms::Compound<int, float> cmp)
{
    if(!excludeLib)
        return false;

    gcms::GCMSSettings::LS_USE_RI        = true;
    gcms::GCMSSettings::LS_RI_DIFF       = 100; // TODO: to config

    std::vector<gcms::LibraryHit<int,float> > hits = excludeLib->getLibraryHits(cmp);

    double mylibScoreCutoff = 0.92; // TODO: to config

    if(hits.size() && hits.at(0).getOverallScore() >=  mylibScoreCutoff) {
        return true;
    }

    return false;
}

std::vector<LabelingDataset *> LabelingDataset::fromXMLFile(std::string file)
{
    std::cout<<"Opening "<<file<<std::endl;

    // read file into buffer
    std::ifstream ifs;
    ifs.open(file.c_str(), std::ios_base::binary|std::ios_base::ate);
    if(!ifs.is_open()) {
        std::cerr<<"Error opening file "<<file<<std::endl;
    }

    std::string xmlDir = file.substr(0, file.find_last_of("\\/") + 1);

    int filesize = ifs.tellg();
    char *xmldata = new char[filesize + 1];
    xmldata[filesize] = 0;
    ifs.seekg(0, std::ios_base::beg);
    ifs.read(xmldata, filesize);
    ifs.close();

    // parse XML
    rapidxml::xml_document<> doc;
    doc.parse<0>(xmldata); // 0 means default parse flags

    rapidxml::xml_node<> *nodeExperiments = doc.first_node("experiments");
    if(!nodeExperiments) throw MIAException("experiments");

    Settings globalSettings;
    globalSettings.cmp_id_mass_filter.insert(std::make_pair (0.0, 100.0));

    rapidxml::xml_node<> *nodeGlobalSettings = nodeExperiments->first_node("globalsettings");
    if(nodeGlobalSettings) {
        parseXMLSettings(nodeGlobalSettings, globalSettings);
    }

    std::vector<LabelingDataset *> LabelingDatasets;

    for (rapidxml::xml_node<> *nodeExperiment = nodeExperiments->first_node("experiment"); nodeExperiment; nodeExperiment = nodeExperiment->next_sibling("experiment")) {
        LabelingDataset *dataset = parseXMLExperiment(nodeExperiment, globalSettings, xmlDir);
        LabelingDatasets.push_back(dataset);
    }

    doc.clear();
    delete[] xmldata;

    return LabelingDatasets;
}

void LabelingDataset::parseXMLSettings(rapidxml::xml_node<> *nodeSettings, Settings &s)
{
    for (rapidxml::xml_node<> *nodeParam = nodeSettings->first_node(); nodeParam; nodeParam = nodeParam->next_sibling()) {

        std::string param = nodeParam->name();

        if(param == "cmp_id_library") {
            s.cmp_id_library = parseXMLGetString(nodeParam);
            if(!Utilities::fileExists(s.cmp_id_library)) {
                std::cerr<<"Compound library "<<s.cmp_id_library<<" is not readable."<<std::endl;
            }
        } else if (param == "nw_gap_penalty") {
            s.nw_gap_penalty = parseXMLGetDouble(nodeParam);
        } else if (param == "cmp_id_score_cutoff") {
            s.cmp_id_score_cutoff = parseXMLGetDouble(nodeParam);
        } else if (param == "lid_maximal_frag_dev") {
            s.lid_maximal_frag_dev = parseXMLGetDouble(nodeParam);
        } else if (param == "mid_distance_cutoff") {
            s.mid_distance_cutoff = parseXMLGetDouble(nodeParam);
        } else if (param == "cmp_id_ri_tol") {
            s.cmp_id_ri_tol = parseXMLGetDouble(nodeParam);
        } else if (param == "lid_sensitivity") {
            s.lid_sensitivity = parseXMLGetDouble(nodeParam);
        } else if (param == "lid_min_frag_num") {
            s.lid_min_frag_num = parseXMLGetInt(nodeParam);
        } else if (param == "nw_exclude_m0") {
            s.nw_exclude_m0 = parseXMLGetBool(nodeParam);
        } else if (param == "labels_max_hits") {
            s.labels_max_hits = parseXMLGetInt(nodeParam);
        } else if (param == "lid_min_m0") {
            s.lid_min_m0 = parseXMLGetDouble(nodeParam);
        } else if (param == "tracer_atom_mass") {
            s.tracer_atom_mass = parseXMLGetDouble(nodeParam);
        } else if (param == "lid_req_r2") {
            s.lid_req_r2 = parseXMLGetDouble(nodeParam);
        } else if (param == "lid_req_label_amount") {
            s.lid_req_label_amount = parseXMLGetDouble(nodeParam);
        } else if (param == "cmp_id_use_ri") {
            s.cmp_id_use_ri = parseXMLGetBool(nodeParam);
        } else if (param == "lid_correction_ratio") {
            s.lid_correction_ratio = parseXMLGetDouble(nodeParam);
        } else if (param == "gcms_pure_factor") {
            s.gcms_pure_factor = parseXMLGetDouble(nodeParam);
        } else if (param == "gcms_impure_factor") {
            s.gcms_impure_factor = parseXMLGetDouble(nodeParam);
        } else if (param == "cmp_matching_ri_tol") {
            s.cmp_matching_ri_tol = parseXMLGetDouble(nodeParam);
        } else if (param == "cmp_matching_score_cutoff") {
            s.cmp_matching_score_cutoff = parseXMLGetDouble(nodeParam);
        } else if (param == "lid_filter_by_conf_interval") {
            s.lid_filter_by_conf_interval = parseXMLGetBool(nodeParam);
        } else if (param == "lid_min_signal_to_noise") {
            s.lid_min_signal_to_noise = parseXMLGetDouble(nodeParam);
        } else if (param == "lid_required_spec_freq") {
            s.lid_required_spec_freq = parseXMLGetDouble(nodeParam);
        } else if (param == "lid_max_mass_isotopomer") {
            s.lid_max_mass_isotopomer = parseXMLGetInt(nodeParam);
        } else {
            std::cerr<<"Unknown settings parameter "<<param<<std::endl;
        }
    }
}

LabelingDataset *LabelingDataset::parseXMLExperiment(rapidxml::xml_node<> *nodeExperiment, Settings settings, std::string curDir)
{
    rapidxml::xml_attribute<> *attr = nodeExperiment->first_attribute("name");
    if(!attr || !(settings.experiment = attr->value()).length()) {
        std::stringstream ss;
        ss<<std::time(0);
        settings.experiment =  ss.str();
        std::cerr<<"Experiment has no name, using "<<settings.experiment<<std::endl;
    }

    rapidxml::xml_node<> *nodeLocalSettings = nodeExperiment->first_node("localsettings");
    if(nodeLocalSettings) {
        parseXMLSettings(nodeLocalSettings, settings);
    }

    rapidxml::xml_node<> *nodeLabeledFiles = nodeExperiment->first_node("labeledfiles");
    if(nodeLabeledFiles) {
        settings.labFiles = parseXMLFileSection(nodeLabeledFiles, true, curDir);
    }

    rapidxml::xml_node<> *nodeUnlabeledFiles = nodeExperiment->first_node("unlabeledfiles");
    if(nodeUnlabeledFiles) {
        settings.unlabFiles = parseXMLFileSection(nodeUnlabeledFiles, true, curDir);
    }

    return new LabelingDataset(settings);
}

std::vector<std::string> LabelingDataset::parseXMLFileSection(rapidxml::xml_node<> *nodeFiles, bool reportNonExistance, std::string curDir)
{
    std::vector<std::string> files;

    for (rapidxml::xml_node<> *nodeFile = nodeFiles->first_node("file"); nodeFile; nodeFile = nodeFile->next_sibling("file")) {
        rapidxml::xml_attribute<> * attr = nodeFile->first_attribute("name");

        if(!attr)
            throw MIAException("experiments");

        std::string fileName = attr->value();

        if(fileName[0] != '/' || fileName[1] != ':') { // is absolute or relative path?
            std::cout<<"Assuming relative path for "<<fileName<<" in "<<curDir<<std::endl;
            fileName = curDir + fileName;
        }

        if(reportNonExistance && !Utilities::fileExists(fileName))
            std::cerr<<"File "<<fileName<<" is not readable."<<std::endl;

        files.push_back(fileName);
    }

    return files;
}

int LabelingDataset::parseXMLGetInt(rapidxml::xml_node<> *node, std::string attr, bool mandatory)
{
    rapidxml::xml_attribute<> *a = node->first_attribute(attr.c_str());

    if(!a) {
        std::stringstream ss;
        ss<<"Attribute "<<attr<< " not found for node "<<node->name()<<std::endl;
        std::cerr<<ss.str();
        //throw LabelingDatasetException(ss.str());
        return 0;
    }

    return atoi(a->value());
}

double LabelingDataset::parseXMLGetDouble(rapidxml::xml_node<> *node, std::string attr, bool mandatory)
{
    rapidxml::xml_attribute<> *a = node->first_attribute(attr.c_str());

    if(!a) {
        std::stringstream ss;
        ss<<"Attribute "<<attr<< " not found for node "<<node->name()<<std::endl;
        std::cerr<<ss.str();
        //throw LabelingDatasetException(ss.str());
        return 0;
    }

    // replace "." decimal separator by the one used in current locale
    std::string tmp = a->value();
    size_t ret = tmp.find_first_of('.');
    if(ret != std::string::npos) {
        lconv *lconv = localeconv();
        tmp[ret] = *(lconv->decimal_point);
    }

    return (float)atof(tmp.c_str());
}

bool LabelingDataset::parseXMLGetBool(rapidxml::xml_node<> *node, std::string attr, bool mandatory)
{
    rapidxml::xml_attribute<> *a = node->first_attribute(attr.c_str());

    if(!a) {
        std::stringstream ss;
        ss<<"Attribute "<<attr<< " not found for node "<<node->name()<<std::endl;
        std::cerr<<ss.str();
        //throw LabelingDatasetException(ss.str());
        return false;
    }

    std::string val = a->value();
    if(val == "true" || val == "1") {
        return true;
    } else if(val == "false" || val == "0") {
        return false;
    }
    std::cerr<<"Attribute "<<attr<< " for node "<<node->name()<<": "<<val<<": invalid value for type bool."<<std::endl;

    return false;
}

std::string LabelingDataset::parseXMLGetString(rapidxml::xml_node<> *node, std::string attr, bool mandatory)
{
    rapidxml::xml_attribute<> *a = node->first_attribute(attr.c_str());

    if(!a) {
        std::stringstream ss;
        ss<<"Attribute "<<attr<< " not found for node "<<node->name()<<std::endl;
        std::cerr<<ss.str();
        //throw LabelingDatasetException(ss.str());
        return "";
    }

    return a->value();
}


void LabelingDataset::removeScoresBelow()
{
    distsCut = removeScoresBelow(dists, 4);
}


std::vector<std::vector<double> > LabelingDataset::removeScoresBelow(std::vector<std::vector<double> > dists, double cutoff) {
    std::vector<std::vector<double> > distsCut;
    // apply cutoff i.e. filter out everything below
    distsCut.resize(dists.size());
    for(int i = 0; i < distsCut.size(); ++i) {
        distsCut[i].resize(dists[i].size());
        for(int j = 0; j < distsCut[i].size(); ++j) {
            if(i == j) // diagonal
                distsCut[i][j] = 1;
            else if (i < j) // do only upper half
                distsCut[i][j] = (dists[i][j] >= 5)?dists[i][j]:0;
        }
    }

    return distsCut;
}


std::vector<std::vector<double> > LabelingDataset::removeScoresAbove(std::vector<std::vector<double> > dists, double cutoff) {
    std::vector<std::vector<double> > distsCut;
    // apply cutoff i.e. filter out everything below
    distsCut.resize(dists.size());
    for(int i = 0; i < distsCut.size(); ++i) {
        distsCut[i].resize(dists[i].size());
        for(int j = 0; j < distsCut[i].size(); ++j) {
            if(i == j) // diagonal
                distsCut[i][j] = 1;
            else if (i < j) // do only upper half
                distsCut[i][j] = (dists[i][j] < cutoff)?dists[i][j]:0;
        }
    }

    return distsCut;
}

/**
 * @brief Normalize vector to highest intensity = 1.
 * @param v Vector to normalize.
 * @return The normalized vector.
 */
std::vector<double> LabelingDataset::basePeakNormalization(const std::vector<double> &v)
{
    double base = 0;
    for(std::vector<double>::const_iterator it = v.begin(); it != v.end(); ++it) {
        base = std::max(base, *it);
    }

    std::vector<double> vv = v;
    for(std::vector<double>::iterator it = vv.begin(); it != vv.end(); ++it) {
        *it /= base;
    }
}

/**
 * @brief Write simple Dot-format file for graphviz visualization.
 *
 * @param mat    Network matrix. Values > 0 represent nodes
 * @param fname  Output filename
 * @param labs   Node labels
 */
void LabelingDataset::distMatsToDot(std::vector<std::vector<std::vector<double> > > distMats, std::string fname, std::vector<std::string> labels)
{
    std::ofstream ofs;
    std::cerr<<"writing " << fname<<std::endl;
    ofs.open(fname.c_str());
    std::vector<std::string> colors;
    colors.push_back("blue");
    colors.push_back("red");
    colors.push_back("darkgreen");
    colors.push_back("orange");
    colors.push_back("pink");

    ofs << "graph G {" << std::endl;

    // write edges for each matrix:
    for(int t = 0; t < distMats.size(); ++t) {
        // each tracer
        std::vector<std::vector<double> > mat = distMats[t];
        for(int i = 0; i < mat[0].size(); ++i) { // row
            for(int j = 0; j < mat.size(); ++j) { // col
                if(mat[i][j] > 0 && j > i) {
                    // i--j [weight=..];
                    double weight = 1 / mat[i][j];
                    //double penwidth = log(1 / mat[i][j]) * 4 / distRange; // maxwidth 4
                    //ofs << i << " -- " << j << " [weight=" << weight << ",penwidth=" << "1" <<"];" << std::endl;
                    ofs << i << " -- " << j << " [penwidth=" << "2" <<",label="<<mat[i][j]<<", color="<<colors[t]<<"];" << std::endl;
                }
            }
        }
    }

    // write nodes only once
    // extended labels
    for(int i = 0; i < labels.size(); ++i) {
        ofs << i << " [";
        ofs << "label=<<TABLE border=\"0\" cellborder=\"0\"><TR><TD COLSPAN=\""<<distMats.size()<<"\">"<<labels[i]<<"</TD></TR><TR>";
        for(int j = 0; j < distMats.size(); ++j) {
            std::stringstream imgname;
            imgname << "midplotsoverlay/"<<i<<"_"<<j<<".png";
            std::ifstream ifs(imgname.str().c_str());
            if(ifs) ofs << "<TD><IMG SRC=\""<<imgname.str()<<"\"/></TD>";
            else ofs << "<TD></TD>";
        }
        ofs << "</TR></TABLE>>,style=filled,fillcolor=\"#ACD9FF\",";
        ofs << "];" << std::endl;
    }

    ofs << "}" << std::endl; // close graph
    ofs.close();
}

void LabelingDataset::distsOverviewHTML(std::string fname) {
    std::ofstream ofs;
    std::cerr<<"writing " << fname<<std::endl;
    ofs.open(fname.c_str());

    ofs << "<html>\n<head>\n<title>dists</title>\n</head>\n<body>\n";
    ofs << "<table border=1>";
    for(int i = 0; i < dists.size(); ++i) {
        for(int j = 0; j < dists[i].size(); ++j) {
            if(i >= j) continue;
            ofs <<"<tr><td><img src=\"midplotssingle/"<<i<<".png\"></td><td>"<<dists[i][j]<<"</td><td><img src=\"midplotssingle/"<<j<<".png\"></td>"<<"<td>";
            if(settings.nw_exclude_m0) { // skip M0
                ofs << nwHTML<double>(basePeakNormalization(std::vector<double>(&(mids[i][1]), &(mids[i][mids[i].size() - 1]))),
                        basePeakNormalization(std::vector<double>(&(mids[j][1]), &(mids[j][mids[j].size() - 1]))),
                        settings.nw_gap_penalty, ofs
                        );
            } else {
                ofs << nwHTML<double>(mids[i], mids[j], settings.nw_gap_penalty, ofs);
            }

            ofs <<"</td></tr>\n";
        }
    }
    ofs << "</table>";
    ofs << "</body>\n";
}

// TODO: use middistancecalculator
void LabelingDataset::doScoring()
{
    // keep some stats on distances:
    double dMin = 0;
    double dMax = 0;
    double dSum = 0;

    // calc needleman-wunsch scores
    dists.resize(mids.size());
    for(int i = 0; i < dists.size(); ++i) {
        dists[i].resize(mids.size());
        for(int j = 0; j < dists[i].size(); ++j) {
            if(i == j) { // diagonal
                dists[i][j] = 1;
            } else if (i < j) { // do only upper half
                if(settings.nw_exclude_m0) { // skip M0
                    dists[i][j] = nw<double>(basePeakNormalization(std::vector<double>(&(mids[i][1]), &(mids[i][mids[i].size() - 1]))),
                            basePeakNormalization(std::vector<double>(&(mids[j][1]), &(mids[j][mids[j].size() - 1]))),
                            settings.nw_gap_penalty
                            );
                } else {
                    dists[i][j] = nw<double>(mids[i], mids[j], settings.nw_gap_penalty);
                }
                // stats
                dMin = dists[i][j] < dMin ? dists[i][j] : dMin;
                dMax = dists[i][j] > dMax ? dists[i][j] : dMax;
                dSum += dists[i][j];
            }
        }
    }
    // problem: nw scores not comparable...
    // normalize by score / meanlength?
    // use euclidean distance after alignment?
    std::cerr<< "Distances ("<<dists.size()<<")\n\tRange: "<<dMin<<" - "<<dMax<<"\n\tMean: "<<dSum/(dists.size()*(dists.size()-1))<<"\n";
}

/**
 * @brief Do the labeled-compound-detection. Sets mia::cmpLab and mia::cmpUnlab and mia::midCompounds.
 */
void LabelingDataset::findLabeledCompounds(labid::LabelIdentificatorProgressListener *listener)
{
    if(lid) {
        delete lid;
        lid = 0;
    }
    std::cout<<"Started label detection thread.\n";
    lid  = new labid::LabelIdentificator(settings.unlabFiles, settings.labFiles);
    lid->setMaximumLabel(1.0 - settings.lid_min_m0);
    lid->setEnsureM0PresenceInLabeledSpec(true);
    lid->setEnsureM1LessThanM0(true);
    lid->setFragmentDetectionPlateauTolerance(1.1);
    lid->setMaxMMinusOne(0.25);
    lid->setMinimumM1(0.01);
    if(listener)
        lid->setProgressListener(listener);

    // DW: No TMS-O-TMS, no < 80
    std::list<std::pair<size_t,size_t> > filter;
    filter.push_back(std::pair<size_t,size_t>(0, 80));
    filter.push_back(std::pair<size_t,size_t>(147, 147));
    lid->setMassFilter(filter);

    /* LabelIdentificator Settings */
    lid->setFilterByConfidenceInterval(settings.lid_filter_by_conf_interval);
    lid->setMinimalSN(settings.lid_min_signal_to_noise);
    lid->setRequiredSpectrumFrequency(settings.lid_required_spec_freq);
    // TODO    lid->setMassFilter(settings.lid_mass_filter);*/
    lid->setRequiredLabelAmount(settings.lid_req_label_amount);
    lid->setRequiredR2(settings.lid_req_label_amount);
    lid->setMinimalNumberOfFragments(settings.lid_min_frag_num);
    lid->setMaximalFragmentDeviation(settings.lid_maximal_frag_dev);
    lid->setCorrectionRatio(settings.lid_correction_ratio);
    /* LabelIdentificator Settings */

    lid->startAnalysis();

    cmpLab = lid->getLabeledCompounds();
    cmpUnlab = lid->getUnlabeledCompounds();

    // Set nicer compound names
    for(std::vector<labid::LabeledCompound*>::iterator it = cmpLab.begin(); it != cmpLab.end(); ++it) {
        labid::LISpectrum* cmp = *it;
        std::stringstream ss;
        ss.precision(6);
        if(cmp->getRetentionIndex() > 0)
            ss<<"Unidentified RI"<<cmp->getRetentionIndex();
        else
            ss<<"Unidentified RT"<<(cmp->getRetentionTime() / 1000.0 / 60.0);
        cmp->setName(ss.str());
    }
    for(std::vector<labid::LISpectrum*>::iterator it = cmpUnlab.begin(); it != cmpUnlab.end(); ++it) {
        labid::LISpectrum* cmp = *it;
        std::stringstream ss;
        ss.precision(6);
        if(cmp->getRetentionIndex() > 0)
            ss<<"Unidentified RI"<<cmp->getRetentionIndex();
        else
            ss<<"Unidentified RT"<<(cmp->getRetentionTime() / 1000.0 / 60.0);
        cmp->setName(ss.str());
    }

    std::cout<<"Finished label detection thread.\n";
}

/**
 * @brief Match compounds against a library.
 * @param comps Compounds to match
 * @param settings Settings for spectrum matching including library path.
 * @return Compounds with name set according to identification.
 */
std::vector<gcms::Compound<int, float> *> LabelingDataset::identifyCompounds(std::set<gcms::Compound<int, float> *> comps, mia::Settings settings) {
    // gcms settings for identification (mapCompoundsToLibrary overrides LS_CUTOFF
    gcms::GCMSSettings::LS_USE_RI        = settings.cmp_id_use_ri;
    gcms::GCMSSettings::LS_RI_DIFF       = settings.cmp_id_ri_tol;
    gcms::GCMSSettings::LS_PURE_FACTOR   = settings.gcms_pure_factor;
    gcms::GCMSSettings::LS_IMPURE_FACTOR = settings.gcms_impure_factor;
    gcms::GCMSSettings::LS_MASS_FILTER   = settings.cmp_id_mass_filter;

    // Library matching
    std::cout<<"Loading library: "<<settings.cmp_id_library<<std::endl;

    gcms::LibrarySearch<int,float> *lib(gcms::LibrarySearch<int,float>::fromDisk(settings.cmp_id_library.c_str()));
    std::vector<gcms::Compound<int, float> *> compsV;

    for(std::set<gcms::Compound<int, float>* >::iterator it = comps.begin(); it != comps.end(); ++it) {

        gcms::Compound<int,float> *cmp = *it;
        std::vector<gcms::LibraryHit<int,float> > hits = lib->getLibraryHits(*cmp);

        std::stringstream ss;

        int hitCount = settings.labels_max_hits; // how many possible identifications are already appended?

        for (std::vector<gcms::LibraryHit<int,float> >::const_iterator hitIt = hits.begin(); hitIt != hits.end() && hitCount--; ++hitIt) {

            if(hitIt->getOverallScore() < settings.cmp_id_score_cutoff) break;

            // copy features from best hit (need e.g. KEGG ids later)
            if(hitIt == hits.begin()) {
                std::map<std::string, std::string> libFeatures = hitIt->getLibraryCompound()->getFeatures();
                libFeatures.insert(cmp->getFeatures().begin(), cmp->getFeatures().end()); // preserve previously set features
                cmp->setFeatures(libFeatures);
            }

            ss<<hitIt->getLibraryCompound()->getName()<<"("<<hitIt->getOverallScore()<<")";
#define LabelingDataset_H_DBG_CMP_ID
#ifdef LabelingDataset_H_DBG_CMP_ID
            std::cout<<"ID"<<cmp->getFeature(COMPOUND_GROUPING_FEATURE)<<" identified as "
                    <<hitIt->getLibraryCompound()->getName()<<"("<<hitIt->getOverallScore()<<")" <<
                      hitIt->getRISimilarityScore() << " " << hitIt->getSpectrumSimilarityScore()<<std::endl;
#endif

        }

        if(!ss.str().length()){
            ss<<"NO MATCH "<<cmp->getRetentionIndex();
        }
        cmp->setName(ss.str());

#ifdef LabelingDataset_H_DBG_CMP_ID
        std::cout<<cmp->getRetentionIndex()<<" "<<cmp->getFeature(COMPOUND_GROUPING_FEATURE)<<" "<<cmp->getName()<<std::endl;
#endif

        compsV.push_back(cmp);
    }
    delete lib;

    return compsV;
}


void LabelingDataset::saveMIDsCsv(std::vector<std::vector<double> > mids, std::vector<std::string> names, std::string fname) {
    if(mids.size() != names.size()) {
        std::cerr<<"Error: Size mismatch\n";
        return;
    }

    std::ofstream ofs;
    ofs.open(fname.c_str());
    for(int i = 0; i < mids.size(); ++i) {// compounds
        std::vector<double> mid = mids[i];
        for(int j = 0; j < mid.size(); ++j) { // isotopomers
            ofs << "\""<< names[i] << "\",";
            ofs << "\""<< j << "\",";
            ofs << "\"" << mid[j] << "\",";
            ofs << std::endl;
        }
    }
    ofs.close();
}

/**
 * @brief Returns the current Settings.
 * @return Settings
 */
const Settings &LabelingDataset::getSettings() const
{
    return settings;
}

/**
 * @brief Set new settings.
 * @param s Settings
 */
void LabelingDataset::setSettings(Settings s)
{
    settings = s;
}

}
