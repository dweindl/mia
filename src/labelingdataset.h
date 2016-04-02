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

#ifndef LABELINGDATASET_H
#define LABELINGDATASET_H

#include <QDataStream>
#include <vector>

#include "labelidentificator.h"

#include "labelingdataset.h"
#include "miaexception.h"
#include "settings.h"
#include "config.h"

#include "../rapidxml/rapidxml.hpp"

namespace mia
{

// class MyLabelIdentificator ?
// TODO: select / deselect detected ions: unselect on filtering, do not remove
//
// Wrap in FilteredLabelingDataset for use with multiple views?
// LabelingDataset(labid::LabelIdentificator)
// LabelingDataset(textfile?)

class LabelingDataset;

/**
 * @brief The LabelingDataset detects labeled compounds from labeled and unlabeled reference spectra based on the labid::LabelIdentificator class.
 * It calculates MID similarities and holds the distance matrices, allows for redetection of given ions.
 */

class LabelingDataset
{

public:

    LabelingDataset();
    LabelingDataset(Settings);

    // TODO getNearestNeighborGraph

    void setExcludeLib(gcms::LibrarySearch<int,float> *excludeLib);
    bool isExcludedMetabolite(const gcms::Compound<int, float> cmp);

    static std::vector<LabelingDataset *> fromXMLFile(std::string file);

    void removeScoresBelow();

    static std::vector<std::vector<double> > removeScoresBelow(std::vector<std::vector<double> > dists, double cutoff);

    static std::vector<std::vector<double> > removeScoresAbove(std::vector<std::vector<double> > dists, double cutoff);

    static std::vector<double> basePeakNormalization(const std::vector<double> &v);

    static void distMatsToDot(std::vector<std::vector<std::vector<double> > > distMats, std::string fname, std::vector<std::string> labels);

    void distsOverviewHTML(std::string fname);

    void doScoring();

    void findLabeledCompounds(labid::LabelIdentificatorProgressListener *listener = 0);

    void saveMIDsCsv(std::vector<std::vector<double> > mids, std::vector<std::string> names, std::string fname);

    void fetchLabeledCompoundsFromDB(int id, std::vector< labid::LabeledCompound*> &lcs);

    labid::LabelIdentificator * getLabelIdentificator() {return lid;}

    //void identifyLabeledCompounds();
    static std::vector<gcms::Compound<int, float> *> identifyCompounds(std::set<gcms::Compound<int, float> *>, Settings settings);

    const Settings &getSettings() const;
    void setSettings(Settings s);

    friend QDataStream &operator << (QDataStream &out, const LabelingDataset*);
    friend QDataStream &operator >> (QDataStream &in, LabelingDataset*&) throw(DeserializationException);

    std::vector<labid::LabeledCompound*> cmpLab;    /**< Detected labeled compounds. */
    std::vector<labid::LISpectrum*> cmpUnlab;       /**< Detected unlabeled compounds. */
    std::vector<std::vector<double> > dists;        /**< Distance matrix. */
    std::vector<std::vector<double> > distsCut;     /**< Adjacency matrix. (Distance matrix after cutoff applied. */

    std::vector<std::vector<double> > mids;         /**< Selected fragment MIDs for network. */
    std::vector<std::vector<double> > midsAll;      /**< All MIDs for debug. */
    std::vector<std::string> nodeLabs;              /**< Labels for network nodes. Compound names if library matching was done. */
    std::vector<std::string> nodeLabsAll;           /**< All labels for debug. */

protected:
    Settings settings;              /**< Settings for the compound detection and distance calculation, ... */
    labid::LabelIdentificator *lid; /**< The LabelIdentificator object used for labeled compound detection. */

private:
    // TODO: move xml functions to mia::Settings
    static void parseXMLSettings(rapidxml::xml_node<> *nodeSettings, Settings &s);
    static LabelingDataset *parseXMLExperiment(rapidxml::xml_node<> *nodeExperiment, Settings settings, std::string curDir = "");
    static std::vector<std::string> parseXMLFileSection(rapidxml::xml_node<> *nodeFiles, bool reportNonExistance = false, std::string curDir = "");
    static int parseXMLGetInt(rapidxml::xml_node<> *node, std::string attr = "value", bool mandatory = true);
    static double parseXMLGetDouble(rapidxml::xml_node<> *node, std::string attr = "value", bool mandatory = true);
    static bool parseXMLGetBool(rapidxml::xml_node<> *node, std::string attr = "value", bool mandatory = true);
    static std::string parseXMLGetString(rapidxml::xml_node<> *node, std::string attr = "value", bool mandatory = true);

    gcms::LibrarySearch<int,float> *excludeLib;
};

}
#endif // LABELINGDATASET_H
