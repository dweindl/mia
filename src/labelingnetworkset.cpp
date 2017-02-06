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

#include <sstream>
#include "labelingnetworkset.h"
#include "misc.h"

namespace mia {

MIDDistanceCalculator *LabelingNetworkSet::distCalc = 0;

LabelingNetworkSet::LabelingNetworkSet()
{
    excludeM0 = 0;
}

LabelingNetworkSet::~LabelingNetworkSet()
{
    for(int i = 0; i < nodes.size(); ++i) {
        delete nodes[i];
    }

    for(int i = 0; i < datasets.size(); ++i) {
        delete datasets[i];
    }
}

void LabelingNetworkSet::exportSelectedMIDs(QTextStream &qout)
{
    std::stringstream out;
    std::string sep = ",";
    std::string quote = "\"";

    // header
    out<<"Metabolite"<<sep<<"RI"<<sep<<"Ions (M0)"<<sep<<"M";
    // header abundances
    for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
        std::string t = datasets[ds]->getSettings().experiment;
        out <<sep<<quote<<t<<quote;
    }
    // header confidence intervals
    for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
        std::string t = "CI " + datasets[ds]->getSettings().experiment;
        out <<sep<<quote<<t<<quote;
    }
    out<<std::endl;

    // data
    for(int n = 0; n < nodes.size(); ++n) { // foreach compound
        NodeCompound *nc = nodes[n];

       // determine max. MID length for this compound
        int midLen = 0;
        for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
            std::string t = datasets[ds]->getSettings().experiment;
            if(nc->hasDataForExperiment(t)) {
                midLen = std::max(midLen, (int) (nc->getSelectedMID(t).size()));
            }
        }

        double ri = nc->getAverageRetentionIndex();


        // list of selected ions m/z
        std::stringstream ions;
        for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
            std::string t = datasets[ds]->getSettings().experiment;
            double ion = 0;

            if(nc->hasDataForExperiment(t))
                ion = nc->getSelectedIon(t);

            ions << ion << " ";
        }

        // abundances and confidence interval
        for(int i = 0; i < midLen; ++i) { // for each M+x
            out<<quote<<nc->getCompoundName()<<quote<<sep<<ri<<sep<<ions.str()<<sep<<i;

            // abundance
            for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
                std::string t = datasets[ds]->getSettings().experiment;
                std::vector<double> mid;

                if(nc->hasDataForExperiment(t))
                    mid = nc->getSelectedMID(t);

                out << sep;
                out << (mid.size() > i ? mid[i] : 0);
            }

            // confidence interval
            for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
                std::string t = datasets[ds]->getSettings().experiment;
                std::vector<double> ci;

                if(nc->hasDataForExperiment(t))
                    ci = nc->getSelectedCI(t);

                out << sep;
                out << (ci.size() > i ? ci[i] : 0);
            }
            out<<std::endl;
        }
    }

    qout<<out.str().c_str();
}

void LabelingNetworkSet::exportAllMIDs(QTextStream &qout)
{
    std::stringstream out;
    std::string sep = ",";
    std::string quote = "\"";

    // header
    out<<"Metabolite"<<sep<<"RI"<<sep<<"Ions (M0)"<<sep<<"M";
    // header abundances
    for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
        std::string t = datasets[ds]->getSettings().experiment;
        out <<sep<<quote<<t<<quote;
    }
    // header confidence intervals
    for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
        std::string t = "CI " + datasets[ds]->getSettings().experiment;
        out <<sep<<quote<<t<<quote;
    }
    out<<std::endl;

    // data
    for(int n = 0; n < nodes.size(); ++n) { // foreach compound
        NodeCompound *nc = nodes[n];

        std::set<int> lions = nc->getAllLabeledIons();
        for(std::set<int>::iterator it = lions.begin(); it != lions.end(); ++it) {
            int ion = *it;

            // determine max. MID length for this ion
            int midLen = 0;
            for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
                std::string t = datasets[ds]->getSettings().experiment;
                if(nc->hasDataForExperiment(t)) {
                    // index of ion
                    const std::vector<float> tmpIons = nc->getLabeledCompound(t)->getLabeledIons();
                    for(int i = 0; i < tmpIons.size(); ++i) {
                        if((int)tmpIons[i] == ion) {
                            midLen = std::max(midLen, (int) (nc->getLabeledCompound(t)->getIsotopomers()[i].size()));
                            break;
                        }
                    }
                }
            }
            double ri = nc->getAverageRetentionIndex();

            // abundances and confidence interval
            for(int i = 0; i < midLen; ++i) { // for each M+x
                out<<quote<<nc->getCompoundName()<<quote<<sep<<ri<<sep<<ion<<sep<<i;

                // abundance
                for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
                    std::string t = datasets[ds]->getSettings().experiment;
                    std::vector<double> mid;

                    if(nc->hasDataForExperiment(t)) {
                        // index of ion
                        const std::vector<float> tmpIons = nc->getLabeledCompound(t)->getLabeledIons();

                        for(int j = 0; j < tmpIons.size(); ++j) {

                            if((int)tmpIons[j] == ion) {
                                mid = nc->getLabeledCompound(t)->getIsotopomers()[j];
                                break;
                            }
                        }
                    }

                    out << sep;
                    out << (mid.size() > i ? mid[i] : 0);
                }

                // confidence interval
                for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
                    std::string t = datasets[ds]->getSettings().experiment;
                    std::vector<double> ci;

                    if(nc->hasDataForExperiment(t)) {
                        // index of ion
                        const std::vector<float> tmpIons = nc->getLabeledCompound(t)->getLabeledIons();
                        for(int j = 0; j < tmpIons.size(); ++j) {
                            if((int)tmpIons[j] == ion) {
                                ci = nc->getLabeledCompound(t)->getConfidenceIntervals()[j];
                                break;
                            }
                        }
                    }

                    out << sep;
                    out << (ci.size() > i ? ci[i] : 0);
                }
                out<<std::endl;
            }
        }
    }

    qout<<out.str().c_str();
}


bool LabelingNetworkSet::nodeHasEdges(int n)
{
    bool connected = false;

    for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment
        // Include this layer?
        if(!datasets[ds]->isVisible())
            continue;

        std::string t = datasets[ds]->getSettings().experiment;
        std::vector<std::vector<double> > dists = distMats[t];

        for(int j = n + 1; j < dists[n].size(); ++j) { // second node
            if(!std::isnan(dists[n][j]) && dists[n][j] <= datasets[ds]->getSettings().mid_distance_cutoff) {
                connected = true;
                break;
            }
        }

        for(int i = n - 1; i >= 0; --i) { // second node
            if(!std::isnan(dists[i][n]) && dists[i][n] <= datasets[ds]->getSettings().mid_distance_cutoff) {
                connected = true;
                break;
            }
        }

    }

    return connected;
}


/**
  Create distance matrix map with the different tracers from the nodes vector
*/

void LabelingNetworkSet::createDistanceMatrices()
{
    std::cout<<"Creating distance matrics... number of nodes: "<< nodes.size()<<std::endl;

    for(int ds = 0; ds < datasets.size(); ++ds) {
        if(ds == 0 || datasets[ds]->getSettings().nw_gap_penalty != datasets[ds - 1]->getSettings().nw_gap_penalty) {
            if(ds) {
                delete LabelingNetworkSet::distCalc;
            }
            LabelingNetworkSet::distCalc = new MIDDistanceCalculator(datasets[ds]->getSettings().nw_gap_penalty); // TODO reuse if same gap penalty
        }

        std::string t = datasets[ds]->getSettings().experiment;
        std::cout<<"### t = "<<t<<"###\n";

        // do scoring
        // keep some stats on distances:
        double dMin = 0;
        double dMax = 0;
        double dSum = 0;

        std::vector<std::vector<double> > dists;
        // calc needleman-wunsch scores
        dists.resize(nodes.size());

        int n1 = 0;
        for(QMap<int, NodeCompound*>::iterator it1 = nodes.begin(); it1 != nodes.end(); ++it1) {
            dists[n1].resize(dists.size());

            NodeCompound *node1 = it1.value();
            std::vector<double> mid1 = std::vector<double>();
            if(node1->hasDataForExperiment(t)) {
                mid1 = node1->getSelectedMID(t);
            }

            int n2 = 0;
            for(QMap<int, NodeCompound*>::iterator it2 = nodes.begin(); it2 != nodes.end(); ++it2) {
                NodeCompound *node2 = it2.value();
                std::vector<double> mid2 = std::vector<double>();
                if(node2->hasDataForExperiment(t)) {
                    mid2 = node2->getSelectedMID(t);
                }
                if(n1 == n2) { // diagonal
                    dists[n1][n2] = 1;
                } else if (n1 < n2) { // do only upper half
                    if(mid1.size() && mid2.size()) {
                        double dist = getDistance(mid1, mid2, excludeM0);

#ifdef MIAMAINWINDOW_H_ENABLE_ZSCORE
                        if(useZScore->checkState() == Qt::Checked) {
                            dist = distCalc->getMonteCarloZScore(dist, mid1.size(), mid2.size());
                        }
#endif

                        dists[n1][n2] = dist;

                        // stats
                        dMin = dist < dMin ? dist : dMin;
                        dMax = dist > dMax ? dist : dMax;
                        dSum += dist;
                    } else {
                        dists[n1][n2] = std::numeric_limits<double>::infinity();
                    }
                }
                ++n2;
            }
            ++n1;
        }
        double dMean = dSum / (dists.size() * (dists.size() - 1));
        std::cout<<"Using "<<distCalc->distanceMeasure<<" / "<<distCalc->distanceNormalization<<std::endl;
        std::cout<< "Distances ("<<dists.size()<<")\n\tRange: "<<dMin<<" - "<<dMax<<"\n\tMean: "<<dMean<<"\n";

        distMats[t]= dists;
        distRanges[t] = std::pair<double, double>(dMin, dMax);
    }

    std::cout<<"Done creating distance matrics..."<<std::endl;
}

void mia::LabelingNetworkSet::matchCompoundsAcrossExperiments(double mylibScoreCutoff, bool useLargestCommonIon)
{
    int cmpID = 0; // set as feature(COMPOUND_GROUPING_FEATURE) to make identifiable; this is also index in bigmap

    qDeleteAll(nodes);
    nodes.clear();

    gcms::GCMSSettings::LS_USE_RI        = true;
    gcms::GCMSSettings::LS_RI_DIFF       = CMP_MATCHING_RI_TOL; // TODO to config dialog TODO use also for labelidentificator!

    gcms::LibrarySearch<int,float> mylib;

    for(int tracerID = 0; tracerID < datasets.size(); ++tracerID){

        mia::NetworkLayer *ds = datasets[tracerID];

        for(std::vector<labid::LabeledCompound*>::iterator it = ds->cmpLab.begin(); it != ds->cmpLab.end(); ++it) {

            labid::LabeledCompound* lc = *it;
            NodeCompound::filterMIDs(*lc, ds->getSettings());

            if(!lc->getLabeledIons().size())
                continue; // skip if no proper label detected

            // check exclude lib
            if(ds->isExcludedMetabolite(*lc)) {
                std::cout<<"Excludelib match: RI"<<lc->getRetentionTime()<<std::endl;
                continue;
            }

            const gcms::LibraryCompound<int, float>* matchingCompound = 0;

            if(tracerID > 0) {
                // check if compounds already there, otherwise add to library
                std::vector<gcms::LibraryHit<int,float> > hits = mylib.getLibraryHits(*lc);

                if(hits.size() && hits.at(0).getOverallScore() >=  mylibScoreCutoff) {
                    // already there, make association
                    matchingCompound = hits.at(0).getLibraryCompound();
                    int prevID = atoi(matchingCompound->getFeature(COMPOUND_GROUPING_FEATURE).c_str());
                    NodeCompound* prevNode = nodes[prevID];
                    // but check that no other association made before for this experiment
                    // TODO: first check if others match better?
                    if(prevNode->hasDataForExperiment(ds->getSettings().experiment)) {
                        matchingCompound = 0; // add as new
                    } else {
                        prevNode->addLabeledCompound(ds->getSettings().experiment, lc);
                    }
                }
            }

            if(!matchingCompound) {
                // first tracer or not yet in library ->  add new

                // set unique ID as feature
                std::stringstream s;
                s<<cmpID;
                lc->addFeature(COMPOUND_GROUPING_FEATURE, s.str().c_str());

                mylib.addCompound(*lc, ds->getSettings().experiment);

                std::string compoundName = lc->getName();
                nodes[cmpID] = new NodeCompound(compoundName);
                nodes[cmpID]->setUseLargestCommonIon(useLargestCommonIon);
                nodes[cmpID]->addLabeledCompound(ds->getSettings().experiment, lc);
                nodes[cmpID]->addFeature(COMPOUND_GROUPING_FEATURE, s.str());
                ++cmpID;
            }
        }
    }
    std::cout<<"Detected "<<mylib.getCompounds().size()<<" different labeled compounds."<<std::endl;

    redetectAllIons();
    filterAndReIndexNodeCompounds();
}

void LabelingNetworkSet::redetectAllIons()
{
    // Recalculate labeling for all fragments detected somewhere
    if(datasets.size() < 2)
        return;

    QMap<std::string, int> datasetMap;
    for(int i = 0; i < datasets.size(); ++i) {
        datasetMap[datasets[i]->getSettings().experiment] = i;
    }

    foreach (NodeCompound *nc, nodes) {
        nc->redetectFragments();
        std::vector<std::string> exps = nc->getExperiments();

        for(std::vector<std::string>::iterator it = exps.begin(); it != exps.end(); ++it) {
            labid::LabeledCompound *lc = nc->getLabeledCompound(*it);
            NodeCompound::filterMIDs(*lc, datasets[datasetMap[*it]]->getSettings());
            if(!lc->getLabeledIons().size()) {
                nc->removeLabeledCompound(lc);
            }
        }
    }

}

void LabelingNetworkSet::filterAndReIndexNodeCompounds()
{
    std::cout<<"Filtering..."<<std::endl;
    // remove "empty" nodecompounds
    QMap<int, NodeCompound*> nodesOld = nodes;
    nodes.clear();

    int i = 0; // assign continuous keys to be compatible with distance matrix indexes. Why used map in the first place?
    for(QMap<int, NodeCompound*>::iterator it = nodesOld.begin(); it != nodesOld.end(); ++it) {

        if((*it)->getExperiments().size()) {
            nodes[i] = it.value();

            std::vector<std::string> exps = nodes[i]->getExperiments();

            for(int e = 0; e < exps.size(); ++e) {
                 labid::LabeledCompound *lc = nodes[i]->getCompound(exps[e]);
                 lc->addFeature(COMPOUND_GROUPING_FEATURE, std::to_string(i));
            }
            nodes[i]->addFeature(COMPOUND_GROUPING_FEATURE, std::to_string(i));
            i++;
        } else {
            std::cout<<"Remove compound with no labeled fragments: "<< (it.value())->getCompoundName()<<std::endl;
            delete it.value();
        }

    }
}

int LabelingNetworkSet::getNumberOfEdges(double variationCutoff, int excludeIfFoundInLessExperiments)
{
    // count edges
    int e = 0;
    for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment

        if(!datasets[ds]->isVisible())
            continue;

        std::string t = datasets[ds]->getSettings().experiment;
        std::vector<std::vector<double> > dists = distMats[t];
        for(int i = 0; i < dists.size(); ++i) { // first node

            if(excludeIfFoundInLessExperiments > 1 && nodes[i]->getExperiments().size() < excludeIfFoundInLessExperiments)
                continue;
            if(datasets.size() > 1 && nodes[i]->getMaxIsotopomerSD() < variationCutoff)
                continue;

            for(int j = i + 1; j < dists[i].size(); ++j) { // second node
                if(excludeIfFoundInLessExperiments > 1 && nodes[j]->getExperiments().size() < excludeIfFoundInLessExperiments)
                    continue;
                if(datasets.size() > 1 && nodes[j]->getMaxIsotopomerSD() < variationCutoff)
                    continue;
                if(dists[i][j] <= datasets[ds]->getSettings().mid_distance_cutoff)
                    ++e;
            }
        }
    }

    return e;
}

void LabelingNetworkSet::matchCompoundsAgainstLibrary(QString libFile, bool overwriteNames)
{
    if(!datasets.size() || ! nodes.size())
        return; // nothing to do

    // Identify compounds from library and set names to @node and @LabeledCompound
    std::cout<<"Identifying "<<std::endl;

    Settings s = datasets[0]->getSettings();
    if(libFile.length())
        s.cmp_id_library = libFile.toStdString();

    // gcms settings for identification (mapCompoundsToLibrary overrides LS_CUTOFF
    gcms::GCMSSettings::LS_USE_RI        = s.cmp_id_use_ri;
    gcms::GCMSSettings::LS_RI_DIFF       = s.cmp_id_ri_tol;
    gcms::GCMSSettings::LS_PURE_FACTOR   = s.gcms_pure_factor;
    gcms::GCMSSettings::LS_IMPURE_FACTOR = s.gcms_impure_factor;
    gcms::GCMSSettings::LS_MASS_FILTER   = s.cmp_id_mass_filter;

    // Library matching
    gcms::LibrarySearch<int,float> *lib;
    std::cout<<"Loading library: "<<s.cmp_id_library<<" ... ";
    try {
        lib = gcms::LibrarySearch<int,float>::fromDisk(s.cmp_id_library.c_str());
        std::cout<<"done."<<std::endl;
    } catch (...) {
        lib = new gcms::LibrarySearch<int,float>();
        std::cout<<"failed."<<std::endl;
    }

    // match first compound against library, set name to all others
    foreach (NodeCompound* nc, nodes.values()) {
        std::vector<std::string> exps = nc->getExperiments();
        std::string label;

        for(int i = 0; i < exps.size(); ++i) {
            if(i == 0) {
                label = generateCompoundLabel(nc->getLabeledCompound(exps[i]), lib, s);

                // keep old label if no match in new lib
                if(!overwriteNames && label.substr(0, 8) == "UNIDENTIFIED")
                    break;

                nc->setCompoundName(label);
                // copy new features from identification
                nc->addFeature("PRECURSOR_KEGG_ID", nc->getCompound(exps[i])->getFeature("PRECURSOR_KEGG_ID"));
            }
            nc->getCompound(exps[i])->setName(label);
        }
    }

    delete lib;

    std::cout<<"Done identifying."<<std::endl;
}


std::string LabelingNetworkSet::generateCompoundLabel(gcms::Compound<int,float> *cmp, gcms::LibrarySearch<int,float> *lib, Settings const &settings)
{
    gcms::GCMSSettings::LS_USE_RI = gcms::GCMSSettings::LS_USE_RI && (cmp->getRetentionIndex() > 0);

    std::vector<gcms::LibraryHit<int,float> > hits = lib->getLibraryHits(*cmp);

    std::stringstream ss;

    int hitCount = settings.labels_max_hits; // how many possible identifications are already appended?

    for (std::vector<gcms::LibraryHit<int,float> >::const_iterator hitIt = hits.begin(); hitIt != hits.end() && hitCount--; ++hitIt) {

        if(hitIt->getOverallScore() < settings.cmp_id_score_cutoff)
            break;

        // copy features from best hit (need e.g. KEGG ids later)
        if(hitIt == hits.begin()) {
            // preserve previously set features
            //std::map<std::string, std::string> libFeatures = hitIt->getLibraryCompound()->getFeatures();
            //libFeatures.insert(cmp->getFeatures().begin(), cmp->getFeatures().end());
            // overwrite old features
            std::map<std::string, std::string> libFeatures = cmp->getFeatures();
            libFeatures.insert(hitIt->getLibraryCompound()->getFeatures().begin(), hitIt->getLibraryCompound()->getFeatures().end());
            cmp->setFeatures(libFeatures);
        }

        ss.precision(4);
        ss<<hitIt->getLibraryCompound()->getName()<<"("<<hitIt->getOverallScore()<<")";
    }

    if(!ss.str().length()){
        ss.precision(6);
        if(cmp->getRetentionIndex() > 0)
            ss<<"UNIDENTIFIED RI "<<cmp->getRetentionIndex();
        else
            ss<<"UNIDENTIFIED RT "<<cmp->getRetentionTime() / 1000.0 / 60.0; // RT in ms->min
    }

    return ss.str();
}


std::vector<LabelingDatasetEdge *> LabelingNetworkSet::getEdges(int excludeIfFoundInLessExperiments, double variationCutoff)
{
    std::vector<LabelingDatasetEdge *> edges;

    for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment

        // Include this layer?
        if(!datasets[ds]->isVisible())
            continue;

        std::string t = datasets[ds]->getSettings().experiment;
        std::vector<std::vector<double> > dists = distMats[t];

        for(int i = 0; i < dists.size(); ++i) { // first node

            if(excludeIfFoundInLessExperiments > 1 && nodes[i]->getExperiments().size() < excludeIfFoundInLessExperiments) {
                continue;
            }

            if(datasets.size() > 1 && nodes[i]->getMaxIsotopomerSD() < variationCutoff)
                continue;


            for(int j = i + 1; j < dists[i].size(); ++j) { // second node

                if(excludeIfFoundInLessExperiments > 1 && nodes[j]->getExperiments().size() < excludeIfFoundInLessExperiments) {
                    continue;
                }

                if(datasets.size() > 1 && nodes[j]->getMaxIsotopomerSD() < variationCutoff)
                    continue;

                if(std::isnan(dists[i][j]) || dists[i][j] > datasets[ds]->getSettings().mid_distance_cutoff)
                    continue;

                LabelingDatasetEdge *e = new LabelingDatasetEdge();
                e->datasetIndex = ds;
                e->node1 = nodes[i];
                e->node2 = nodes[j];
                e->distance = dists[i][j];
                edges.push_back(e);
            }
        }
    }
    return edges;
}

std::map<int, NodeCompound *> LabelingNetworkSet::getNodesInGraph(bool showUnconnectedNodes, bool hideLessVarying, double variationCutoff, bool hideFoundInLessExperiments, int excludeIfFoundInLessExperiments)
{
    std::map<int, NodeCompound *> visNodes;

    for(int n = 0; n < nodes.size(); ++n) {

        if(!showUnconnectedNodes && !nodeHasEdges(n)) {
            continue; // not show if node has no edges
        }

        NodeCompound *nc = nodes[n];

        double variation = nc->getMaxIsotopomerSD();

        if(datasets.size() > 1 && hideLessVarying && variation < variationCutoff)
            continue;

        if(datasets.size() > 1 && hideFoundInLessExperiments && nodes[n]->getExperiments().size() < excludeIfFoundInLessExperiments)
            continue;

        visNodes[n] = nc;
    }

    return visNodes;
}


void LabelingNetworkSet::getMinMaxDistances(double &overallMin, double &overallMax)
{
    // collect distance info for edge line scaling
    overallMin = std::numeric_limits<double>::max();
    overallMax = 0;

    std::vector<std::vector<double> > minDistances(distMats[datasets[0]->getSettings().experiment]);
    std::vector<std::vector<double> > maxDistances(distMats[datasets[0]->getSettings().experiment]);
    for(int i = 0; i < minDistances.size(); ++i) {
        minDistances[i] = std::vector<double>(minDistances.size(), std::numeric_limits<double>::max());
        maxDistances[i] = std::vector<double>(minDistances.size(), 0);
    }

    for(int ds = 0; ds < datasets.size(); ++ds) { // each experiment

        // Include this layer?
        if(!datasets[ds]->isVisible())
            continue;

        std::string t = datasets[ds]->getSettings().experiment;
        std::vector<std::vector<double> > &dists = distMats[t];

        for(int i = 0; i < dists.size(); ++i) { // first node

            for(int j = i + 1; j < dists[i].size(); ++j) { // second node

                if(std::isnan(dists[i][j]) || dists[i][j] > datasets[ds]->getSettings().mid_distance_cutoff)
                    continue;

                minDistances[i][j] = std::min(minDistances[i][j], dists[i][j]);
                maxDistances[i][j] = std::max(maxDistances[i][j], dists[i][j]);

                overallMin = std::min(overallMin, dists[i][j]);
                overallMax = std::max(overallMax, dists[i][j]);
            }
        }
    }

}

QMap<int, NodeCompound *> LabelingNetworkSet::getNodeCompounds() const
{
    return nodes;
}

void LabelingNetworkSet::setUseLargestCommonIon(bool newUseCommonIon)
{
    // Need to tell NodeCompounds
    foreach (NodeCompound* c, nodes.values()) {
        c->setUseLargestCommonIon(newUseCommonIon);
    }
}

QList<NetworkLayer *> LabelingNetworkSet::getDatasets()
{
    return datasets;
}

NetworkLayer *LabelingNetworkSet::getDataset(int idx)
{
    return datasets[idx];
}

void LabelingNetworkSet::addDataset(NetworkLayer *ds)
{
    datasets.push_back(ds);
}

int LabelingNetworkSet::getSize()
{
    return datasets.size();
}

void LabelingNetworkSet::removeDataset(NetworkLayer *ds)
{
    datasets.removeOne(ds);
    createDistanceMatrices();
}

void LabelingNetworkSet::removeAllDatasets()
{
    qDeleteAll(datasets);
    datasets.clear();
}

void LabelingNetworkSet::setDistanceCutoff(double cutoff)
{
    for(int ds = 0; ds < datasets.size(); ++ds) {
        Settings s = datasets[ds]->getSettings();
        s.mid_distance_cutoff = cutoff;
        datasets[ds]->setSettings(s);
    }
}

void LabelingNetworkSet::setRelativeDistanceCutoff(double cutoff)
{
    for(int ds = 0; ds < datasets.size(); ++ds) {
        Settings s = datasets[ds]->getSettings();
        std::pair<double, double> distRange = distRanges[s.experiment];
        double min = distRange.first;
        double max = distRange.second;
        double newCutoff = min + cutoff / 100.0 * (max - min);
        std::cout<<"Set '"<<s.experiment<<"' cutoff to "<<newCutoff<<" rel cutoff is "<<cutoff<<std::endl;

        s.mid_distance_cutoff = newCutoff;

        datasets[ds]->setSettings(s);
    }
}

void LabelingNetworkSet::setExcludeM0(int excludeM0)
{
    if(this->excludeM0 != excludeM0) {
        this->excludeM0 = excludeM0;
        createDistanceMatrices();
    }
}

double LabelingNetworkSet::getDistance(std::vector<double> mid1, std::vector<double> mid2, int excludeM0)
{
    switch(excludeM0) {
    case 1:
        return distCalc->getMIDDistance(std::vector<double>(&(mid1[1]), &(mid1[mid1.size()])),
                std::vector<double>(&(mid2[1]), &(mid2[mid2.size()])));
    case 2:
        return distCalc->getMIDDistance(basePeakNormalization(std::vector<double>(&(mid1[1]), &(mid1[mid1.size()]))),
                basePeakNormalization(std::vector<double>(&(mid2[1]), &(mid2[mid2.size()]))));
    case 3:
        return distCalc->getMIDDistance(sumNormalization(std::vector<double>(&(mid1[1]), &(mid1[mid1.size()]))),
                sumNormalization(std::vector<double>(&(mid2[1]), &(mid2[mid2.size()]))));
    default:
        return distCalc->getMIDDistance(mid1, mid2);
        //dist = distCalc->getMIDDistance(basePeakNormalization(mid1), basePeakNormalization(mid2));
    }
}

}
