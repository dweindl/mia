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

#ifndef LABELINGNETWORKSET_H
#define LABELINGNETWORKSET_H

#include <QObject>
#include <QTextStream>
#include <map>

#include "nodecompound.h"
#include "networklayer.h"
#include "middistancecalculator.h"

namespace mia {

class LabelingDatasetEdge {
public:
    int datasetIndex;
    NodeCompound *node1;
    NodeCompound *node2;
    double distance;
};


/**
 * @brief The LabelingNetworkSet class holds several LabelingDatasets, matches all compounds across these datasets, invokes reanalysis of undetected
 * labeled fragments and manages the overlay of compound networks in its LabelingDatasets.
 */
class LabelingNetworkSet : public QObject
{
    Q_OBJECT

public:
    // TODO: refactor: subclass Labelidentificator
    // addLabelingNetwork();
    // setActive(tracer);
    // getGraph() // as dot?

    LabelingNetworkSet();
    ~LabelingNetworkSet();

    void exportSelectedMIDs(QTextStream &qout);

    void exportAllMIDs(QTextStream &qout);

    bool nodeHasEdges(int n);

    void createDistanceMatrices(); /** Setup the distance matrices */

    void matchCompoundsAcrossExperiments(double mylibScoreCutoff, bool useLargestCommonIon);

    void redetectAllIons();

    void filterAndReIndexNodeCompounds();

    int getNumberOfEdges(double variationCutoff, int excludeIfFoundInLessExperiments);

    void matchCompoundsAgainstLibrary(QString libFile, bool overwriteNames);

    std::string generateCompoundLabel(gcms::Compound<int, float> *cmp, gcms::LibrarySearch<int,float> *lib, Settings const &settings);

    std::vector<LabelingDatasetEdge *> getEdges(int excludeIfFoundInLessExperiments, double variationCutoff);

    std::map<int, NodeCompound *> getNodesInGraph(bool showUnconnectedNodes,
                                                  bool hideLessVarying, double variationCutoff,
                                                  bool hideFoundInLessExperiments, int excludeIfFoundInLessExperiments
                                                  );

    void getMinMaxDistances(double &overallMin, double &overallMax);

    QMap<int, NodeCompound*> getNodeCompounds() const;

    void setUseLargestCommonIon(bool newUseCommonIon);

    QList<NetworkLayer*> getDatasets();
    NetworkLayer* getDataset(int idx);
    void addDataset(NetworkLayer *ds);

    int getSize();

    void removeDataset(NetworkLayer *ds);
    void removeAllDatasets();

    void setDistanceCutoff(double cutoff);
    void setRelativeDistanceCutoff(double cutoff);

    void setExcludeM0(int excludeM0);

    double getDistance(std::vector<double> mid1, std::vector<double> mid2, int excludeM0);

    static MIDDistanceCalculator *distCalc;

signals:

public slots:

private:
    std::map<std::string, std::vector<std::vector<double> > > distMats; /** Distance matrices */
    QMap<int, NodeCompound*> nodes; /** All the different compounds found in any experiment, index is the ID-feature of the compound */
    QList<NetworkLayer*> datasets; /** The "raw" data from the different experiments */
    std::map<std::string, std::pair<double, double> > distRanges; /** Distance matrices (min, max) */
    int excludeM0;
};

}
#endif // LABELINGNETWORKSET_H
