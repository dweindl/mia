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

#ifndef NODECOMPOUND_H
#define NODECOMPOUND_H
#include<vector>
#include<iostream>
#include<map>
#include"labeledcompound.h"
#include"lispectrum.h"
#include <QDataStream>
#include <QMap>
#include "settings.h"
#include "config.h"

namespace mia {

class NodeCompound;

/**
 * @class NodeCompound
 * @brief The NodeCompound class represent one specific compound and holds LabeledCompounds from different LabelingDatasets from different tracers/conditions.
 *
 */
class NodeCompound
{

public:

    NodeCompound();
    NodeCompound(std::string name);

    void addLabeledCompound(std::string experiment, labid::LabeledCompound *lc);
    labid::LabeledCompound *getLabeledCompound(std::string experiment);
    labid::LabeledCompound *getLabeledCompound(QString experiment);
    void removeLabeledCompound(labid::LabeledCompound *);
    void addUnlabeledCompound(std::string experiment, labid::LISpectrum *ls);
    labid::LISpectrum* getUnlabeledCompound(std::string experiment);
    labid::LabeledCompound *getCompound(std::string experiment);

    std::vector<double> getSelectedMID(std::string t);
    double getSelectedIon(std::string t);
    const std::vector<double> getSelectedCI(std::string t);
    double getSelectedR2(std::string t);
    float getLargestCommonIon();
    int getSelectedIndex(std::string t);
    int getMMinusNIon(std::string t, int loss);
    double getMaxIsotopomerSD();
    double getMinANOVAPvalue();
    double getANOVAPvalueForMassIsotopomer(int m);
    double getANOVAPvalue(std::vector<double> means, std::vector<double> sds, int n);
    double getAverageRetentionIndex();

    std::vector<std::string> getExperiments(); // The different tracer names // unique!
    bool hasDataForExperiment(std::string);
    std::string getCompoundName();
    void setCompoundName(std::string name);
    void addFeature(std::string name, std::string value);
    std::string getFeature(std::string name);
    std::string toString();
    void redetectFragments();
    labid::LabeledCompound *detectFragments(labid::LISpectrum &ul_comp, labid::LISpectrum &l_comp, const std::list< std::pair < size_t , size_t > > & fragments);

    std::set<int> getAllLabeledIons();

    static void filterMIDs(labid::LabeledCompound &lc, Settings s);
    static std::vector<double> removeTailingAbundances(const std::vector<double> &mid, double threshold);

    // Serialization
    friend QDataStream &operator<<(QDataStream &out, const NodeCompound &nc);
    friend QDataStream &operator>>(QDataStream &in, NodeCompound &nc);
    friend QDataStream &operator<<(QDataStream &out, const NodeCompound *nc);
    friend QDataStream &operator>>(QDataStream &in, NodeCompound *&nc);

    template<typename T, typename U> static bool sortByFirstPairValue(const std::pair<T, U>& i, const std::pair<T, U>& j){
        return i.first < j.first;
    }

    bool getUseLargestCommonIon() const;
    void setUseLargestCommonIon(bool value);

private:
    float largestCommonIon;
    int selectedIon;
    bool useLargestCommonIon;
    double maxSD;
    QString compoundName;               /**< Compound label. */
    QList<QString> experimentsLab;      /**< List of the tracers/conditions in which this compound was identified as labeled. */
    QList<labid::LabeledCompound*> lcs; /**< The label information from all LabelingDatasets where the compound was identified as labeled. */
    QList<QString> experimentsUnlab;    /**< List of the tracers/conditions in which this compound was identified as NOT labeled. */
    QList<labid::LISpectrum*> lis;      /**< The compound from the not-label-detected LabelingDatasets. */
    QMap<QString, QString> features;    /**< Holds additional (key->value) pairs. */

    //QList<labid::LabeledCompound> lcs2;
};

}

#endif // NODECOMPOUND_H
