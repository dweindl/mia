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

#include <limits>
#include <algorithm>
#include <gsl/gsl_cdf.h>

#include "nodecompound.h"
#include "labeledcompound.h"
#include "labelingdataset.h"
#include "misc.h"
#include "midsolver.h"
#include "labidexception.h"
#include "serializationqt.h"

namespace mia {

NodeCompound::NodeCompound()
{
    compoundName = "";
    maxSD = -1;

    largestCommonIon = -1;
    selectedIon = -1;
}

NodeCompound::NodeCompound(std::string name)
{
    compoundName = QString::fromStdString(name);
    maxSD = -1;

    largestCommonIon = -1;
    selectedIon = -1;
}

void NodeCompound::addLabeledCompound(std::string experiment, labid::LabeledCompound *lc)
{
    //std::cerr<<"addLabeledCompound: "<<compoundName.toStdString()<<" "<<lc->getRetentionIndices()[0]<<std::endl;
    lcs.push_back(lc);
    //lcs2.push_back(labid::LabeledCompound(*lc));
    experimentsLab.push_back(QString::fromStdString(experiment));

    // copy features TODO: no-overwrite!
    std::map<std::string, std::string> fts = lc->getFeatures();
    for(std::map<std::string, std::string>::const_iterator it = fts.begin(); it != fts.end(); ++it) {
        features[QString::fromStdString(it->first)] = QString::fromStdString(it->second);
    }

    maxSD = -1;
    largestCommonIon = -1;
    selectedIon = -1;
}

labid::LabeledCompound *NodeCompound::getLabeledCompound(std::string experiment)
{
    int idx = experimentsLab.indexOf(QString::fromStdString(experiment));
    return lcs[idx];
}

labid::LabeledCompound *NodeCompound::getLabeledCompound(QString experiment)
{
    int idx = experimentsLab.indexOf(experiment);
    return lcs[idx];
}

void NodeCompound::removeLabeledCompound(labid::LabeledCompound *lc)
{
    int idx = lcs.indexOf(lc);
    lcs.removeAt(idx);
    experimentsLab.removeAt(idx);

    maxSD = -1;
    largestCommonIon = -1;
    selectedIon = -1;
}

void NodeCompound::addUnlabeledCompound(std::string experiment, labid::LISpectrum *ls)
{
    lis.push_back(ls);
    experimentsUnlab.push_back(QString::fromStdString(experiment));

    maxSD = -1;
    selectedIon = -1;
}

labid::LISpectrum *NodeCompound::getUnlabeledCompound(std::string experiment)
{
    int idx = experimentsUnlab.indexOf(QString::fromStdString(experiment));
    return lis[idx];
}

labid::LabeledCompound *NodeCompound::getCompound(std::string experiment)
{
    int idx = experimentsLab.indexOf(QString::fromStdString(experiment));

    if(idx > -1)
        return lcs[idx];

    idx = experimentsUnlab.indexOf(QString::fromStdString(experiment));

    return dynamic_cast<labid::LabeledCompound *>(lis[idx]);
}

std::vector<double> NodeCompound::getSelectedMID(std::string t)
{
    labid::LabeledCompound* lc = getCompound(t);
    // TODO mia::pickMIDs

    // find best ion  // TODO: consider intensity and R2
    int idx = getSelectedIndex(t);
    return lc->getIsotopomers().at(idx);
    //return lc->getMIDOfLargestFragment(0.95);
}

double NodeCompound::getSelectedIon(std::string t)
{
    labid::LabeledCompound* lc = getCompound(t);
    return lc->getLabeledIons().at(getSelectedIndex(t));
    //return lc->getLabeledIons().at(lc->getLargestSignificantIonIndex(0.95));
}

const std::vector<double> NodeCompound::getSelectedCI(std::string t)
{
    labid::LabeledCompound* lc = getCompound(t);
    return lc->getConfidenceIntervals().at(getSelectedIndex(t));
}

double NodeCompound::getSelectedR2(std::string t)
{
    labid::LabeledCompound* lc = getCompound(t);
    return lc->getR2s().at(getSelectedIndex(t));
//    return lc->getR2s().at(lc->getLargestSignificantIonIndex(0.95));
}

/**
 * @brief NodeCompound::getLargestCommonIon
 * @return highest common m/z
 */
float NodeCompound::getLargestCommonIon()
{
    if(largestCommonIon > 0) {
        return largestCommonIon;
    }

    std::vector<float> oldIons = lcs[0]->getLabeledIons(); // ions already sorted
    std::vector<float> commonIons(oldIons.size());
    foreach (labid::LabeledCompound* lc, lcs) {
        std::vector<float> curIons = lc->getLabeledIons();
        std::vector<float>::iterator isit = std::set_intersection(oldIons.begin(), oldIons.end(), curIons.begin(), curIons.end(), commonIons.begin());
        commonIons.resize(isit - commonIons.begin());
        if(!commonIons.size())
            break;
        oldIons = commonIons;
    }

    // average r^2
    std::vector<float> meanR2(commonIons.size(), 0);
    foreach (labid::LabeledCompound* lc, lcs) {
        std::vector<float> curIons = lc->getLabeledIons();

        for(int i = 0; i < commonIons.size(); ++i) {
            int curIdx = std::find(curIons.begin(), curIons.end(), commonIons[i]) - curIons.begin();
            meanR2[i] += lc->getR2s()[curIdx] / lcs.size();
        }
    }

    // choose largest > 0.95
    int idx = commonIons.size() - 1;
    for(; idx > 0; --idx) {
        if(meanR2[idx] >= 0.95)
            break;
    }

    if(idx == -1) {
        return commonIons[commonIons.size() - 1];
    }

    largestCommonIon = commonIons[idx];

    return largestCommonIon;
}

int NodeCompound::getSelectedIndex(std::string t)
{
    labid::LabeledCompound* lc = getCompound(t);

    if(selectedIon < 0) {
        if(useLargestCommonIon) {
            selectedIon = getLargestCommonIon();
        } else {
            try{
                //selectedIon = lc->getLargestSignificantIon(0.95);
                selectedIon = lc->getLabeledIons().size() - 1;
            } catch(labid::LabIDException e) {
                // bug in ntfd where R2 can be < 0, which will
                // throw and exception in getLargestSignificantIon
            }
        }
    }

    std::vector<float> ions = lc->getLabeledIons();
    for(int i = 0; i < ions.size(); ++i) {
        if(ions[i] == selectedIon) {
            return i;
        }
    }

    return ions.size() - 1;
}

/**
 * @brief Get the m/z of the [M-n] ion if exists, otherwise 0
 * @param t Experiment
 * @param loss n
 * @return
 */
int NodeCompound::getMMinusNIon(std::string t, int loss)
{
    labid::LabeledCompound* lc = getCompound(t);

    std::list<std::list<std::pair<float,int> > > dummy;
    std::list<gcms::Fragment> fragments = lc->getSpectrumFragments(dummy);

    float largestFragment;

    // check isotope clusters

    for(std::list<gcms::Fragment>::reverse_iterator it = fragments.rbegin(); it != fragments.rend(); ++it) {

        if(it == fragments.rbegin()) {
             largestFragment = it->getIon();
        } else {
            float curFrag = it->getIon();

            if(curFrag == largestFragment - loss)
                return curFrag;
            else if(curFrag < largestFragment - loss)
                break;
        }
    }

    // or see if at least M0 for higher mass peak is present
    if(lc->getIntensity(largestFragment + loss, 0.2))
        return largestFragment;

    return 0;
}

double NodeCompound::getMaxIsotopomerSD()
{
    if(maxSD >= 0)
        return maxSD;

    std::vector<std::string> exps = getExperiments();
    if(exps.size() == 1)
        return 0;

    // to adjust color to labeling variation
    int expCount = 1;
    int isoCount = 0;

    while(expCount) {
        expCount = 0;
        int i;
        double sum = 0;

        // calc mean M_n
        for(i = 0; i < exps.size(); ++i) {
            const std::vector<double> &mid = getSelectedMID(exps[i]);
            if(mid.size() > isoCount) {
                sum += mid[isoCount];
                ++expCount;
            }
        }

        double mean = sum / expCount;

        // calc SD M_n
        sum = 0;
        expCount = 0;
        for(i = 0; i < exps.size(); ++i) {
            const std::vector<double> &mid = getSelectedMID(exps[i]);
            double mi = 0;
            if(mid.size() > isoCount) {
                mi = mid[isoCount];
                ++expCount;
            }
            sum += (mi - mean) * (mi - mean);
        }

        double sd = sqrt(sum / i);

        //maxSD = std::max(maxSD, sd / mean);
        maxSD = std::max(maxSD, sd);
        ++isoCount;
    }

    return maxSD;
}

double NodeCompound::getMinANOVAPvalue()
{
    std::vector<std::string> exps = getExperiments();
    if(exps.size() == 1)
        return 1;

    // max MID length
    int maxLen = 0;
    for(int i = 0; i < exps.size(); ++i) {
        int curSize = getSelectedMID(exps[i]).size() + 1;
        maxLen = std::max(maxLen, curSize);
    }

    int minP = 1;

    for(int m = 0; m < maxLen; ++m) {
        minP = minP < getANOVAPvalueForMassIsotopomer(m) ? minP : getANOVAPvalueForMassIsotopomer(m);
    }

    return minP;
}

double NodeCompound::getANOVAPvalueForMassIsotopomer(int m)
{
    // experiment means and standard deviations
    std::vector<double> means;
    std::vector<double> sds;

    std::vector<std::string> exps = getExperiments();

    //int oldN = 0;
    int minNumFiles = 999999;

    for(int i = 0; i < exps.size(); ++i) {

        int midSize = getSelectedMID(exps[i]).size();
        if(midSize < m + 1)
            continue;

        // current group number of observations
        int N = getLabeledCompound(exps[i])->getLabeledSpecCount()
                * getLabeledCompound(exps[i])->getUnLabeledSpecCount();

        int curMinNumFiles = std::max(getLabeledCompound(exps[i])->getLabeledSpecCount(),
                                      getLabeledCompound(exps[i])->getUnLabeledSpecCount());
        minNumFiles = std::min(curMinNumFiles, minNumFiles);

        //if(oldN > 0 && N != oldN) {
            //std::cerr<<"Cannot perform ANOVA for different numbers of replicates."<<std::endl;
            //return 1;
        //}

        means.push_back(getSelectedMID(exps[i])[m]);

        int df = (N - 1) * midSize;
        double t = gsl_cdf_tdist_Pinv(0.9, df);
        double ci = getSelectedCI(exps[i])[m];
        sds.push_back(ci / t);

        //oldN = N;
    }

    return getANOVAPvalue(means, sds, minNumFiles);
}

double NodeCompound::getANOVAPvalue(std::vector<double> means, std::vector<double> sds, int n)
{
    if(n < 2)
        return 1;

    int k = means.size();

    if(k < 2)
        return 1;

    assert(k == sds.size());

    int df1 = n - 1;
    int df2 = n * k - df1 - 1;

    double meanMean = 0;
    double sumVar = 0;
    for(int i = 0; i < k; ++i) {
        meanMean += means[i];
        sumVar += sds[i] * sds[i];
    }
    meanMean /= k;

    double s_x = 0;
    for(int i = 0; i < k; ++i) {
        s_x += (means[i] - meanMean) * (means[i] - meanMean);
    }
    s_x /= k - 1;

    double f = n * s_x / (sumVar / k);

    double p = gsl_cdf_fdist_Q(f, df1, df2);

    return p;
}

double NodeCompound::getAverageRetentionIndex()
{
    double ri = 0;
    double riCount = 0;
    foreach(QString key, experimentsLab){
        ri += getLabeledCompound(key)->getRetentionIndex();
        ++riCount;
    }

    return ri /= riCount;
}

std::vector<std::string> NodeCompound::getExperiments()
{
    std::vector<std::string> v;

    foreach(QString key, experimentsLab){
      v.push_back(key.toStdString());
    }
    // include unlabeled
    foreach(QString key, experimentsUnlab){
      v.push_back(key.toStdString());
    }
    return v;
}

bool NodeCompound::hasDataForExperiment(std::string s)
{
    return experimentsLab.contains(QString::fromStdString(s));
}

std::string NodeCompound::getCompoundName()
{
    return compoundName.toStdString();
}

void NodeCompound::setCompoundName(std::string name)
{
    compoundName = QString::fromStdString(name);
}

void NodeCompound::addFeature(std::string name, std::string value)
{
    features[QString::fromStdString(name)] = QString::fromStdString(value);
}

std::string NodeCompound::getFeature(std::string name)
{
    return features[QString::fromStdString(name)].toStdString();
}

std::string NodeCompound::toString()
{
    std::stringstream ss;
    ss << "**********"<<std::endl;
    ss << "Name: " << compoundName.toStdString() << std::endl;
    ss << "Experiments: ";
    foreach(QString exp, experimentsLab) {
        ss << exp.toStdString() << "\t";
    }
    ss << std::endl;
    for(int i = 0; i < experimentsLab.size(); ++i) {
        labid::LabeledCompound* lc = lcs[i];
        ss<<lc->getName()<<"\t"<<*(lc->getRetentionIndices().begin())<<std::endl;
        ss<<"\tTotal signal:"<<lc->getTotalSignal()<<std::endl;
        ss<<"\t";
        printVec<float>(lc->getLabeledIons(), ss);
    }
    ss << "Features:"<<std::endl;
    for(QMap<QString, QString>::Iterator it = features.begin(); it != features.end(); ++it) {
        ss<<"\t"<<it.key().toStdString()<<"\t"<<it.value().toStdString()<<std::endl;
    }
    ss << "**********"<<std::endl;
    return ss.str();
}

void NodeCompound::redetectFragments()
{
    std::cout<<"Redetecting labeled fragments: "<<compoundName.toStdString()<<"\n";

    // find all ions labeled in any experiment
    std::list<std::pair<size_t, size_t> > frags;
    foreach(labid::LabeledCompound *lc, lcs) {
        const std::vector<float> ions = lc->getLabeledIons();
        const std::vector<std::vector<double> > mids = lc->getIsotopomers();
        for(int i = 0; i < ions.size(); ++i) {
            size_t lower = ions[i];
            size_t upper = lower + mids[i].size();

            // if not exists, add to list (?or alter range)
            std::list<std::pair<size_t, size_t> >::iterator it;
            for(it = frags.begin(); it != frags.end(); ++it) {
                if(it->first == lower) {
                    it->second = std::max(upper, it->second);
                    break;
                }
            }
            if(it == frags.end()) {
                frags.push_back(std::pair<size_t, size_t>(lower, upper));
            }
        }
    }

    frags.sort(sortByFirstPairValue<int,int>);

    // if not present redetect
    for(QList<labid::LabeledCompound*>::iterator it = lcs.begin(); it != lcs.end(); ++it){
        //std::cout<<"Redetect "<<compoundName.toStdString()<<std::endl;
        labid::LabeledCompound lc = labid::LabeledCompound(**it);
        labid::LISpectrum lis = labid::LISpectrum(*lc.getLabeledReferenceCompound());
        *it = detectFragments(lc, lis, frags); // delete old lc?

        // copy features
        (*it)->setFeatures(lc.getFeatures());
    }

    // TODO: also in "unlabeled" compounds
    for(QList<labid::LISpectrum*>::iterator it = lis.begin(); it != lis.end(); ++it){
        //std::cout<<"Redetect "<<compoundName.toStdString()<<std::endl;
        labid::LabeledCompound lc = labid::LabeledCompound(*(dynamic_cast<labid::LabeledCompound *>(*it)));
        labid::LISpectrum lis = labid::LISpectrum(*lc.getLabeledReferenceCompound());
        *it = detectFragments(lc, lis, frags); // delete old lc?

        // copy features
        (*it)->setFeatures(lc.getFeatures());
    }

}

labid::LabeledCompound* NodeCompound::detectFragments(labid::LISpectrum &ul_comp, labid::LISpectrum &l_comp, const std::list<std::pair<size_t, size_t> > &fragments)
{
    // Code from labid::LabelIdentificator::analyseFragments
    //*it = labid::LabelIdentificator::analyseFragments(**it, *((*it)->getLabeledReferenceCompound()), frags, 0); // can make static?? no, notwork

    //init specs
    std::list<std::vector<double> > ul_specs = ul_comp.getNormalizedSourceSpectra();
    std::list<std::vector<double> > l_specs = l_comp.getNormalizedSourceSpectra();

    //init result containers
    std::vector<float> labeled_ions;
    std::vector<std::vector<double> > isotopomers;
    std::vector<std::vector<double> > confidence;
    std::vector<double> scores;
    std::vector<double> r2s;

    //iterate over all fragments
    for ( std::list<std::pair<size_t,size_t> >::const_iterator it=fragments.begin();it!=fragments.end();++it )
    {


        size_t start=it->first;
        size_t end=it->second;

        //create fragment vectors for MID calculation
        std::vector<std::vector<double> > frags_unlab;
        for(std::list<std::vector<double> >::const_iterator it2=ul_specs.begin();it2!=ul_specs.end();++it2)
        {
            const std::vector<double>& spec=*it2;

            if(spec.size() <= end)
                continue;

            std::vector<double> rc_frag(end-start+1, 0.0);
            for(size_t i=start;i<=end && i<spec.size();i++)
            {
                rc_frag.at(i-start)=spec.at(i);
            }
            frags_unlab.push_back(rc_frag);

        }

        std::vector<std::vector<double> > frags_lab;
        for(std::list<std::vector<double> >::const_iterator it2=l_specs.begin();it2!=l_specs.end();++it2)
        {
            const std::vector<double>& spec=*it2;

            if(spec.size() <= end)
                continue;

            std::vector<double> rc_frag(end-start+1, 0.0);
            for(size_t i=start;i<=end && i<spec.size();i++)
            {
                rc_frag.at(i-start)=spec.at(i);
            }
            frags_lab.push_back(rc_frag);

        }

        //solve system
        tools::MIDSolver mid_solv(frags_lab, frags_unlab, 0.010934);

        double r2=mid_solv.getR2();
        // 			cout<<"R^2: "<<r2<<std::endl;

        //calculate confidence intervals of coefficients
        std::vector<double> cis=mid_solv.getCIs();

        //calculate sum and create result list
        std::vector<double> mids=mid_solv.getMIDs();
        double sum=0.0;
        double label_amount=0.0;
        double m0=mids.at(0);
        for ( size_t n=0;n<mids.size();n++ )
        {
            double value=mids.at(n);
            sum+=std::fabs ( value );
            if ( n>0 && value>0 )
            {
                label_amount+=value;
            }
        }

        bool conf_ok=false;
        std::vector<double>::const_iterator it_c=cis.begin();
        for ( std::vector<double>::const_iterator it_s=mids.begin();it_s!=mids.end();++it_s, it_c++ )
        {
            const double& rc_val=*it_s;
            const double& rc_conf=*it_c;

            if ( it_c==cis.begin() )
            {
                continue; //skip M0 isotopomer
            }
            if ( rc_val>=0 && ( rc_conf==-1 || rc_val-rc_conf>0 || rc_val+rc_conf<0 ) )
            {
                conf_ok=true;
            }
        }

        double max_label = 0.55; // DW
        max_label = 1;
        //if ( r2>=min_r2 && sum>=1.0-sum_thr && sum<=1.0+sum_thr && label_amount>=min_label && m0<=1.0-min_label &&  conf_ok && m0 >= 1 - max_label)
        {
            l_comp.addFragment(start, end);
            ul_comp.addFragment(start, end);
            labeled_ions.push_back ( start );
            isotopomers.push_back ( mids  );
            confidence.push_back ( cis );
            scores.push_back ( sum  );
            r2s.push_back (r2  );
        }
    }

    labid::LabeledCompound* result=new labid::LabeledCompound( &l_comp, &ul_comp );
    result->setLabeledIons ( labeled_ions );
    result->setIsotopomers ( isotopomers );
    result->setConfidenceIntervals ( confidence );
    result->setScores ( scores );
    result->setR2s ( r2s );
    return result;
}

std::set<int> NodeCompound::getAllLabeledIons()
{
    std::set<int> s;

    for(int i = 0; i < lcs.size(); ++i) {
        labid::LabeledCompound* lc = lcs[i];
        const std::vector< float > ions = lc->getLabeledIons();
        for(int j = 0; j < ions.size(); ++j) {
            s.insert(ions[j]);
        }
    }

    return s;
}

void NodeCompound::filterMIDs(labid::LabeledCompound &lc, Settings s)
{
    // Remove ions with too high mass isotopomers
    std::vector<float> ions = lc.getLabeledIons();
    std::vector<std::vector<double> > mids = lc.getIsotopomers();
    std::vector<double> r2s = lc.getR2s();
    std::vector<std::vector<double> > cis = lc.getConfidenceIntervals();

    std::vector<float> ions2;
    std::vector<std::vector<double> > mids2;
    std::vector<double> r2s2;
    std::vector<std::vector<double> > cis2;
    for(int i = 0; i < ions.size(); ++i) {
        // filter r2
        if(r2s[i] < s.lid_req_r2)
            continue;

        // remove if max(mid) * tracer mass increment > m/z
        if(mids[i][mids[i].size() - 1] * s.tracer_atom_mass > ions[i]
                && fabs(mids[i][mids[i].size() - 1]) > cis[i][cis[i].size() - 1])
            continue;

        // Filter compounds with less M0 than the tracer
        if(mids[i].size() <= s.lid_max_mass_isotopomer + 1
                && mids[i][0] >= s.lid_min_m0
                && mids[i].size() * s.tracer_atom_mass <= ions[i] ) { // mass filter: make sure m/z < M_max * tracer_isotope_mass
            ions2.push_back(ions[i]);

            // remove tailing values with abundance < 0.01
            mids2.push_back(removeTailingAbundances(mids[i], 0.01)); // TODO add to settings
            //mids2.push_back(mids[i]);
            r2s2.push_back(r2s[i]);
            cis2.push_back(cis[i]);
            cis2[cis2.size() - 1].resize(mids2[mids2.size() - 1].size());
        }

        // filter sum
        double absSum = 0;
        for(std::vector<double>::iterator miIt = mids[i].begin(); miIt != mids[i].end(); ++miIt)
            absSum += std::fabs(*miIt);
        if(absSum - 1 > s.lid_maximal_frag_dev)
            continue;
    }
    lc.setLabeledIons(ions2);
    lc.setIsotopomers(mids2);
    lc.setConfidenceIntervals(cis2);
    lc.setR2s(r2s2);
}

std::vector<double> NodeCompound::removeTailingAbundances(const std::vector<double> &mid, double threshold)
{
    bool tail = 1; // still in tail (nothing > threshold yet)
    std::vector<double> res;
    for(std::vector<double>::const_reverse_iterator it = mid.rbegin(); it < mid.rend(); ++it) {
        if(!tail || *it > threshold) {
            tail = 0;
            res.push_back(*it);
        }
    }
    std::reverse(res.begin(), res.end());

    return res;
}
bool NodeCompound::getUseLargestCommonIon() const
{
    return useLargestCommonIon;
}

void NodeCompound::setUseLargestCommonIon(bool value)
{
    useLargestCommonIon = value;

    selectedIon = -1;
}

QDataStream &operator <<(QDataStream &out, const NodeCompound &nc)
{
    out<<QString("NodeCompound");

    static const uint version = 1; // TODO

    out << nc.compoundName;
    out << nc.features;
    out << nc.lcs.size();
    for(int i = 0; i < nc.lcs.size(); ++i) {
        out << nc.experimentsLab[i] << nc.lcs[i];
    }

    return out;
}

QDataStream &operator >>(QDataStream &in, NodeCompound &nc)
{
    QString magic;
    in >> magic;
    Q_ASSERT(magic == "NodeCompound");

    in >> nc.compoundName;
    in >> nc.features;

    int size;
    in >> size;
    for(int i = 0; i < size; ++i) {
        QString s;
        in >> s;
        nc.experimentsLab.push_back(s);

        labid::LabeledCompound *lc;
        in >> lc;
        nc.lcs.push_back(lc);
    }

    return in;
}

QDataStream &operator <<(QDataStream &out, const NodeCompound *nc)
{
    out << *nc;
    return out;

}

QDataStream &operator >>(QDataStream &in, NodeCompound *&nc)
{
    nc = new NodeCompound;
    in >> *nc;
    return in;
}

}
