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

#include <cstdlib>
#include <iomanip>
#include <assert.h>

#include <QtConcurrent>

#include "alg/ap.h"
#include "alg/statistics.h"

#include "src/misc.h"
#include "middistancecalculator.h"

namespace mia {

// initialize distance measures
MIDDistanceCalculator::DISTANCE_MEASURE MIDDistanceCalculator::distanceMeasure = MIDDistanceCalculator::D_EUCLIDEAN;
MIDDistanceCalculator::DISTANCE_NORMALIZATION MIDDistanceCalculator::distanceNormalization = MIDDistanceCalculator::DN_SUM;

double MIDDistanceCalculator::gapPenalty = 0.3;
std::map<std::pair<int, int>, std::pair<double, double> > MIDDistanceCalculator::MCModels;

MIDDistanceCalculator::MIDDistanceCalculator(double gapPenalty)
{
     gapPenalty = gapPenalty;
}

void MIDDistanceCalculator::setDistanceMeasure(MIDDistanceCalculator::DISTANCE_MEASURE d)
{
    // doesnt make sense, is static
    distanceMeasure = d;
}

double MIDDistanceCalculator::getMIDDistance(std::vector<double> mid1, std::vector<double> mid2)
{
    return getMIDDistance(nw(mid1, mid2), mid1.size(), mid2.size());
}

double MIDDistanceCalculator::getMIDDistance(std::vector<double> mid1, std::vector<double> mid2, double gapPenalty)
{
    return getMIDDistance(nw<double>(mid1, mid2, gapPenalty), mid1.size(), mid2.size());
}

double MIDDistanceCalculator::getMIDDistance(std::pair<std::vector<double>, std::vector<double> > const &aligned, size_t origSize1, size_t origSize2)
{
    const std::vector<double> &alV1 = aligned.first;
    const std::vector<double> &alV2 = aligned.second;

    // printAlignedVectors(alV1, alV2);

    double dist;

    switch(distanceMeasure) {
    case D_EUCLIDEAN:
        dist = getEuclideanDistance(alV1, alV2);
        break;
    case D_CANBERRA:
        dist = getCanberraDistance(alV1, alV2);
        break;
    case D_MANHATTAN:
        dist = getManhattanDistance(alV1, alV2);
        break;
    case D_COSINE:
        dist = getCosineCorrelation(alV1, alV2);
        break;
    case D_CUSTOM:
        dist = getCustomDistance(alV1, alV2);
        break;
    default:
        dist = getEuclideanDistance(alV1, alV2);
    }

    double divideBy = 1;

    origSize1 = origSize1?alV1.size():origSize1;
    origSize2 = origSize2?alV2.size():origSize2;

    switch(MIDDistanceCalculator::distanceNormalization) {
    case DN_MAX:
        divideBy = std::max(origSize1, origSize2);
        break;
    case DN_MIN:
        divideBy = std::min(origSize1, origSize2);
        break;
    case DN_PROD:
        divideBy = origSize1 * origSize2;
        break;
    case DN_SUM:
        divideBy = origSize1 + origSize2;
        break;
    }

    // log(alV2.size())

    return fabs(dist / divideBy);
}

/**
 * @brief MIDDistanceCalculator::createMonteCarloModel
 * @param len1
 * @param len2
 * @param size
 * @return Pair of mean and standard deviation.
 */
std::pair<double,double> MIDDistanceCalculator::createMonteCarloModel(int len1, int len2, int size)
{
    alglib::real_1d_array dists;
    dists.setlength(size);

    // Perform sampling in separate threads
    int sizePerThread = size / QThread::idealThreadCount();
    QVector<MonteCarloHelper*> threads;
    for(int t = 0; t < QThread::idealThreadCount(); ++t) {
        int num = (t == QThread::idealThreadCount())?size - sizePerThread * t : sizePerThread; // take the rest in the last round
        MonteCarloHelper *mch = new MonteCarloHelper(&dists[sizePerThread * t], num, len1, len2, gapPenalty);
        mch->start();
        threads.push_back(mch);
    }

    bool allFinished = 0; // wait for finished
    while(!allFinished && threads.size()) {
        allFinished = 0;
        for(int t = 0; t < threads.size(); ++t) {
            MonteCarloHelper *mch = threads[t];
            if(mch->isFinished()) {
                delete mch;
                threads.remove(t);
                allFinished = 1;
            } else {
                allFinished = 0;
                mch->wait();
            }
        }
        // std::cerr<<"threads "<<threads.size()<<std::endl;
    }

    for(int i = 0; i < size; ++i)
        if(dists[i]> 1)
            std::cout<<dists[i]<<std::endl;

    double mean, variance;
    double tmp1, tmp2;
    alglib::samplemoments(dists, mean, variance, tmp1, tmp2);
    // std::cout<<"MC: "<<len1<<","<<len2<<" m = "<<mean<<", sd = "<<sqrt(variance)<<std::endl;

    return std::pair<double,double>(mean, sqrt(variance));
}

void MIDDistanceCalculator::normalize(std::vector<double>::iterator itBegin, std::vector<double>::iterator itEnd, double sum)
{
    double oldSum = MIDDistanceCalculator::sum(itBegin, itEnd);
    while(itBegin != itEnd) {
        *itBegin /= oldSum * sum;
        ++itBegin;
    }
}

double MIDDistanceCalculator::sum(std::vector<double>::iterator itBegin, std::vector<double>::iterator itEnd)
{
    double sum = 0;
    while(itBegin != itEnd) {
        sum += *itBegin;
        ++itBegin;
    }
    return sum;
}

std::vector<double> MIDDistanceCalculator::getNormalizedRandomVector(int size, int to)
{
    std::vector<double> v(size);
    double sum = 0;
   // do {
        sum = 0;
        size = v.size();
        while(size) {
            v[--size] = rand();
            sum += v[size];
        }
   // } while(v[0] / sum < 0.45);
    normalize(v.begin(), v.end(), to);
    return v;
}

std::pair<std::vector<double>, std::vector<double> > MIDDistanceCalculator::nw(std::vector<double> v1, std::vector<double> v2)
{
    return MIDDistanceCalculator::nw<double>(v1, v2, gapPenalty);
}


double MIDDistanceCalculator::getMonteCarloZScore(double distance, int size1, int size2)
{
    int s1 = std::min(size1, size2);
    int s2 = std::max(size1, size2);

    // mean already calculated?
    std::map<std::pair<int, int>, std::pair<double, double> >::iterator it = MCModels.find(std::pair<int, int>(s1, s2));
    if(it == MCModels.end()) {
        // no, calculate
        MCModels[std::pair<int, int>(s1, s2)] = createMonteCarloModel(s1, s2);
    }
    std::pair<double, double> p = MCModels[std::pair<int, int>(s1, s2)];

    return (distance - p.first) / p.second; // z-score
}

/**
 * @brief Generate random MID vectors, align them and write score to given destination. See MonteCarloHelper::MonteCarloHelper.
 */
void MonteCarloHelper::run()
{
    // std::cerr<<"Starting "<<len1<<","<<len2<<" "<<&arr<<" "<<number<<std::endl;
    for(int i = 0; i < number; ++i) {
        std::vector<double> v1 = MIDDistanceCalculator::getNormalizedRandomVector(len1, 1);
        std::vector<double> v2 = MIDDistanceCalculator::getNormalizedRandomVector(len2, 1);

        arr[i] = MIDDistanceCalculator::getMIDDistance(MIDDistanceCalculator::nw<double>(v1, v2, gapPenalty), v1.size(), v2.size());
    }
    // std::cerr<<"Finished "<<len1<<","<<len2<<" "<<&arr<<" "<<number<<std::endl;
}

template<class T>
std::pair<std::vector<T>, std::vector<T> > MIDDistanceCalculator::nw(std::vector<T> v1, std::vector<T> v2, double gapPenalty) {
    //         j  0     1     2     3    4
    //    i       0   V2[0] V2[1] V2[2] V2[3]
    //    0  0    0    gp   2gp    3gp  4gp
    //    1 V1[0] gp
    //    2 V1[1] 2gp
    //    3 V1[2] 3gp

    // the lower the score the better

    std::vector<std::vector<double> > scoreMat;
    std::vector<std::vector<char> > tracebackMat;
    initNWMatrices<T>(v1.size(), v2.size(), gapPenalty, scoreMat, tracebackMat);

    double curGapPenalty = gapPenalty;

    // do actual stuff
    for(int i = 1; i <= v1.size(); ++i) { // row
        for(int j = 1; j <= v2.size(); ++j) { // col
            double sright, sdown, sdiag; // scores

            // had gap before? make consecutive gaps less expensive
            //if((j > 1 && scoreMat[i][j-2] < scoreMat[i-1][j-1]) || (scoreMat[i-1][j-1] < scoreMat[i-1][j-1])) curGapPenalty = gapPenalty / 10;
            //else curGapPenalty = gapPenalty;
            // go right?
            sright = scoreMat[i][j - 1] + curGapPenalty;

            // make tailing gaps less expensive:
            if(i == v1.size()) sright = scoreMat[i][j - 1];

            // had gap before? make consecutive gaps less expensive
            //if((scoreMat[i-1][j-1] < scoreMat[i-1][j-1]) || (i > 1 && scoreMat[i-2][j-1] < scoreMat[i-1][j-1])) curGapPenalty = gapPenalty / 10;
            //else curGapPenalty = gapPenalty;
            // go down?
            sdown = scoreMat[i - 1][j] + curGapPenalty;
            // make tailing gaps less expensive:
            if(j == v2.size()) sdown = scoreMat[i - 1][j];

            // go diagonal?
            sdiag = scoreMat[i - 1][j - 1] + nwScoreMID(v1[i - 1], v2[j - 1]); // -1 because of initial 0 in matrix

            // least expensive path?
            if(sdown < sright) {
                if(sdown <= sdiag) {
                    scoreMat[i][j] = sdown;
                    tracebackMat[i][j] = '|';
                } else {
                    scoreMat[i][j] = sdiag;
                    tracebackMat[i][j] = '\\';
                }
            } else {
                if(sright <= sdiag) {
                    scoreMat[i][j] = sright;
                    tracebackMat[i][j] = '-';
                } else {
                    scoreMat[i][j] = sdiag;
                    tracebackMat[i][j] = '\\';
                }
            }
            //printMat(scoreMat);
            //printMat(tracebackMat);
            //std::cerr<<std::endl;
        }
    }

    // retrace best alignment, start at bottom right
    std::vector<T> alV1, alV2;
    int i = v1.size(); // matrix is v1.size + 1
    int j = v2.size();
    int gaps = 0;
    while(i > 0 || j > 0) { // i: row, v1 ; j: col, v2
        switch(tracebackMat[i][j]) {
        case '\\':
            alV1.push_back(v1[--i]);
            alV2.push_back(v2[--j]);
            break;
        case '|':
            alV1.push_back(v1[--i]);
            alV2.push_back(0);
            ++gaps;
            break;
        case '-':
            alV1.push_back(0);
            alV2.push_back(v2[--j]);
            ++gaps;
            break;
        default:
            std::cerr<<"MIDDistanceCalculator::nw: severe problem.\n";
            abort();
        }
    }
    std::reverse(alV1.begin(), alV1.end());
    std::reverse(alV2.begin(), alV2.end());

    return std::pair<std::vector<T>, std::vector<T> >(alV1, alV2);
}

template<class T>
void MIDDistanceCalculator::initNWMatrices(int size1, int size2, double gapPenalty, std::vector<std::vector<T> > &scoreMat, std::vector<std::vector<char> > &tracebackMat)
{
    // resize, fill left and top with gap penalty
    scoreMat.resize(size1 + 1); // need 0
    tracebackMat.resize(scoreMat.size());
    for(int i = 0; i <= size1; ++i) {
        scoreMat[i].resize(size2 + 1);
        scoreMat[i][0] = gapPenalty * i;

        tracebackMat[i].resize(scoreMat[i].size());
        tracebackMat[i][0] = '|';
    }
    for(int j = 1; j <= size2; ++j) {
        scoreMat[0][j] = gapPenalty * j;
        tracebackMat[0][j] = '-';
    }
    tracebackMat[0][0] = 'N';
}

template<class T>
void MIDDistanceCalculator::printAlignedVectors(const std::vector<T> &v1, const std::vector<T> &v2)
{
    assert(v1.size() == v2.size());
    std::cout<<std::setw(6) << std::setprecision(3) << std::resetiosflags(std::ios::right);
    std::cout<<"Seq1: ";
    for(int i = 0; i < v1.size(); ++i) {
        std::cout<<v1[i]<<"  ";
    }
    std::cout<<std::endl;
    std::cout<<"Seq2: ";
    for(int i = 0; i < v2.size(); ++i) {
        std::cout<<v2[i]<<"  ";
    }
    std::cout<<std::endl;
}

}
