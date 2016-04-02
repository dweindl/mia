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

#ifndef MIDDISTANCECALCULATOR_H
#define MIDDISTANCECALCULATOR_H

#include <vector>
#include <map>

#include <QThread>

#include "alg/ap.h"

namespace mia {

class MIDDistanceCalculator;
class MonteCarloHelper;

/**
 * @brief The MIDDistanceCalculator class does MID alignment and calculates the difference score.
 */
class MIDDistanceCalculator
{
public:

    enum DISTANCE_MEASURE {
        D_EUCLIDEAN,
        D_CANBERRA,
        D_MANHATTAN,
        D_COSINE,
        D_CUSTOM,
    };

    enum DISTANCE_NORMALIZATION{
        DN_ONE,
        DN_SUM,
        DN_PROD,
        DN_MAX,
        DN_MIN
    };

    static DISTANCE_MEASURE distanceMeasure;
    static DISTANCE_NORMALIZATION distanceNormalization;

    MIDDistanceCalculator(double gapPenalty);

    void setDistanceMeasure(DISTANCE_MEASURE d);

    double getMIDDistance(std::vector<double> mid1, std::vector<double> mid2);
    static double getMIDDistance(std::vector<double> mid1, std::vector<double> mid2, double gapPenalty);

    static double getMIDDistance(const std::pair<std::vector<double>, std::vector<double> > &aligned, size_t origSize1 = 0, size_t origSize2 = 0);

    // z-score functions need checking!
    static std::pair<double,double> createMonteCarloModel(int len1, int len2, int size = MCMsize);
    static double getMonteCarloZScore(double distance, int size1, int size2);

    static void normalize(std::vector<double>::iterator itBegin, std::vector<double>::iterator itEnd, double sum = 1);
    static double sum(std::vector<double>::iterator itBegin, std::vector<double>::iterator itEnd);
    static std::vector<double> getNormalizedRandomVector(int size, int sum = 1);

    // Needleman-Wunsch
    template<class T> static std::pair<std::vector<T>, std::vector<T> > nw(std::vector<T> v1, std::vector<T> v2, double gapPenalty);
    std::pair<std::vector<double>, std::vector<double> > nw(std::vector<double> v1, std::vector<double> v2);

    template<class T> static void initNWMatrices(int size1, int size2, double gapPenalty, std::vector<std::vector<T> > &scoreMat, std::vector<std::vector<char> > &tracebackMat);
    template<class T> static void printAlignedVectors(std::vector<T> const &v1, std::vector<T> const &v2);

private:
    static const int MCMsize = 1000; /**< Number of MID pairs to generate. */
    static double gapPenalty; /**< Gap penalty for Needleman-Wunsch-alignment. */
    static std::map<std::pair<int, int>, std::pair<double, double> > MCModels; /**< Monte-Carlo models by MID size. (size1, size2) -> (mean, standard deviation); size1 <= size2. */
};

/**
 * @brief The MonteCarloHelper class creates random mass isotopomer distrubutions for Monte-Carlo-based distance cutoff. Needed for multi-threading.
 */
class MonteCarloHelper : public QThread {

    Q_OBJECT

public:
    MonteCarloHelper(double *arr, int number, int len1, int len2, double gapPenalty)
        : arr(arr), number(number), len1(len1), len2(len2), gapPenalty(gapPenalty) {}

private:
    void run();

    MIDDistanceCalculator *parent; /**< The instanciating object. */
    double *arr;        /**< Pointer to double array to put the @b number MIDs into. */
    int number;         /**< Number of MID pairs to generate. */
    int len1;           /**< Length of first MID vector. */
    int len2;           /**< Length of second MID vector. */
    double gapPenalty;  /**< Gap penalty for Needleman-Wunsch-alignment. */
};

}
#endif // MIDDISTANCECALCULATOR_H
