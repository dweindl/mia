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

#ifndef MISC_H
#define MISC_H
#define MISC_H_DEBUG 1

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "librarysearch.h"
#include "labelidentificator.h"
#include "compound.h"
#include "library.h"
#include "librarycompound.h"

std::vector<gcms::Compound<int, float>* > compoundsFromFiles(std::vector<std::string> files);

gcms::LibrarySearch<int, float> *libraryFromCompoundFiles(std::vector<std::string> files);

std::vector<gcms::LibraryCompound<int, float>* > libraryCompoundsFromFiles(std::vector<std::string> files);

//double getCustomDistance(const std::vector<double> &v, const std::vector<double> &w);
double getCustomDistance(const std::vector<double> &v, const std::vector<double> &w, const std::vector<double> &vCI = std::vector<double>(), const std::vector<double> &wCI = std::vector<double>());
double getCanberraDistance(const std::vector<double> &v, const std::vector<double> &w);
double getManhattanDistance(const std::vector<double> &v, const std::vector<double> &w);
double getEuclideanDistance(const std::vector<double> &v, const std::vector<double> &w);
double getCosineCorrelation(const std::vector<double> &v, const std::vector<double> &w);

std::vector<double> basePeakNormalization(std::vector<double> const &v);
std::vector<double> sumNormalization(std::vector<double> const &v);

template <class T> void printVec(std::vector<T> v, std::ostream &os = std::cerr) {
    os<<v.size()<<" elements: ";
    for(int i = 0; i < v.size(); ++i) {
        os<<v[i]<<"\t";
    }
    os<<"\n";
}

template<class T> T max(std::vector<std::vector<T> > vv)  {
    T max = vv[0][0];
    for(int i = 0; i < vv[0].size(); ++i) { // row
        for(int j = 0; j < vv.size(); ++j) { // col
            T m = vv[i][j];
            if(m > max) max = m;
        }
    }
    return max;
}

template<class T> T min(std::vector<std::vector<T> > vv)  {
    T min = INFINITY;
    for(int i = 0; i < vv[0].size(); ++i) { // row
        for(int j = 0; j < vv.size(); ++j) { // col
            T m = vv[i][j];
            if(m < min) min = m;
        }
    }
    return min;
}

template<class T> T minAboveZero(std::vector<std::vector<T> > vv)  {
    T min = 0;
    for(int i = 0; i < vv[0].size(); ++i) { // row
        for(int j = 0; j < vv.size(); ++j) { // col
            T m = vv[i][j];
            if(m && m < min) min = m;
        }
    }
    return min;
}

/**
 * @brief Write simple Dot-format file for graphviz visualization.
 *
 * @param mat    Network matrix. Values > 0 represent nodes
 * @param fname  Output filename
 * @param labs   Node labels
 */
void matToDot(const std::vector<std::vector<double> > &mat, std::string fname, std::vector<std::string> labs);

/**
 * @brief Write matrix to comma separated value file.
 *
 * @param mat       Data matrix
 * @param fname     Output filename
 * @param colnames  Column labels
 * @param rownames  Row labels
 */
void matToCsv(const std::vector<std::vector<double> > &mat, std::string fname, std::vector<std::string> colnames, std::vector<std::string> rownames);

/**
 * @brief Get the maximum of 3 values.
 *
 * @param v1 Value 1
 * @param v2 Value 2
 * @param v3 Value 3
 * @return T The maximum value of v1, v2, v3.
 */
template<class T> T max(T v1, T v2, T v3) {
    return (v1 >= v2)?((v1 >= v3)?v1:v3):((v2 >= v3)?v2:v3);
}

/**
 * @brief Get the minimum of 3 values.
 *
 * @param v1 Value 1
 * @param v2 Value 2
 * @param v3 Value 3
 * @return T The minimum value of v1, v2, v3.
 */
template<class T> T min(T v1, T v2, T v3) {
    return (v1 <= v2)?((v1 <= v3)?v1:v3):((v2 <= v3)?v2:v3);
}

/**
 * @brief (Dis-)similarity score for two mass-isotopomer abundances used for Needleman-Wunsch-scoring.
 *
 * @param v1 Value 1
 * @param v2 Value 2
 * @return float The score
 */
float nwScoreMID(float v1, float v2);

/**
 * @brief Output matrix to ostream
 *
 * @param mat The matrix
 * @param os The ostream (defaults to std::cerr)
 */
void printMat(const std::vector<std::vector<double> > &mat, std::ostream &os = std::cerr);
void printMat(const std::vector<std::vector<char> > &mat, std::ostream &os = std::cerr);

/**
 * @brief Print aligned sequences.
 *
 * @param mat Needleman-Wunsch-matrix
 * @param v1 Sequence 1
 * @param v2 Sequence 1
 */
template<class T>
void printAlign(const std::vector<std::vector<double> > mat, std::vector<T> v1, std::vector<T> v2) {
    // retrace best alignment, start at bottom right
    std::vector<double> alV1, alV2;
    int i = v1.size(); // matrix is v1.size + 1
    int j = v2.size();
    int gaps = 0;
    while(i > 0 || j > 0) { // i: row, v1 ; j: col, v2
            double sright, sdown, sdiag; // scores (go left, right, diag)
            sright = j > 0?mat[i][j - 1]:1000; // if border set 1000: hilariously big score, cant go this way
            sdown = i > 0?mat[i - 1][j]:1000;
            sdiag = i*j > 0?mat[i - 1][j - 1]:1000;
            double mn = min<double>(sright, sdown, sdiag);
            if(mn == sdiag) {
                alV1.push_back(v1[--i]);
                alV2.push_back(v2[--j]);
            } else if (mn == sdown) {
                alV1.push_back(v1[--i]);
                alV2.push_back(0);
                ++gaps;
            } else if (mn == sright) {
                alV1.push_back(0);
                alV2.push_back(v2[--j]);
                ++gaps;
            } else {
                std::cerr<<"severe nw() problem";
            }
    }

    // output in reverse
    std::stringstream al;
    al << std::setw(6) << std::setprecision(2) << std::resetiosflags(std::ios::right);
    std::vector<double>::reverse_iterator it = alV1.rbegin();
    while(it < alV1.rend()) al << *(it++) << "\t\t";
    al << std::endl;
    it = alV2.rbegin();
    while(it < alV2.rend()) al << *(it++) << "\t\t";
    al << std::endl;
    std::cerr<<al.str();
}

/**
 * @brief Print aligned sequences.
 *
 * @param mat Needleman-Wunsch-matrix
 * @param v1 Sequence 1
 * @param v2 Sequence 1
 */
template<class T>
void printAlignHTML(const std::vector<std::vector<double> > mat, std::vector<T> v1, std::vector<T> v2, std::ofstream &ofs) {
    // retrace best alignment, start at bottom right
    std::vector<double> alV1, alV2;
    int i = v1.size(); // matrix is v1.size + 1
    int j = v2.size();
    int gaps = 0;
    while(i > 0 || j > 0) { // i: row, v1 ; j: col, v2
            double sright, sdown, sdiag; // scores (go left, right, diag)
            sright = j > 0?mat[i][j - 1]:1000; // if border set 1000: hilariously big score, cant go this way
            sdown = i > 0?mat[i - 1][j]:1000;
            sdiag = i*j > 0?mat[i - 1][j - 1]:1000;
            double mn = min<double>(sright, sdown, sdiag);
            if(mn == sdiag) {
                alV1.push_back(v1[--i]);
                alV2.push_back(v2[--j]);
            } else if (mn == sdown) {
                alV1.push_back(v1[--i]);
                alV2.push_back(-999);
                ++gaps;
            } else if (mn == sright) {
                alV1.push_back(-999);
                alV2.push_back(v2[--j]);
                ++gaps;
            } else {
                std::cerr<<"prob printAlignHTML";
            }
    }

    // output in reverse
    ofs <<"<br><table border=1><tr>";
    ofs << std::setw(6) << std::setprecision(4) << std::resetiosflags(std::ios::right);
    std::vector<double>::reverse_iterator it = alV1.rbegin();
    while(it < alV1.rend()) ofs << "<td>"<<*(it++) << "</td>";
    ofs << "</tr><tr>"<<std::endl;
    it = alV2.rbegin();
    while(it < alV2.rend()) ofs << "<td>"<<*(it++) << "</td>";
    ofs <<"</table>"<<std::endl;
    ofs <<"Introduced "<<gaps<<" gaps";

}

/**
 * @brief Generic Needleman&Wunsch sequence alignment algorithm (Needleman, S. B. & Wunsch, C. D. A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol, 1970, 48, 443-453)
 *
 * @param v1 Sequence 1
 * @param v2 Sequence 1
 * @param gapPenalty Penalty for gap insertion
 * @return T Score for best alignment, normalized by sequence length
 */
template<class T> T nw(std::vector<T> v1, std::vector<T> v2, double gapPenalty) {
    //         j  0     1     2     3    4
    //    i       0   V2[0] V2[1] V2[2] V2[3]
    //    0  0    0    gp   2gp    3gp  4gp
    //    1 V1[0] gp
    //    2 V1[1] 2gp
    //    3 V1[2] 3gp

    // the lower the score the better
    std::vector<std::vector<double> > mat;

    // resize, fill left and top with gap penalty
    mat.resize(v1.size() + 1); // need 0
    for(int i = 0; i <= v1.size(); ++i) {
        mat[i].resize(v2.size() + 1);
        mat[i][0] = gapPenalty * i;
    }
    for(int j = 1; j <= v2.size(); ++j) {
        mat[0][j] = gapPenalty * j;
    }

    double curGapPenalty = gapPenalty;

    // do actual stuff
    for(int i = 1; i <= v1.size(); ++i) { // row
        for(int j = 1; j <= v2.size(); ++j) { // col
            double sright, sdown, sdiag; // scores

            // had gap before? make consecutive gaps less expensive
            if((j > 1 && mat[i][j-2] < mat[i-1][j-1]) || (mat[i-1][j-1] < mat[i-1][j-1])) curGapPenalty = gapPenalty / 10;
            else curGapPenalty = gapPenalty;
            // go right?
            sright = mat[i][j - 1] + curGapPenalty;

            // had gap before? make consecutive gaps less expensive
            if((mat[i-1][j-1] < mat[i-1][j-1]) || (i > 1 && mat[i-2][j-1] < mat[i-1][j-1])) curGapPenalty = gapPenalty / 10;
            else curGapPenalty = gapPenalty;
            // go down?
            sdown = mat[i - 1][j] + curGapPenalty;

            // go diagonal?
            sdiag = mat[i - 1][j - 1] + nwScoreMID(v1[i - 1], v2[j - 1]); // + 1 because of initial 0 in matrix
            mat[i][j] = min<double>(sright, sdown, sdiag);

            // remember last action
        }
    }

#if MISC_H_DEBUG > 1
    // Debug
    std::cerr<<"vecs"<<std::endl;
    std::vector<double>::const_iterator it = v1.begin();//template??
    std::cerr<<v1.size()<<": ";
    while(it < v1.end()) std::cerr<<*(it++)<<"\t\t";
    std::cerr<<std::endl;
    it = v2.begin();
    std::cerr<<v2.size()<<": ";
    while(it < v2.end()) std::cerr<<*(it++)<<"\t\t";
    std::cerr<<std::endl;
    std::cerr<<"mat"<<std::endl;
    printMat(mat);
    std::cerr<<"al"<<std::endl;
    printAlign(mat, v1, v2);
#endif //MISC_H_DEBUG
std::cerr<<std::endl;printAlign(mat, v1, v2);std::cerr<<std::endl;
    // retrace best alignment, start at bottom right
    std::vector<double> alV1, alV2;
    int i = v1.size(); // matrix is v1.size + 1
    int j = v2.size();
    int gaps = 0;
    while(i > 0 || j > 0) { // i: row, v1 ; j: col, v2
            double sright, sdown, sdiag; // scores (go left, right, diag)
            sright = j > 0?mat[i][j - 1]:1000; // if border set 1000: hilariously big score, cant go this way
            sdown = i > 0?mat[i - 1][j]:1000;
            sdiag = i*j > 0?mat[i - 1][j - 1]:1000;
            double mn = min<double>(sright, sdown, sdiag);
            if(mn == sdiag) {
                alV1.push_back(v1[--i]);
                alV2.push_back(v2[--j]);
            } else if (mn == sdown) {
                alV1.push_back(v1[--i]);
                alV2.push_back(0);
                ++gaps;
            } else if (mn == sright) {
                alV1.push_back(0);
                alV2.push_back(v2[--j]);
                ++gaps;
            } else {
                std::cerr<<"severe nw() problem";
            }
    }
    //return getEuclideanDistance(alV1, alV2);
    return getCosineCorrelation(alV1, alV2) / (alV1.size() * alV2.size());
//    return mat[v1.size()][v2.size()] / std::max(v1.size(), v2.size());
}


template<class T> T nwHTML(std::vector<T> v1, std::vector<T> v2, double gapPenalty, std::ofstream &ofs) {
    // the lower the score the better
    std::vector<std::vector<double> > mat;

    // resize, fill left and top with gap penalty
    mat.resize(v1.size() + 1); // need 0
    for(int i = 0; i <= v1.size(); ++i) {
        mat[i].resize(v2.size() + 1);
        mat[i][0] = gapPenalty * i;
    }
    for(int j = 1; j <= v2.size(); ++j) {
        mat[0][j] = gapPenalty * j;
    }


    // v2 header
    ofs<<"<table border=1><tr><td></td><td><b>0</b></td>"<<std::endl;
    std::vector<double>::const_iterator it = v2.begin();
    while(it < v2.end()) ofs<<"<td><b>"<<*(it++)<<"</b></td>";
    ofs<<"</tr>"<<std::endl;

    double curGapPenalty = gapPenalty;

    // do actual stuff
    for(int i = 0; i <= v1.size(); ++i) { // row
        ofs<<"<tr><td><b>"<<(i==0?0:v1[i - 1])<<"</b></td>";
        for(int j = 0; j <= v2.size(); ++j) { // col
            if(i*j > 0) { // i,j== 0; only for printing

            double sright, sdown, sdiag; // scores

            // had gap before? make consecutive gaps less expensive
            if((j > 1 && mat[i][j-2] < mat[i-1][j-1]) || (mat[i-1][j-1] < mat[i-1][j-1])) curGapPenalty = gapPenalty / 10;
            else curGapPenalty = gapPenalty;
            // go right?
            sright = mat[i][j - 1] + curGapPenalty;

            // had gap before? make consecutive gaps less expensive
            if((mat[i-1][j-1] < mat[i-1][j-1]) || (i > 1 && mat[i-2][j-1] < mat[i-1][j-1])) curGapPenalty = gapPenalty / 10;
            else curGapPenalty = gapPenalty;
            // go down?
            sdown = mat[i - 1][j] + curGapPenalty;

            // go diagonal?
            sdiag = mat[i - 1][j - 1] + nwScoreMID(v1[i - 1], v2[j - 1]); // - 1 because of initial 0 in matrix
            mat[i][j] = min<double>(sright, sdown, sdiag);

            }
            ofs<<"<td>"<<mat[i][j]<<"</td>";
        }
        ofs << "</tr>\n";
    }
    ofs<<"</table>\n";

    printAlignHTML(mat, v1, v2, ofs);
    // use bottom right as similarity score; normalize by mid length
    //return mat[v1.size() - 1][v2.size() - 1] / std::max(v1.size(), v2.size());
    return mat[v1.size()][v2.size()];
}
#endif // MISC_H
