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
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <assert.h>
#include <limits>

#include "librarysearch.h"
#include "labelidentificator.h"
#include "compound.h"
#include "library.h"
#include "librarycompound.h"

std::vector<gcms::Compound<int, float>* > compoundsFromFiles(std::vector<std::string> files) {
    std::vector<gcms::Compound<int, float>* > comps; // list with compounds from all files
    for(std::vector<std::string>::iterator it = files.begin(); it < files.end(); ++it) {
        std::vector<gcms::Compound<int, float>* > tmp = gcms::Compound<int, float>::fromDisk(it->c_str());
        for(int i = 0; i < tmp.size(); ++i) {
            comps.push_back(tmp.at(i));
        }
    }
    return comps;
}

double getEuclideanDistance(const std::vector<double> &v, const std::vector<double> &w) {
    assert(v.size() == w.size());
    double d = 0; // distance
    for(int i = 0; i < v.size(); ++i) {
        d += (v[i] - w[i]) * (v[i] - w[i]);
    }
    return sqrt(d);
}


gcms::LibrarySearch<int, float> *libraryFromCompoundFiles(std::vector<std::string> files) {
    gcms::LibrarySearch<int, float> *lib;
    for(std::vector<std::string>::iterator it = files.begin(); it < files.end(); ++it) {
        std::vector<gcms::Compound<int, float>* > tmp = gcms::Compound<int, float>::fromDisk(it->c_str());
        lib->addCompoundSpectra(tmp);
    }
    return lib;
}


std::vector<gcms::LibraryCompound<int, float>* > libraryCompoundsFromFiles(std::vector<std::string> files) {
    std::vector<gcms::LibraryCompound<int, float>* > comps; // list with compounds from all files
    for(std::vector<std::string>::iterator it = files.begin(); it < files.end(); ++it) {
        std::vector<gcms::LibraryCompound<int, float>* > tmp = gcms::LibraryCompound<int, float>::fromDisk(it->c_str());
        for(int i = 0; i < tmp.size(); ++i) {
            comps.push_back(tmp.at(i));
        }
    }
    return comps;
}


/**
 * @brief Write simple Dot-format file for graphviz visualization.
 *
 * @param mat    Network matrix. Values > 0 represent nodes
 * @param fname  Output filename
 * @param labs   Node labels
 */
void matToDot(const std::vector<std::vector<double> > &mat, std::string fname, std::vector<std::string> labs) {
    std::ofstream ofs;
    std::cout<<"writing " << fname<<std::endl;
    ofs.open(fname.c_str());
    //double distRange = abs(1 / max<double>(mat) - 1 / min<double>(mat));

    ofs << "graph G {" << std::endl;
    for(int i = 0; i < mat[0].size(); ++i) { // row
        for(int j = 0; j < mat.size(); ++j) { // col
            if(mat[i][j] > 0 && j > i) {
                // i--j [weight=..];
                double weight = 1 / mat[i][j];
                //double penwidth = log(1 / mat[i][j]) * 4 / distRange; // maxwidth 4
                //ofs << i << " -- " << j << " [weight=" << weight << ",penwidth=" << "1" <<"];" << std::endl;
                ofs << i << " -- " << j << " [penwidth=" << "1" <<",label="<<mat[i][j]<<"];" << std::endl;
            }
        }
    }

    // extended labels
    for(int i = 0; i < labs.size(); ++i) {
        ofs << i << " [";
        //ofs << "label=\"" << labs[i] << "\",style=filled,fillcolor=\"#ACD9FF\",";
        //ofs << "image=\"midplotssingle/"<<i<<".png\"";
        ofs << "label=<<TABLE border=\"0\" cellborder=\"0\"><TR><TD>"<<labs[i]<<"</TD></TR><TR><TD><IMG SRC=\"midplotssingle/"<<i<<".png\"/></TD></TR></TABLE>>,style=filled,fillcolor=\"#ACD9FF\",";
        ofs << "];" << std::endl;
    }
    ofs << "}" << std::endl;
    ofs.close();
}


/**
 * @brief Write matrix to comma separated value file.
 *
 * @param mat       Data matrix
 * @param fname     Output filename
 * @param colnames  Column labels
 * @param rownames  Row labels
 */
void matToCsv(const std::vector<std::vector<double> > &mat, std::string fname, std::vector<std::string> colnames, std::vector<std::string> rownames) {
    std::ofstream ofs;
    std::cout<<"writing " << fname;
    ofs.open(fname.c_str());
    if(colnames.size()) {
        ofs << ",";
        for(int i = 0; i < colnames.size(); ++i) ofs << "\"" << colnames[i] << "\",";
        ofs << std::endl;
    }
    for(int i = 0; i < mat[0].size(); ++i) { // row
        if(rownames.size()) ofs <<  "\"" << rownames[i] << "\",";
        for(int j = 0; j < mat.size(); ++j) { // col
            ofs << "\"" << mat[i][j] << "\",";
        }
        ofs << std::endl;
    }
    ofs.close();
}

/**
 * @brief (Dis-)similarity score for two mass-isotopomer abundances used for Needleman-Wunsch-scoring.
 *
 * @param v1 Value 1
 * @param v2 Value 2
 * @return float The score
 */
float nwScoreMID(float v1, float v2) {
    return std::fabs(v1 - v2);
}


/**
 * @brief Output matrix to ostream
 *
 * @param mat The matrix
 * @param os The ostream (defaults to std::cerr)
 */
void printMat(const std::vector<std::vector<double> > &mat, std::ostream &os) {
    os<<std::setw(6) << std::setprecision(4) << std::resetiosflags(std::ios::right);
    for(int i = 0; i < mat.size(); ++i) { // row
        for(int j = 0; j < mat[i].size(); ++j) { // col
            os<<mat[i][j]<<"\t\t";
        }
        os<<std::endl;
    }
}

void printMat(const std::vector<std::vector<char> > &mat, std::ostream &os) {
    for(int i = 0; i < mat.size(); ++i) { // row
        for(int j = 0; j < mat[i].size(); ++j) { // col
            os<<mat[i][j]<<" ";
        }
        os<<std::endl;
    }
}


void printMat(const std::vector<std::vector<bool> > &mat, std::ostream &os) {
    for(int i = 0; i < mat.size(); ++i) { // row
        for(int j = 0; j < mat[i].size(); ++j) { // col
            os<<(mat[i][j]?"1":"0")<<" ";
        }
        os<<std::endl;
    }
}


double getCosineCorrelation(const std::vector<double> &v, const std::vector<double> &w)
{
    assert(v.size() == w.size());
    double dot = 0;
    for(int i = 0; i < v.size(); ++i)
        dot += v[i] * w[i];
    return dot;
}


double getCanberraDistance(const std::vector<double> &v, const std::vector<double> &w)
{
    assert(v.size() == w.size());
    double d = 0; // distance
    for(int i = 0; i < v.size(); ++i) {
        d += fabs(v[i] - w[i]) / (fabs(v[i]) + fabs(w[i]));
    }
    return d;
}


double getManhattanDistance(const std::vector<double> &v, const std::vector<double> &w)
{
    assert(v.size() == w.size());
    double d = 0; // distance
    for(int i = 0; i < v.size(); ++i) {
        d += fabs(v[i] - w[i]);
    }
    return d;
}


double getCustomDistance(const std::vector<double> &v, const std::vector<double> &w, const std::vector<double> &vCI, const std::vector<double> &wCI)
{
    assert(v.size() == w.size());
    double d = 0; // distance
 /*
    for(int i = 0; i < v.size(); ++i) {
        d += fabs(v[i] - w[i]) * vCI[i] * wCI[i];
    }
    */
    /*for(int i = 0; i < v.size(); ++i) {
        d += sqrt(v[i] * w[i]);
    }*/

    //d = -log(d);

/*
    for(int i = 0; i < v.size(); ++i) {
        if(v[i] < 0.01 || w[i] < 0.01) // i.e. not present (-> summand would become one)
            continue;
        d += fabs(v[i] - w[i]) / (fabs(v[i]) + fabs(w[i]));
    }
    */
/*
    double m;
    for(int i = 0; i < v.size(); ++i) {
        m = std::min(v[i], w[i]);
        d += m < 0 ? 0 : m;
    }

    d = 1/d;
    if(d > std::numeric_limits<double>::max())
        d = 5;
        */

    // pearson:
    // means
    double m_v = 0, m_w = 0;
    for(int i = 0; i < v.size(); ++i) {
        m_v += v[i];
        m_w += w[i];
    }
    m_v /= v.size();
    m_w /= w.size();
    // sds
    double sd_v = 0, sd_w = 0;
    for(int i = 0; i < v.size(); ++i) {
        sd_v += (m_v - v[i]) * (m_v - v[i]);
        sd_w += (m_w - w[i]) * (m_w - w[i]);
    }
    sd_v = sqrt(sd_v / v.size());
    sd_w = sqrt(sd_w / w.size());

    for(int i = 0; i < v.size(); ++i) {
        d += (v[i] - m_v ) / sd_v * (w[i] - m_w ) / sd_w;
    }

    d = 1 - d / v.size();

    return d;
}


std::vector<double> basePeakNormalization(const std::vector<double> &v)
{
    std::vector<double> res(v.size());
    if(v.size()) {
        // find max
        double max = v[0];
        for(int i = 0; i < v.size(); ++i) {
            max = std::max(max, v[i]);
        }
        // normalize
        for(int i = 0; i < v.size(); ++i) {
            res[i] = v[i] / max;
        }
    }
    return res;
}

std::vector<double> sumNormalization(const std::vector<double> &v)
{
    std::vector<double> res(v.size());
    if(v.size()) {
        // find sum
        double sum = 0;
        for(int i = 0; i < v.size(); ++i) {
            sum += v[i];
        }
        // normalize
        for(int i = 0; i < v.size(); ++i) {
            res[i] = v[i] / sum;
        }
    }
    return res;
}
