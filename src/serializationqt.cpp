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

#include "serializationqt.h"

#include <QVector>

namespace mia {

QDataStream &operator <<(QDataStream &out, const mia::NetworkLayer *layer)
{
    out << QString("Layer");

    out << layer->settings << layer->cmpLab << layer->cmpUnlab << layer->dists << layer->distsCut;
    out << layer->mids << layer->midsAll << layer->nodeLabs << layer->nodeLabsAll;

    out << QString("LID");

    // serialize labid::Labelidentificator
    std::ostringstream ss;
    layer->lid->toStream(ss);
    std::string str = ss.str();
    uint strlen = str.length();
    out << strlen;
    out.writeRawData(str.c_str(), strlen);

    out << QString("/Layer");

    return out;
}

QDataStream &operator >>(QDataStream &in, mia::NetworkLayer *&layer) throw(DeserializationException)
{
    QString magic;
    in >> magic;
    if(magic != "Layer")
        throw(DeserializationException("NetworkLayer -- Layer"));

    mia::Settings s;
    in >> s;
    layer = new mia::NetworkLayer(s);

    in >> layer->cmpLab >> layer->cmpUnlab >> layer->dists >> layer->distsCut;
    in >> layer->mids >> layer->midsAll >> layer->nodeLabs >> layer->nodeLabsAll;

    // read labid::Labelidentificator
    in >> magic;
    if(magic != "LID")
        throw(DeserializationException("NetworkLayer -- LID"));
    uint size;
    in >> size;

    char *buf = (char*)malloc(size*sizeof(char));
    in.readRawData(buf, size);

    std::stringstream ss;
    ss.write(buf, size);
    free(buf);
    layer->lid = labid::LabelIdentificator::fromStream(ss, SERIALIZATIONQT_H_MD_VERSION);

    in >> magic;
    if(magic != "/Layer")
        throw(DeserializationException("NetworkLayer -- /Layer"));

    return in;
}

QDataStream &operator <<(QDataStream &out, const labid::LabeledCompound * lc)
{
    out<<QString("LabeledCompound");

    std::ostringstream ss;
    lc->toStream(ss);
    std::string str = ss.str();
    out << (uint)(str.length());
    out.writeRawData(str.c_str(), str.length());

    out<<QString("/LabeledCompound");

    return out;
}

QDataStream &operator >>(QDataStream &in, labid::LabeledCompound *&lc) throw(DeserializationException)
{
    QString magic;
    in >> magic;
    if(magic != "LabeledCompound")
        throw(DeserializationException("LabeledCompound -- LabeledCompound"));

    uint bufsize;
    in >> bufsize;
    char buf[bufsize];
    in.readRawData(buf, bufsize);
    std::stringstream ss;
    ss.write(buf, bufsize);
    lc = labid::LabeledCompound::fromStream(ss, SERIALIZATIONQT_H_MD_VERSION);

    in >> magic;
    if(magic != "/LabeledCompound")
        throw(DeserializationException("LabeledCompound -- LabeledCompound"));

    return in;
}


QDataStream &operator <<(QDataStream &out, const std::vector<labid::LabeledCompound *> lcs)
{
    uint size = lcs.size();
    out << size;
    for(std::vector<labid::LabeledCompound *>::const_iterator it = lcs.begin(); it != lcs.end(); ++it) {
        out << *it;
    }

    return out;
}

QDataStream &operator >>(QDataStream &in, std::vector<labid::LabeledCompound *> &lcs) throw(DeserializationException)
{
    lcs.clear();

    uint size1;
    in >> size1;

    lcs.reserve(size1);

    while(size1--) {
        labid::LabeledCompound *lc;
        in >> lc;
        lcs.push_back(lc);
    }

    return in;
}

QDataStream &operator <<(QDataStream &out, const std::vector<labid::LISpectrum *> lis)
{
    uint size = lis.size();
    out << size;

    for(std::vector<labid::LISpectrum *>::const_iterator it = lis.begin(); it != lis.end(); ++it) {
        out << *it;
    }

    return out;
}

QDataStream &operator >>(QDataStream &in, std::vector<labid::LISpectrum *> &lis) throw(DeserializationException)
{
    lis.clear();

    uint size1;
    in >> size1;
    lis.reserve(size1);

    while(size1--) {
        labid::LISpectrum *li;
        in>>li;
        lis.push_back(li);
    }

    return in;
}


QDataStream &operator <<(QDataStream &out, labid::LISpectrum const* lis)
{
    out << QString("LISpectrum");

    std::ostringstream ss;
    lis->toStream(ss);
    std::string str = ss.str();
    uint strlen = str.length();
    out << strlen;
    out.writeRawData(str.c_str(), strlen);

    return out;
}

QDataStream &operator >>(QDataStream &in, labid::LISpectrum* &lis) throw(DeserializationException)
{
    QString magic;
    in >> magic;
    if(magic != "LISpectrum")
        throw(DeserializationException("LISpectrum -- LISpectrum"));

    uint bufsize;
    in >> bufsize;

    char *buf = (char*)malloc(bufsize*sizeof(char));
    in.readRawData(buf, bufsize);
    std::stringstream ss;
    ss.write(buf, bufsize);
    free(buf);

    lis = labid::LISpectrum::fromStream(ss, SERIALIZATIONQT_H_MD_VERSION);

    return in;
}

QDataStream &operator <<(QDataStream &out, const std::vector<std::vector<double> > v)
{
    uint size1 = v.size();
    out << size1;
    for(std::vector<std::vector<double> >::const_iterator it1 = v.begin(); it1 != v.end(); ++it1){
        QVector<double> v1 = QVector<double>::fromStdVector(*it1);
        out << v1;
    }
    return out;
}

QDataStream &operator >>(QDataStream &in, std::vector<std::vector<double> > &v) throw(DeserializationException)
{
    v.clear();

    uint size1;
    in >> size1;
    v.reserve(size1);
    while(size1--) {
        QVector<double> v1;
        in >> v1;
        v.push_back(v1.toStdVector());
    }
    return in;
}

QDataStream &operator <<(QDataStream &out, const std::vector<std::string> v)
{
    uint size = v.size();
    out << size;
    for(std::vector<std::string>::const_iterator it = v.begin(); it != v.end(); ++it){
        out << QString::fromStdString(*it);
    }
    return out;
}

QDataStream &operator >>(QDataStream &in, std::vector<std::string> &v) throw(DeserializationException)
{
    v.clear();

    uint size;
    in >> size;
    while(size--) {
        QString qs;
        in>>qs;
        v.push_back(qs.toStdString());
    }
    return in;
}

QDataStream &operator <<(QDataStream &out, const gcms::Compound<int, float> *c)
{
    out<<QString("Compound");

    std::ostringstream ss;
    c->toStream(ss);

    std::string str = ss.str();
    uint strlen = str.length();
    out << strlen;
    out.writeRawData(str.c_str(), strlen);

    return out;
}

QDataStream &operator >>(QDataStream &in, gcms::Compound<int, float> *c) throw(DeserializationException)
{
    QString magic;
    in >> magic;
    if(magic != "Compound")
        throw(DeserializationException("Compound -- Compound"));

    uint bufsize;
    in >> bufsize;
    char *buf = (char*)malloc(bufsize*sizeof(char));
    in.readRawData(buf, bufsize);
    std::stringstream ss;
    ss.write(buf, bufsize);
    free(buf);

    c = gcms::Compound<int,float>::fromStream(ss, SERIALIZATIONQT_H_MD_VERSION);

    return in;
}

QDataStream &operator <<(QDataStream &out, const gcms::LibraryCompound<int, float> *c)
{
    out<<QString("LibraryCompound");

    std::ostringstream ss;
    c->toStream(ss);

    std::string str = ss.str();
    uint strlen = str.length();
    out << strlen;
    out.writeRawData(str.c_str(), strlen);

    return out;
}

QDataStream &operator >>(QDataStream &in, gcms::LibraryCompound<int, float> *c) throw(DeserializationException)
{
    QString magic;
    in >> magic;
    if(magic != "LibraryCompound")
        throw(DeserializationException("LibraryCompound -- LibraryCompound"));

    uint bufsize;
    in >> bufsize;
    char *buf = (char*)malloc(bufsize*sizeof(char));
    in.readRawData(buf, bufsize);
    std::stringstream ss;
    ss.write(buf, bufsize);
    free(buf);

    c = gcms::LibraryCompound<int,float>::fromStream(ss, SERIALIZATIONQT_H_MD_VERSION);

    return in;
}

QDataStream &operator <<(QDataStream &out, const mia::Settings s)
{
    out<<QString("Settings");

    out << QString::fromStdString(s.experiment) << s.cmp_id_use_ri << s.cmp_id_ri_tol << s.cmp_id_score_cutoff
        << s.labels_max_hits << s.gcms_pure_factor << s.gcms_impure_factor
           //<< std::set<std::pair<float,float> > cmp_id_mass_filter;
        << s.cmp_matching_ri_tol << s.cmp_matching_score_cutoff
        << QString::fromStdString(s.cmp_id_library) << s.unlabFiles << s.labFiles << s.lid_filter_by_conf_interval
        << s.lid_min_signal_to_noise << s.lid_required_spec_freq << s.lid_req_label_amount
        << s.lid_req_r2 << s.lid_min_frag_num << s.lid_sensitivity << s.lid_maximal_frag_dev
        << s.lid_correction_ratio << s.nw_gap_penalty<< s.nw_exclude_m0
        << s.mid_distance_cutoff;
    return out;
}

QDataStream &operator >>(QDataStream &in, mia::Settings &s) throw(DeserializationException)
{
    QString magic;
    in >> magic;
    if(magic != "Settings")
        throw(DeserializationException("Settings -- Settings"));

    QString t, lib;
    in >> t >> s.cmp_id_use_ri >> s.cmp_id_ri_tol >> s.cmp_id_score_cutoff;
    in >> s.labels_max_hits >> s.gcms_pure_factor >> s.gcms_impure_factor;
           //>> std::set<std::pair<float,float> > cmp_id_mass_filter;
    in >> s.cmp_matching_ri_tol >> s.cmp_matching_score_cutoff;
    in >> lib >> s.unlabFiles;
    in >> s.labFiles;
    in >> s.lid_filter_by_conf_interval >> s.lid_min_signal_to_noise >>
          s.lid_required_spec_freq >> s.lid_req_label_amount
       >> s.lid_req_r2 >> s.lid_min_frag_num >> s.lid_sensitivity >> s.lid_maximal_frag_dev
        >> s.lid_correction_ratio >> s.nw_gap_penalty>> s.nw_exclude_m0
        >> s.mid_distance_cutoff;
    s.experiment = t.toStdString();
    s.cmp_id_library = lib.toStdString();

    return in;
}
}
