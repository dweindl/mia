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

#ifndef SERIALIZATIONQT_H
#define SERIALIZATIONQT_H

#include <vector>
#include <sstream>

#include <QDataStream>

#include "labeledcompound.h"
#include "compound.h"

#include "settings.h"
#include "miaexception.h"
#include "labelingdataset.h"
#include "networklayer.h"

namespace mia {

#define SERIALIZATIONQT_H_MD_VERSION 4

/** Serialization functions to use with QDataStreams */

QDataStream &operator << (QDataStream &out, const std::vector<labid::LabeledCompound*>);
QDataStream &operator >> (QDataStream &in, std::vector<labid::LabeledCompound*>&) throw(DeserializationException);

QDataStream &operator << (QDataStream &out, const labid::LabeledCompound*);
QDataStream &operator >> (QDataStream &in, labid::LabeledCompound*&) throw(DeserializationException);

QDataStream &operator << (QDataStream &out, const std::vector<labid::LISpectrum*>);
QDataStream &operator >> (QDataStream &in, std::vector<labid::LISpectrum*>&) throw(DeserializationException);

QDataStream &operator << (QDataStream &out, const labid::LISpectrum*);
QDataStream &operator >> (QDataStream &in, labid::LISpectrum*&) throw(DeserializationException);

QDataStream &operator << (QDataStream &out, const std::vector<std::vector<double> >);
QDataStream &operator >> (QDataStream &in, std::vector<std::vector<double> >&) throw(DeserializationException);

QDataStream &operator << (QDataStream &out, const std::vector<std::string>);
QDataStream &operator >> (QDataStream &in, std::vector<std::string>&) throw(DeserializationException);

QDataStream &operator << (QDataStream &out, const gcms::Compound<int, float>*);
QDataStream &operator >> (QDataStream &in, gcms::Compound<int, float>*) throw(DeserializationException);

QDataStream &operator << (QDataStream &out, const gcms::LibraryCompound<int, float>*);
QDataStream &operator >> (QDataStream &in, gcms::LibraryCompound<int, float>*) throw(DeserializationException);

QDataStream &operator << (QDataStream &out, const mia::Settings);
QDataStream &operator >> (QDataStream &in, mia::Settings&) throw(DeserializationException);

}
#endif // SERIALIZATIONQT_H
