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

#ifndef HDFWRITER_H
#define HDFWRITER_H

#include <H5Cpp.h>
#include <string>
#include <vector>

class HDFStringWriter
{
public:
    HDFStringWriter(H5::DataSet dataset, H5::DataType datatype
                    , H5::DataSpace dataspace, H5::DataSpace memspace)
                  : m_dataset (dataset), m_datatype (datatype)
                  , m_dataspace (dataspace), m_memspace (memspace)
                  , m_pos () {}

     void operator ()(std::vector<std::string>::value_type const & v);
     static void writeVector (H5::Group group, std::string ds, std::vector<std::string> const & v);
private:
  H5::DataSet m_dataset;
  H5::DataType m_datatype;
  H5::DataSpace m_dataspace;
  H5::DataSpace m_memspace;
  int m_pos;

};

#endif // HDFWRITER_H
