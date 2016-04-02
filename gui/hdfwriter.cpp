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

#include "hdfwriter.h"
#include <algorithm>


void HDFStringWriter::operator ()(const std::vector<std::string>::value_type &v)
{
    // Select the file position, 1 record at position 'pos'
    hsize_t count[] = { 1 } ;
    hsize_t offset[] = { m_pos++ } ;
    m_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    const char * s = v.c_str ();
    m_dataset.write(&s, m_datatype
                    , m_memspace
                    , m_dataspace
                    , H5P_DEFAULT
                             );
}

void HDFStringWriter::writeVector(H5::Group group, std::string ds, const std::vector<std::string> &v)
{
    hsize_t     dims[] = { v.size ()  } ;
    H5::DataSpace dataspace(1, dims);
//   hsize_t s = distMats.size();
  //  H5::DataSpace dataspace(1, &s);
   // H5::DataSet dsStr = h5f.createDataSet("/Dists/Names", H5::PredType::ALPHA_U8, dataspace);

    dims[0] = 1;
    H5::DataSpace memspace(1, dims);
    H5::DataType datatype = H5::DataType(H5T_STRING, H5T_VARIABLE);

//    hid_t datatype = H5Tcopy (H5T_C_S1);
//    H5Tset_size (datatype, H5T_VARIABLE);

    H5::DataSet dataset = group.createDataSet(ds, datatype, dataspace);
    //hid_t dataset = H5Dcreate1 (group, "files", datatype
      //                         , dataspace, H5P_DEFAULT);

    //
    // Select the "memory" to be written out - just 1 record.
    hsize_t offset[] = { 0 } ;
    hsize_t count[] = { 1 } ;
     memspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    std::for_each (v.begin ()
        , v.end ()
        , HDFStringWriter (dataset, datatype, dataspace, memspace));

    dataset.close();
    dataspace.close();
    memspace.close();
    datatype.close();
}


