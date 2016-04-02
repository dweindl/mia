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

/* Adapted from NTFD source code */

#ifndef COMPOUNDDETECTOR_H
#define COMPOUNDDETECTOR_H

#include <QThread>
#include <QString>

#include<vector>

#include"compound.h"
#include"gcmsdiskscan.h"
#include"abstractpeakdetector.h"

class CompoundDetector
{
      //Q_OBJECT
public:
      CompoundDetector (const QString& file, gcms::AbstractPeakDetector<int,float>*, bool redetect);
      ~CompoundDetector();

      void run();
      QString getFileName()const;
      QString getErrorMessage() const;
      std::vector<gcms::Compound<int,float>*> getCompounds();

private:
      const gcms::GCMSDiskScan<int,float>* scan;
      gcms::AbstractPeakDetector<int,float>* det;
      bool baseline;
      QString error;
      std::vector<gcms::Compound<int,float>*> comps;
    //  bool comps_returned;
      QString file;
      bool redetect;
};

#endif // COMPOUNDDETECTOR_H
