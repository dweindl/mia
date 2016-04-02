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

/* Adapted from NTFD and MetaboliteDetector source code */

#include "compounddetector.h"

#include<QFileInfo>
#include<QDir>

#include"gcmsanalyzer.h"
#include"gcmsmemoryscan.h"
#include"gcmsdiskscan.h"
#include"ricalculator.h"


CompoundDetector::CompoundDetector (const QString& file, gcms::AbstractPeakDetector<int,float>* det, bool redetect )
            :
            scan(0),
            det ( det ),
            //comps_returned(false),
            file(file),
            redetect(redetect)
{
}

CompoundDetector::~CompoundDetector ()
{

      /*if(!comps_returned){
            //delete compounds
            for(std::vector<gcms::Compound<int,float>*>::const_iterator it=comps.begin();it!=comps.end();++it)
            {
                  delete *it;
            }
            comps.clear();
      }*/
}

void CompoundDetector::run()
{
      try
      {            
            QFileInfo info(file);            
            QDir dir(info.absoluteDir());

            std::cout<<"Redetect:"<<redetect<<std::endl;
            std::cout<<"File path:"<<info.completeBaseName().toStdString()<<std::endl;
            std::cout<<"Deconvolution width:"<< gcms::GCMSSettings::AN_DECONVOLUTION_WIDTH<<std::endl;
            std::cout<<"Peak height:"<<gcms::GCMSSettings::AN_MIN_PEAK_HEIGHT<<std::endl;
            std::cout<<"Peak threshold begin:"<<gcms::GCMSSettings::AN_PEAK_THRESHOLD_BEGIN<<std::endl;
            std::cout<<"Peak threshold end:"<< gcms::GCMSSettings::AN_PEAK_THRESHOLD_END<<std::endl;

            if(!redetect && dir.exists(info.completeBaseName()+".cmp"))
                return;

            gcms::GCMSDiskScan<int,float>* scan;
            QFile f (file + ".bin"); //QFile f ( file +tr ( ".bin" ) );
            f.open ( QIODevice::ReadOnly );
            qint64 size=f.size();
            f.close();
            std::cout<<"INFORMATION: file size: "<<size<<std::endl;
            if ( size<=150000000 )
            {
                std::cout<<"INFORMATION: Using GCMSMemoryScan."<<std::endl;
                scan=new gcms::GCMSMemoryScan<int,float> ( file.toStdString().c_str());
            }
            else
            {
                scan=new gcms::GCMSDiskScan<int,float> ( file.toStdString().c_str() );
                std::cout<<"INFORMATION: Using GCMSDiskScan."<<std::endl;
            }

            scan->setBaselineCorrectionEnabled ( baseline );

            gcms::GCMSAnalyzer<int,float> an ( *scan, det );

            an.detectAllPeaks ( true, 10, 0.05 ); //TODO: SETTINGS
            an.deconvoluteChromatogram();
            comps=an.getCompounds();

            QString file_out ( file );
            file_out+=".cmp";
            gcms::Compound<int,float>::toDisk ( comps,file_out.toStdString().c_str() );

            for ( std::vector<gcms::Compound<int,float>* >::const_iterator it=comps.begin();it!=comps.end();++it )
            {
                delete *it;
            }
            comps.clear();

            if(scan)
            {
                delete scan;
            }
    }
    catch ( ... )
    {
        error= "Error detecting compounds " ;
    }
}

QString CompoundDetector::getFileName()const
{
      return file;
}

std::vector<gcms::Compound<int,float>*> CompoundDetector::getCompounds()
{
     // comps_returned=true;
      return comps;
}

QString CompoundDetector::getErrorMessage() const
{
      return error;
}

