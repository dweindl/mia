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

#include "mitrendplot.h"
#include<qwt_plot_curve.h>
#include<qwt_series_data.h>
namespace nwrecon {

/*
MITrendplot::MITrendplot() : QwtPlot()
{

}

MITrendplot::MITrendplot(NodeCompound *nc, QWidget* parent = 0) : QwtPlot(QwtText("Test"), parent)
{
    // TODO: only if same ion present
    QVector<std::vector<double> > mids;
    std::vector<std::string> exps = nc->getExperiments();

    int maxMI = 0;
    for(int i = 0; i < exps.size(); ++i) {
        std::vector<double> mid = nc->getSelectedMID(exps[i]);
        mids.push_back(mid);
        maxMI = std::max<int>(maxMI, mid.size());
    }

    QVector<QwtPlotCurve*> curves;
    for(int i = 0; i < maxMI; ++i) {
        curves[i] == new QwtPlotCurve(QString("M%1").arg(i));
        QwtPointSeriesData *sd = new QwtPointSeriesData();
        QVector<QPointF>* points = new QVector<QPointF>;

        int exp = 0;
        foreach (std::vector<double> mid, mids) {
            ++exp;
            if(i < mid.size()) {
                points->push_back(QPointF(exp, mid[i]));
            }
        }
        sd->setSamples(*points);
        curves[i]->setData(sd);
        curves[i]->attach(this);
    }
    replot();
}*/
}
