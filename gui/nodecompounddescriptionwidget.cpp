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

#include "nodecompounddescriptionwidget.h"
#include <stdexcept>
#include "midplot.h"
#include <QPrinter>
#include <QtSvg>
//#include "libraryspectrumplot.h"

namespace mia {

NodeCompoundDescriptionWidget::NodeCompoundDescriptionWidget(NodeCompound *nc, QWidget *parent) : QScrollArea(parent), nc(nc)
{
    setWindowTitle(QString::fromStdString(nc->getCompoundName()));

    w = new QWidget(this);

    QGridLayout *layout = new QGridLayout(w);

    QToolBar *toolBar = new QToolBar(w);
    QAction *actExportImage = new QAction(QIcon(":/gui/icons/video-x-mng.png"), tr("Export Image"), this);
    toolBar->addAction(actExportImage);
    layout->addWidget(toolBar);
    connect(actExportImage, SIGNAL(triggered()), this, SLOT(exportWindow()));

    QString text = QString("<h1>%1</h1>").arg(QString::fromStdString(nc->getCompoundName()));

    std::vector<std::string> exps = nc->getExperiments();
    text.append(QString("<table style='padding:2px'><tr>"
                        "<td><b>Dataset</b></td>"
                        "<td><b>RI</b></td>"
                        "<td><b>RT</b></td>"
                        "<td><b>Intensity</b></td>"
                        "</tr>"));

    for(int i = 0; i < exps.size(); ++i) {
        labid::LabeledCompound *lc = nc->getCompound(exps[i]);

        text.append(QString("<tr><td>%1</td><td align=\"right\">%2</td><td align=\"right\">%3</td><td align=\"right\">%4</td></tr>").arg(QString::fromStdString(exps[i]),
                                                                                                                  QString::number(lc->getRetentionIndex(), 'f', 2),
                                                                                                                  QString::number(lc->getRetentionTime() / 1000 / 60, 'f', 2),
                                                                                                                  QString::number(lc->getTotalSignal(), 'f', 0) // TODO check if that is normalized from labid
                                                                                                                  ));
        // Intensities of unlabeled compounds
        std::list<const gcms::Compound<int, float>*> cmps = lc->getSourceSpectra();
        for(std::list<const gcms::Compound<int, float>* >::iterator it = cmps.begin(); it != cmps.end(); ++it) {
            const gcms::LibraryCompound<int, float>* cmp = dynamic_cast<const gcms::LibraryCompound<int, float>*>(*it);
            if(!cmp) continue;
            text.append(QString("<tr><td>%1</td><td align=\"right\">%2</td><td align=\"right\">%3</td><td align=\"right\">%4</td></tr>").arg(QString::fromStdString(""),
                                                                                                                      QString::number(cmp->getRetentionIndex(), 'f', 2),
                                                                                                                      QString::number(cmp->getRetentionTime() / 1000 / 60, 'f', 2),
                                                                                                                      QString::number(cmp->getTotalSignal(), 'f', 0) // TODO check if that is normalized from labid
                                                                                                                      ));
        }
    }
    text.append("</table><br>");

    /*
            // MIDs
            int selIdx = nc->getSelectedIndex(exp);
            tt.append("<br><table border=1><tr>");
            for(int i = 0; i < lc->getLabeledIons().size(); ++i) { // all ions

                if(i == selIdx)
                    tt.append(QString("<td style=\"color:red\"><u><b>%1</b> (%2)</u><table>").arg(QString::number(lc->getLabeledIons()[i]), QString::number(lc->getR2s()[i], 'g', 3)));
                else
                    tt.append(QString("<td><u><b>m/z %1</b></u> (R^2%2)<table>").arg(QString::number(lc->getLabeledIons()[i]), QString::number(lc->getR2s()[i], 'g', 3)));

                std::vector<double> mids = lc->getIsotopomers()[i];
                for(int m = 0; m < mids.size(); ++m) { // all isotopomers
                    tt.append(QString("<tr><td>M<sub>%1</sub></td><td>%2</td></tr>").arg(QString::number(m), QString::number(mids[m], 'g', 3)));
                }
                tt.append("</table></td>");
            }
            tt.append("</tr></table></center>");
*/

    QLabel *lab = new QLabel(text);
    lab->setTextFormat(Qt::RichText);
    layout->addWidget(lab);

    layout->addWidget(getMIDPlotsWidget());

/*    metabolitedetector::LibrarySpectrumPlot *specPlot = new metabolitedetector::LibrarySpectrumPlot("", this);
 #   specPlot->setPlotPeakCaptionEnabled(true);
    specPlot->setCentroidButtonEnabled(false);
    specPlot->setPlotSliderEnabled(false);
    specPlot->setZoomEnabled(true);
    layout->addWidget(specPlot);*/
    // connect(lab, SIGNAL(clicked()), this, SLOT(close()));
    w->setLayout(layout);
    this->setWidget(w);
}

QWidget *NodeCompoundDescriptionWidget::getMIDPlotsWidget()
{
    QWidget *w = new QWidget(this);
    QGridLayout *layout = new QGridLayout(w);

    std::vector<std::string> exps = nc->getExperiments();
    std::set<int> lions = nc->getAllLabeledIons();

    int row = 0;
    // plots
    for(std::set<int>::iterator it = lions.begin(); it != lions.end(); ++it) {
        int ion = *it;


        // header
        for(int i = 0; i < exps.size(); ++i) {
            layout->addWidget(new QLabel(QString::fromStdString(exps[i]), w), row, i + 1);
        }
        ++row;

        int col = 0;
        layout->addWidget(new QLabel(QString::number(ion, 'f', 0), w), row, col++);

        for(int i = 0; i < exps.size(); ++i) {
            labid::LabeledCompound *lc = nc->getCompound(exps[i]);
            const std::vector< float > ions = lc->getLabeledIons();

            std::vector<float>::const_iterator findIt = std::find(ions.begin(), ions.end(), ion);
            if(findIt != ions.end()) {
                int idx = findIt - ions.begin();
                QVector<double> isotopomers = QVector<double>::fromStdVector(lc->getIsotopomers()[idx]);
                MIDPlot* mp = new MIDPlot(isotopomers, "", lc->getR2s()[idx], ion, this);
                mp->setCI(QVector<double>::fromStdVector(lc->getConfidenceIntervals()[idx]));

                // color full backbone plots
                if(nc->getMMinusNIon(exps[i], 15) == ion) {
                    mp->setBarColor(Qt::yellow);
                } else if(nc->getMMinusNIon(exps[i], 57) == ion) {
                    mp->setBarColor(Qt::red);
                }

                layout->addWidget(mp, row, col);
            }
            ++col;
        }
        ++row;
    }

    return w;
}

void NodeCompoundDescriptionWidget::close()
{
//    destroy();
//    delete this;
}

void NodeCompoundDescriptionWidget::exportWindow()
{
    QString filename = QFileDialog::getSaveFileName(this, "Export svg", "", "SVG files (*.svg);;All files (*)");
    if(filename.isNull()) return;

    // A4 export
    double widthInch = 8.27;
    double heightInch = 11.69;

    QSvgGenerator svgGen;
    svgGen.setFileName(filename);
    svgGen.setResolution(QDesktopWidget().physicalDpiX());
    svgGen.setSize(QSize(svgGen.resolution() * widthInch, svgGen.resolution() * heightInch));
    svgGen.setViewBox(QRect(0, 0, svgGen.width(), svgGen.height()));
    svgGen.setTitle(tr("MIA export"));
    svgGen.setDescription(tr("TODO: Put settings here..."));
    w->render(&svgGen, QPoint(), QRegion(), QWidget::DrawChildren);
}

}
