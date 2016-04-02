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

#include <QHBoxLayout>
#include <QToolTip>
#include <QHelpEvent>
#include <QSettings>
#include <QImageReader>
#include <QNetworkAccessManager>
#include <QPainter>

//#include "qwt_plot_multi_barchart.h"
//#include "qwt_samples.h"
#include "nodewidget.h"
//#include "mitrendplot.h"

namespace mia {

/**
 * @brief Default constructor.
 * @param parent The parent widget.
 */
NodeWidget::NodeWidget() : QGraphicsObject()
{

}

NodeWidget::NodeWidget(NodeCompound *nc, QList<QColor> experimentColors, bool multiExperiment) : QGraphicsObject(), nc(nc), experimentColors(experimentColors)
{
    // TODO use multiExperiment to draw compound only once here, and remove from midplots
    // MID plot TODO: plot all :: 3D barplot?
    // TODO: multirows rowNum = sqrt(exps.size())

    double width = 0, height = 0;

    if(nc->getFeature("PRECURSOR_KEGG_ID").size()) {
        // guess image size and hope we'll find one later
        width = 100;
        height = 250;
    }

    std::vector<std::string> exps = nc->getExperiments();

    bool plotUnlabeledCompounds = true;

    for(int i = 0; i < exps.size(); ++i) {
        labid::LabeledCompound *lc = nc->getLabeledCompound(exps[i]);
        if(lc == (labid::LabeledCompound*)-1) {
            if(plotUnlabeledCompounds) {
                lc = nc->getCompound(exps[i]);
            } else {
                continue;
            }
        }
        QVector<double> mids = QVector<double>::fromStdVector(nc->getSelectedMID(exps[i]));

        float ion = nc->getSelectedIon(exps[i]);
        MIDPlot *mp = new MIDPlot(mids, QString::fromStdString(lc->getName()), nc->getSelectedR2(exps[i]), ion, 0);
        mp->setCI(QVector<double>::fromStdVector(nc->getSelectedCI(exps[i])));
        midPlots.push_back(mp);
        mp->setBarColor(experimentColors[experimentColors.size() >  1 ? (i % experimentColors.size()) : 0]);
        setToolTip(toolTip() + genToolTipText(exps[i]));

        // color full backbone plots
        if(nc->getMMinusNIon(exps[i], 15) == ion) {
            //mp->setBgColor(Qt::gray);
            mp->setBorderColor(Qt::magenta);
        } else if(nc->getMMinusNIon(exps[i], 57) == ion) {
            //mp->setBgColor(Qt::gray);
            mp->setBorderColor(Qt::magenta);
        }

        // adjust color to labeling
        //mp->setBgColor(QColor(Qt::lightGray).lighter().darker(100 + 500 * (1 - mids[0])));
        mp->setBgColor(QColor(Qt::lightGray).lighter().darker(100 + 500 * (mids[0] - 0.5)));
    }


    plotLayout = (PLOT_LAYOUT) QSettings().value("nodewidget_layout_direction", LAYOUT_VERTICAL).toInt();

    switch (plotLayout) {

    case LAYOUT_HORIZONTAL:
        for(int i = 0; i < midPlots.size(); ++i) {
            MIDPlot *mp = midPlots[i];
            width += mp->sizeHint().width();
            height = qMax(height, (double)mp->sizeHint().height());
        }
        break;

    case LAYOUT_VERTICAL:
        for(int i = 0; i < exps.size(); ++i) {
            MIDPlot *mp = midPlots[i];
            height += mp->sizeHint().height();
            width = qMax(width, (double)mp->sizeHint().width());
        }

        break;
    case LAYOUT_SQUARE:
        1;
    }


    boundRect = QRectF(0, 0, width, height);

    //TODO barplot? qwt barplot: // do only if all ions are the same
    //QwtPlotMultiBarChart *bc = new QwtPlotMultiBarChart("Testplot");
    //QwtSetSample ss = QwtSetSample();
    //bc->setSamples();
    //layout->addWidget(mp);

    // Trendline for isotopomers
    if(exps.size() > 1) {
        //layout->addWidget(new MITrendplot(nc, this));
    }

#ifdef MIA_WITH_METABOBASE
    // TODO display all ions midplot tooltip
    if(QSettings().value("nodewidget_show_structure", LAYOUT_VERTICAL).toBool())
        addStructure(QString::fromStdString(nc->getFeature("PRECURSOR_KEGG_ID"))); // TODO or use METABOBASE_COMPOUND and find ALL precursors
#endif
}

NodeWidget::~NodeWidget()
{
    for(int i = 0; i < midPlots.size(); ++i) {
        delete midPlots[i];
    }
}

//bool NodeWidget::event(QEvent *event)
//{
/*    if (event->type() == QEvent::ToolTip) {
        QHelpEvent *helpEvent = static_cast<QHelpEvent *>(event);
        int index = itemAt(helpEvent->pos());
        if (index != -1) {
            QToolTip::showText(helpEvent->globalPos(), );
        } else {
            QToolTip::hideText();
            event->ignore();
        }

        return true;
    }*/
    //return QGraphicsItem::event(event);
//}

QRectF NodeWidget::boundingRect() const
{
    return boundRect;
}

void NodeWidget::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    // plot drawing position
    QPoint pos(0, 0);

#ifdef MIA_WITH_METABOBASE
    // the structure image
    if(structureImage.width()) {
        pos.setX((boundRect.width() - structureImage.width()) / 2);
        painter->drawPixmap(pos, structureImage);
        pos.setX(0);
        pos.setY(structureImage.height());
    }
#endif

    // the MID plots
    for(int i = 0; i < midPlots.size(); ++i) {
        MIDPlot *mp = midPlots[i];
        mp->render(painter, pos);

        QSize s = mp->sizeHint();

        switch (plotLayout) {
        case LAYOUT_HORIZONTAL:
            pos.setX(pos.x() + s.width()); // move pos right
            break;
        case LAYOUT_VERTICAL:
            pos.setY(pos.y() + s.height()); // move pos down
            break;
        case LAYOUT_SQUARE:
            // TODO
            break;
        }
    }
}

#ifdef MIA_WITH_METABOBASE
/**
 * @brief Add chemical structure to plot
 * @param id KEGG or metabobase ID
 */
void NodeWidget::addStructure(QString id)
{
    if(!id.length()) return;

    QNetworkAccessManager *nam;
    nam = new QNetworkAccessManager(this);

    connect(nam, SIGNAL(finished(QNetworkReply*)), this, SLOT(structureReply(QNetworkReply*)));

    QUrl url(QString("XXXXX/mddb/struc_svg.php?cid=").append(id));

    QNetworkReply *reply = nam->get(QNetworkRequest(url));
}

/**
 * @brief Reply for the image http get request.
 * @param reply
 */
void NodeWidget::structureReply(QNetworkReply *reply)
{
    if(reply->error() == QNetworkReply::NoError) {
        QImageReader imageReader(reply);
        structureImage = QPixmap::fromImage(imageReader.read());
    }
}

#endif

/**
 * @brief Generate summary text for the compound in the given experiment to be displayed as tooltip.
 * @param exp Which experiment.
 * @return Summary text.
 */
QString NodeWidget::genToolTipText(std::string exp)
{
    QString tt = QString("<center><b>%1 - %2</b><br>").arg(QString::fromStdString(nc->getCompoundName()), QString::fromStdString(exp));
    labid::LabeledCompound* lc = nc->getLabeledCompound(exp);
    if(lc == (labid::LabeledCompound*)-1) {
        return tt.append("<br>Unlabeled<br>");
    }

    // Retention indices of all samples
    tt.append("RIs: ");
    std::list<double> ris = lc->getRetentionIndices();
    for(std::list<double>::iterator it = ris.begin(); it != ris.end(); ++it) {
        if(it != ris.begin()) tt.append(", ");
        tt.append(QString::number(*it, 'f', 2));
    }

    // MIDs
    int selIdx = nc->getSelectedIndex(exp);
    tt.append("<br><table border=1><tr>");
    for(int i = 0; i < lc->getLabeledIons().size(); ++i) { // all ions

        if(i == selIdx)
            tt.append(
                        QString("<td style=\"color:red\"><u><b>m/z %1</b></u><br><i>R</i><sup>2</sup> = %2<table>").arg(
                            QString::number(lc->getLabeledIons()[i]),
                            QString::number(lc->getR2s()[i], 'g', 3)
                            )
                        );
        else
            tt.append(
                        QString("<td><u><b>m/z %1</b></u><br><i>R</i><sup>2</sup> = %2<table>").arg(
                            QString::number(lc->getLabeledIons()[i]),
                            QString::number(lc->getR2s()[i], 'g', 3)
                            )
                        );

        std::vector<double> mids = lc->getIsotopomers()[i];
        for(int m = 0; m < mids.size(); ++m) { // all isotopomers
            // std::cout<<nc->getCompoundName()<<" "<<nc->getANOVAPvalueForMassIsotopomer(m)<<" "<<QString::number(nc->getANOVAPvalueForMassIsotopomer(m), 'g', 3).toStdString()<<std::endl;
            tt.append(
                        QString("<tr><td>M<sub>%1</sub> = </td><td>%2</td></tr>").arg(
                            QString::number(m),
                            QString::number(mids[m], 'g', 3)
                            // QString::number(nc->getANOVAPvalueForMassIsotopomer(m), 'g', 3) //  p=%3
                            )
                        );
        }
        tt.append("</table></td>");
    }
    tt.append("</tr></table></center>");

    return tt;
}
}
