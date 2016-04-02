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

#ifndef NODEWIDGET_H
#define NODEWIDGET_H

#include <QWidget>
#include <QGraphicsObject>
#include <QEvent>
#ifdef MIA_WITH_METABOBASE
  #include <QNetworkReply>
#endif
#include "src/nodecompound.h"
#include "midplot.h"

namespace mia  {

/**
 * @brief The NodeWidget class does all the plotting for MID graphs, ... for a node.
 */
class NodeWidget : public QGraphicsObject
{
    Q_OBJECT

public:

    enum PLOT_LAYOUT {
        LAYOUT_VERTICAL,
        LAYOUT_HORIZONTAL,
        LAYOUT_SQUARE
    }; /**< Show plots in row, column or multirow*/

    explicit NodeWidget();
    explicit NodeWidget(NodeCompound *nc, QList<QColor> experimentColors, bool multiExperiment = false);

    ~NodeWidget();

    //bool event(QEvent *event);

    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
#ifdef MIA_WITH_METABOBASE
    void addStructure(QString url);

public slots:
    void structureReply(QNetworkReply*reply);
#endif

private:
    QRectF boundRect;
    QString genToolTipHeader();
    QString genToolTipText(std::string);
    std::vector<MIDPlot*> midPlots;
    NodeCompound *nc; /**< The label information for this node's compound. */
    QList<QColor> experimentColors; /**< Colors for the different experiments/conditions. */
    PLOT_LAYOUT plotLayout;

#ifdef MIA_WITH_METABOBASE
    QPixmap structureImage;
#endif

};
}
#endif // NODEWIDGET_H
