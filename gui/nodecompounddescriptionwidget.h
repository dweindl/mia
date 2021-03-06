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

#ifndef NODECOMPOUNDDESCRIPTIONWIDGET_H
#define NODECOMPOUNDDESCRIPTIONWIDGET_H

#include <QtGui>
#include <QtWidgets>

#include "../src/nodecompound.h"

namespace mia {

class NodeCompoundDescriptionWidget : public QScrollArea
{
    Q_OBJECT
public:
    NodeCompoundDescriptionWidget(NodeCompound *nc, QWidget *parent = 0);

    QWidget* getMIDPlotsWidget();
signals:

public slots:
    void close();
    void exportWindow();

private:
    QWidget *w;
    NodeCompound *nc;
    QGraphicsScene *scene;
};

}
#endif // NODECOMPOUNDDESCRIPTIONWIDGET_H
