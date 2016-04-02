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

#ifndef NWVIEW_H
#define NWVIEW_H

#include <QFrame>
#include <QGraphicsView>
#include <QWheelEvent>

#include "graphvizqt.h"
#include "nodewidget.h"
#include "nodecompounddescriptionwidget.h"

namespace mia {

class NWView;

/**
 * @brief The NWGView class shows a graphvizqt::Graph and handles its events.
 */

class NWGView : public QGraphicsView
{
    Q_OBJECT
public:
    NWGView(NWView *v, QWidget *parent);
    NWGView(NWView *v, graphvizqt::Graph *graph, QWidget *parent);

    void init();

    ~NWGView();

    //void centerOn(NodeWidget *w);

protected:
#ifndef QT_NO_WHEELEVENT
    void wheelEvent(QWheelEvent *);
#endif
    void mouseDoubleClickEvent(QMouseEvent *e);
    void mouseMoveEvent(QMouseEvent *);
    void mouseReleaseEvent(QMouseEvent *);
    virtual void keyPressEvent(QKeyEvent *e);
    //void mousePressEvent(QMouseEvent *);


private:
    NWView *view;           /**< The parent widget. */
    graphvizqt::Graph *g;   /**< The graph to be drawn in the scene. */
    QPoint lastMousePos;    /**< Mouse position for scrolling. */

    double zoomStep;
};

/**
 * @brief The NWView class holds a NWGView and handles the zooming.
 */
class NWView : public QFrame
{
    Q_OBJECT
public:
    explicit NWView(QWidget *parent);
    explicit NWView(QGraphicsScene *gs, graphvizqt::Graph *g, QWidget *parent);

    QGraphicsView *view() const;

signals:
    
public slots:
    void zoomIn(int level = 1);
    void zoomOut(int level = 1);
    void exportImage();
    void resetView();

private slots:
    void setupMatrix();

private:
    void init();
    NWGView *gv;
    double zoom;
    QMatrix originalMatrix;
    QPointF originalCenter;
};

}
#endif // NWVIEW_H
