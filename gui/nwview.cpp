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

#include<QtWidgets>
#include<QSvgGenerator>
#include<iostream>
#include "nwview.h"
#include "src/nodecompound.h"
#include "midplot.h"
#include "graphvizmia.h"

namespace mia {

#ifndef QT_NO_WHEELEVENT
void NWGView::wheelEvent(QWheelEvent *e)
{
    if (e->modifiers() & Qt::ControlModifier) {
        if (e->delta() > 0)
            view->zoomIn(zoomStep);
        else
            view->zoomOut(zoomStep);
        e->accept();
    } else {
        QGraphicsView::wheelEvent(e);
    }
}
#endif

void NWGView::mouseDoubleClickEvent(QMouseEvent *e)
{
    QGraphicsItem *gi = itemAt(e->pos());

    if(gi) {
        QVariant v = gi->data(g->getDataIndex("mynode")); // use for tooltip and so on
        if(!v.isNull()) {
            if(v.canConvert<void *>()) {
                mynode *mn = reinterpret_cast<mynode *>(v.value<void *>());
                NodeCompound *nc = mn->nc;
// TODO: Neighborplotwindow


                std::cerr<<nc->toString()<<std::endl<<std::endl;
            } else {
                std::cerr<<"no convert "<< v.typeName() <<std::endl;
            }
        }
    }
}

void NWGView::mouseMoveEvent(QMouseEvent *e)
{
    if(e->buttons() & Qt::MiddleButton) {
        QPoint center = rect().center();
        const QPoint newCenter = center - e->globalPos() + lastMousePos;
        //std::cerr<<"Center "<<center.x()<< " " <<center.y()<<std::endl;
        //std::cerr<<"NCenter "<<newCenter.x()<< " " <<newCenter.y()<<std::endl;
        centerOn(mapToScene(newCenter));
        e->accept();
    }
    e->ignore();
    lastMousePos = e->globalPos();
}

void NWGView::mouseReleaseEvent(QMouseEvent *e)
{
    static NodeCompound *lastSelectedNodeCompound;

    //foreach (QGraphicsItem *gi, scene()->items(QRectF(loc.x(), loc.y(), 1, 1))) {
        QGraphicsItem *gi = itemAt(e->pos());
        if(gi) {
            //TODO: check for midplot to hide onclick  std::cerr<<"midplot"<<std::endl;
            QVariant v = gi->data(g->getDataIndex("mynode")); // use for tooltip and so on
            if(!v.isNull()) {
                if(v.canConvert<void *>()) {
                    mynode *mn = reinterpret_cast<mynode *>(v.value<void *>());
                    NodeCompound *nc = mn->nc;
                    std::cerr<<nc->toString()<<std::endl<<std::endl;

                    NodeCompoundDescriptionWidget *desc = new NodeCompoundDescriptionWidget(nc);
                    //desc->move(loc.x(), loc.y());
                    desc->show();
                    //QGraphicsProxyWidget *p = scene()->addWidget(desc);


                    if(lastSelectedNodeCompound) {
                        // show distance between current and previously selected component
                        std::vector<std::string> ex = nc->getExperiments();
                        for(int i = 0; i < ex.size(); ++i) {


                        }
                    }
                    lastSelectedNodeCompound = nc;
                } else {
                    std::cerr<<"no convert "<< v.typeName() <<std::endl;
                }
            }
       // }
        }
}

void NWGView::keyPressEvent(QKeyEvent *e)
{
    if(e->key() == Qt::Key_Plus) {
        view->zoomIn(zoomStep);
    } else if(e->key() == Qt::Key_Minus) {
        view->zoomOut(zoomStep);
    } else {
        QGraphicsView::keyPressEvent(e);
    }
}

NWView::NWView(QWidget *parent) :
    QFrame(parent)
{
    gv = new NWGView(this, this);
    init();
}

NWView::NWView(QGraphicsScene *gs, graphvizqt::Graph *g, QWidget *parent)
{
    gv = new NWGView(this, g, this);
    gv->setScene(gs);
    init();
}


QGraphicsView *NWView::view() const
{
    return static_cast<QGraphicsView *>(gv);
}

void NWView::zoomIn(int level)
{
    zoom = level;
    setupMatrix();
}

void NWView::zoomOut(int level)
{
    zoom = 1.0/level;
    setupMatrix();
}

void NWView::exportImage()
{
    // TODO remember folder in settings
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
    QPainter painter(&svgGen);
    painter.begin(&svgGen);
    gv->scene()->render(&painter);
    painter.end();

}

void NWView::resetView()
{
    // original scene region, original zoom level
    zoom = 1.0;
    QRectF itemsBB = gv->scene()->itemsBoundingRect();
    gv->centerOn(itemsBB.center());
    gv->setMatrix(originalMatrix);
    gv->fitInView(itemsBB, Qt::KeepAspectRatio);
}


void NWView::setupMatrix()
{
    // std::cerr<<"Zoom: "<<zoom<<std::endl;
    gv->scale(zoom, zoom);
}

void NWView::init()
{
    gv->setDragMode(QGraphicsView::RubberBandDrag);
    gv->setRenderHint(QPainter::Antialiasing);

    zoom = 1;

    QToolBar *toolBar = new QToolBar(this);

    QAction *actZoomOriginal = new QAction(QIcon(":/gui/icons/zoom-original.png"), tr("Reset Zoom"), this);
    toolBar->addAction(actZoomOriginal);
    QAction *actRefresh = new QAction(QIcon(":/gui/icons/view-refresh.png"), tr("Refresh view"), this);
    toolBar->addAction(actRefresh);
    QAction *actExportImage = new QAction(QIcon(":/gui/icons/video-x-mng.png"), tr("Export Image"), this);
    toolBar->addAction(actExportImage);

    QGridLayout *layout = new QGridLayout(this);
    layout->addWidget(toolBar);
    layout->addWidget(gv);
    setLayout(layout);

    connect(actExportImage, SIGNAL(triggered()), this, SLOT(exportImage()));
    connect(actRefresh, SIGNAL(triggered()), gv, SLOT(invalidateScene()));
    connect(actZoomOriginal, SIGNAL(triggered()), this, SLOT(resetView()));
    originalMatrix = gv->matrix();
    gv->mapToScene(rect().center());
}


NWGView::NWGView(NWView *v, QWidget *parent) : QGraphicsView(parent), view(v)
{
    init();
}

NWGView::NWGView(NWView *v, graphvizqt::Graph *graph, QWidget *parent) : QGraphicsView(parent), view(v), g(graph)
{
    init();
}

void NWGView::init()
{
    zoomStep = 2;
    setFocusPolicy(Qt::ClickFocus);
}

NWGView::~NWGView()
{
}


}
