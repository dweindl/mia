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

#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <math.h>

#include <QGraphicsTextItem>
#include <QDesktopWidget>
#include <QString>
#include <QtWidgets>

#include "graphvizqt.h"
#include "src/nodecompound.h"

// not needed for linux. there config6 file enough?
// adapted from graphviz dot: graphviz-2.38.0/cmd/dot/dot_builtins.c
extern gvplugin_library_t gvplugin_dot_layout_LTX_library;
extern gvplugin_library_t gvplugin_neato_layout_LTX_library;
extern gvplugin_library_t gvplugin_core_LTX_library;

lt_symlist_t lt_preloaded_symbols[] = {
        { "gvplugin_dot_layout_LTX_library", (void*)(&gvplugin_dot_layout_LTX_library) },
        { "gvplugin_neato_layout_LTX_library", (void*)(&gvplugin_neato_layout_LTX_library) },
        { "gvplugin_core_LTX_library", (void*)(&gvplugin_core_LTX_library) },
        { 0, 0 }
};
// end graphviz


namespace graphvizqt {

GraphVizQT::GraphVizQT()
{
}

const float Graph::EDGE_STYLE_BOLD_SCALING = 3;

char Graph::EDGE_ATTR_LABEL[] = "label";
char Graph::EDGE_ATTR_COLOR[] = "color";
char Graph::EDGE_ATTR_STYLE[] = "style";
char Graph::EDGE_ATTR_LABELFONTSIZE[] = "labelfontsize";
char Graph::EDGE_ATTR_LABELFONTCOLOR[] = "labelfontcolor";
char Graph::EDGE_ATTR_TOOLTIP[] = "edgetooltip";
char Graph::EDGE_ATTR_PENWIDTH[] = "penwidth";

char Graph::NODE_ATTR_COLOR[] = "color";
char Graph::NODE_ATTR_FILLCOLOR[] = "fillcolor";
char Graph::NODE_ATTR_FONTCOLOR[] = "fontcolor";
char Graph::NODE_ATTR_WIDTH[] = "width";
char Graph::NODE_ATTR_HEIGHT[] = "height";

/**
 * @brief Init empty graph.
 */
Graph::Graph()
{
    g = 0;
    //gvc = gvContext(); // for  linux
    gvc = gvContextPlugins(lt_preloaded_symbols, 0); // "0" !!! for static dot on windows
    init();

}

/**
 * @brief Init graph from graph-viz dot-string.
 * @param dot Graph definition in dot-format.
 */
Graph::Graph(std::string dot)
{
    g = 0;
    gvc = gvContext();
    char c[dot.size() + 1];
    strcpy(c, dot.c_str());
    if(!(g = agmemread(c))){
        std::cerr<<"Could not read graph."<<std::endl;
    }
    init();
}

/**
 * @brief Init graph from graph-viz dot-string.
 * @param dot Graph definition in dot-format.
 */
Graph::Graph(QString dot)
{
    g = 0;
    gvc = gvContext();
    char c[dot.size() + 1];
    strcpy(c, dot.toStdString().c_str());
    if(!(g = agmemread(c))){
        std::cerr<<"Could not read graph."<<std::endl;
    }
    init();
}

Graph::~Graph()
{
}

/**
 * @brief See addNode.
 */
void *Graph::addNode(QString name, const char *recordType, size_t recordTypeSize, QStyle *style)
{
    return addNode(name.toStdString(), recordType, recordTypeSize, style);
}

/**
 * @brief Add a new node to graph.
 * @param name Node identifier.
 * @param recordType Optional: Together with recordTypeSize,
 * allows to bind user-type record to the node. Can be retrieved later. See GraphViz agbindrec for more info. Will also be added to node GraphicsItem as QVariant if graph is drawn
 * on QGraphicsScene. @sa getDataIndex
 * @param recordTypeSize Optional: Size of the record type.
 * @param style
 * @return If recordTypeSize was provided, pointer to recordType. Otherwise 0.
 */
void *Graph::addNode(std::string name, const char *recordType, size_t recordTypeSize, QStyle *style)
{
    // std::cerr<<"Adding node <"<<name<<">"<<std::endl;

    Agnode_t *n = agnode(g, const_cast<char *>(name.c_str()), TRUE);

    char defaultAttr[] = "NoName";
    agsafeset(n, EDGE_ATTR_LABEL, const_cast<char *>(name.c_str()), defaultAttr);// set name as preliminary label // should happen by default?!
    //agsafeset(n, "fontsize", "24", "24");// default fontsize; TODO set globally somewhere

    if(recordTypeSize) {
        if(userData.indexOf(recordType) < 0) userData.push_back(recordType);
        return agbindrec(n, const_cast<char *>(recordType), recordTypeSize, FALSE);
    }

    return 0;
}

/**
 * @brief Set a QWidget to be drawn within the node.
 * @param node Node identifier.
 * @param w The widget to add.
 */
void Graph::setNodeGraphicsItem(std::string node, QGraphicsItem *w)
{
    Agnode_t *n = getNode(node);
    nodeWidgets[n] = w;
    agsafeset(n, NODE_ATTR_WIDTH,  const_cast<char *>(QString::number(M_SQRT2*2*w->boundingRect().width() / physicalDpiX).toStdString().c_str()), const_cast<char *>("")); // need to set in inch
    agsafeset(n, NODE_ATTR_HEIGHT, const_cast<char *>(QString::number(M_SQRT2*2*w->boundingRect().height() / physicalDpiY).toStdString().c_str()), const_cast<char *>(""));
}

/**
 * @brief Convenience function for void Graph::addEdge(std::string nn1, std::string nn2, std::string en, QStyle *style).
 */
void Graph::addEdge(QString nn1, QString nn2, QString en, QStyle *style)
{
    addEdge(nn1.toStdString(), nn2.toStdString(), en.toStdString(), style);
}

/**
 * @brief Adds an undirected edge to the graph.
 * @param nn1 Identifier of first node.
 * @param nn2 Identifier of second node.
 * @param en Identifier of the edge to be added.
 * @param style Optional style for the edge.
 */
void Graph::addEdge(std::string nn1, std::string nn2, std::string en, QStyle *style)
{
    Agnode_t *n1 = getNode(nn1);
    if(!n1)
        throw GraphVizQTException("Invalid node:", nn1);

    Agnode_t *n2 = getNode(nn2);
    if(!n2)
        throw GraphVizQTException("Invalid node:", nn2);

    // std::cerr<<"Adding edge <"<<nn1<<"> to <"<<nn2<<">"<<std::endl;

    Agedge_t * e = agedge(g, n1, n2, const_cast<char *>(en.c_str()), TRUE);
    char *defaultAttr = const_cast<char *>(en.c_str());
    agsafeset(e, EDGE_ATTR_LABEL, const_cast<char *>(en.c_str()), defaultAttr);
}

/**
 * @brief Convenience function for Graph::getNodeRecord(node_t *n, char *recordType)
 */
void *Graph::getNodeRecord(QString nn, char *recordType)
{
    return getNodeRecord(nn.toStdString(), recordType);
}

/**
 * @brief Convenience function for Graph::getNodeRecord(node_t *n, char *recordType)
 */
void *Graph::getNodeRecord(std::string nn, char *recordType)
{
    return getNodeRecord(getNode(nn), recordType);
}

/**
 * @brief Retrieve the user record of type recordType associated with node n. @sa addNode().
 * @param n The node.
 * @param recordType Which user-record?
 * @return Pointer to user-record.
 */
void *Graph::getNodeRecord(node_t *n, char *recordType)
{
    return aggetrec(n, recordType, 0); // TODO : check moveToFront param
}

/**
 * @brief Create new undirected graph.
 */
void Graph::newGraph()
{
    // stupid workaround; can't recreate graph, thus no initialization in constructor || or need different constructor
    char name[] = "network";
    g = agopen(name, Agundirected, 0);
    // TODO clear stuff?
}

/**
 * @brief Get number of nodes in this graph.
 * @return Number of nodes.
 */
int Graph::nNodes()
{
    return agnnodes(g);
}

/**
 * @brief Get number of edges in this graph.
 * @return Number of edges.
 */
int Graph::nEdges()
{
    return agnedges(g);
}

/**
 * @brief Draw this graph onto a QGraphicsScene.
 * @param gs Scene to draw on.
 */
void Graph::draw(QGraphicsScene *gs)
{
    if(!g)
        return;

    if(!layouted)
        layout();

    Agnode_t *n;
    Agedge_t *e;
    for (n = agfstnode(g); n; n = agnxtnode(g, n)) {
        for(e = agfstout(g, n); e; e = agnxtout(g, e)) {
            drawEdge(gs, e);
        }
        drawNode(gs, n);
    }
}

/**
 * @brief Layout the graph using GraphViz engines.
 * !!Do not call before previous layouting process has finished (seg fault)!!
 * TODO: lock
 * @param engine Layout engine to use.
 * @return See gvLayout.
 */
int Graph::layout()
{
    std::cout<<"Layouting using " << layoutEngine.toStdString()<<"... ";

    int res = gvLayout(gvc, g, layoutEngine.toStdString().c_str());
    // dpi = QString(agget(g, "dpi")).toDouble(); // is zero?+
    // std::cerr<<"dpi: "<<dpi;

    std::cout<<"done."<<std::endl;

    return res;
}

/**
 * @brief Render the graph to image file.
 * @param format Image file format to use (svg, pdf, png, ...).
 * @param filename Filename for the output.
 * @return See gvRender.
 */
int Graph::render(QString format, QString filename)
{
    //TODO need to check if layouted?
    FILE *fp = fopen(filename.toStdString().c_str(), "w");
    if(fp) {
        gvRender(gvc, g, format.toStdString().c_str(), fp);
    } else {
        std::cerr<<"Could not open "<<filename.toStdString() << "\n";
    }
    fclose(fp);
}

/**
 * @brief Reset graph and read from dot-file
 * @param filename File to be read.
 * @bug works only if no graph was created before.
 */
void Graph::readFromDotFile(QString filename)
{

    if(g) agclose(g);
    FILE *fp = fopen(filename.toStdString().c_str(), "r");
    std::cout<<"Reading "<<filename.toStdString()<<"... ";
    if(!(g = agread(fp, NULL))){
        std::cerr<<"Could not read graph."<<std::endl;
    }
    fclose(fp);
    std::cout<<"done."<<std::endl;
}

/**
 * @brief Convenience function for Graph::setGraphAttribute(QString a, QString v)
 */
void Graph::setGraphAttribute(QString a, QString v)
{
    setGraphAttribute(a.toStdString(), v.toStdString());
}

/**
 * @brief Set GraphViz attributes on graph. See GraphViz doc for more information. Not all attributes supported yet.
 * @param a Attribute.
 * @param v Value.
 */
void Graph::setGraphAttribute(std::string a, std::string v)
{
    agsafeset(g, const_cast<char *>(a.c_str()), const_cast<char *>(v.c_str()), const_cast<char *>(v.c_str()));//const_cast<char *>(defaultGraphAttributes[a].c_str())); // TODO default?
}

/**
 * @brief Get bounding box for node n.
 * @param n Node.
 * @return Bounding box.
 */
QRectF Graph::getNodeRectF(Agnode_t *n)
{
    // TODO check if layouted?
    pointf c = ND_coord(n); // in points
    double w = ND_width(n); // in inches
    double h = ND_height(n); // in inches

#ifndef _WINDOWS_
    QRectF rf = QRectF(c.x, c.y, w/2 * physicalDpiX, h/2 * physicalDpiY);
#endif
#ifdef _WINDOWS_
    QRectF rf = QRectF(c.x, c.y, w * 0.75 * physicalDpiX, h * 0.75 * physicalDpiY); // magic factor 0.75... otherwise gap between nodes and edges
#endif

    // coords are topleft of bounding box
    rf.moveTo(rf.x() - 0.5 * rf.width(), rf.y() - 0.5 * rf.height());
    return rf;
    // TODO edgelabels
}

/**
 * @brief Draw the given edge onto a QGraphicsScene.
 * @param gs The scene.
 * @param e The edge.
 */
void Graph::drawEdge(QGraphicsScene *gs, Agedge_t *e)
{
    // edge label
    char *label = agget(e, EDGE_ATTR_LABEL);
    char *tooltip = agget(e, EDGE_ATTR_TOOLTIP);

    splines *spl = ED_spl(e);
    // color
    QPen pen = QPen(getColor(e, EDGE_ATTR_COLOR));

    // width
    if(double penwidth = QString(agget(e, EDGE_ATTR_PENWIDTH)).toDouble())
        pen.setWidthF(penwidth);
    // pen.setWidthF(pen.widthF() * 4);

    // edge style
    char *style = agget(e, EDGE_ATTR_STYLE);
    if(style) {
        std::string sstyle = style;
        if(sstyle == "bold") {
            pen.setWidthF(EDGE_STYLE_BOLD_SCALING * pen.widthF());
        } else if(sstyle == "dotted") {
            pen.setStyle(Qt::DotLine);
        } else if(sstyle == "dashed") {
            pen.setStyle(Qt::DashLine);
        // } else if(sstyle == "solid") { default
        }
    }

    // Connections are combination of splines
    for(int sc = 0; sc < spl->size; ++sc) {
        // graphViz bezier holds starting point + n*3 points for cubic bezier
        bezier b = spl->list[sc];
        assert(b.size % 3 == 1);

        if(sc == 0) { // draw label to center of first spline
            QPointF labelPos;

            if(b.size % 2) {
                int idx = (b.size - 1)/ 2;
                labelPos.setX(b.list[idx].x);
                labelPos.setY(b.list[idx].y);
            } else {
                int idx = b.size / 2;
                labelPos.setX((b.list[idx].x + b.list[idx - 1].x) / 2);
                labelPos.setY((b.list[idx].y + b.list[idx - 1].y) / 2);
            }

            // Edge label
            if(label) {
                QGraphicsSimpleTextItem *labelTextItem = new QGraphicsSimpleTextItem(QString(label));
                labelTextItem->setPos(labelPos);
                QFont labFont = labelTextItem->font();
                if(int fontsize = QString(agget(e, EDGE_ATTR_LABELFONTSIZE)).toInt())
                    labFont.setPointSize(fontsize);
                labelTextItem->setPen(QPen(getColor(e, EDGE_ATTR_LABELFONTCOLOR)));
                labelTextItem->setToolTip(tooltip);
                gs->addItem(labelTextItem);
            }
        }

        QPainterPath spl;
        spl.moveTo(b.list[0].x, b.list[0].y); // starting point
        ++b.list;

        for(int i = 1; i < b.size; i += 3) {
            //pointf sp = b.sp; // starting point, (0,0) if undirected
            spl.cubicTo(b.list[0].x, b.list[0].y,
                        b.list[1].x, b.list[1].y,
                        b.list[2].x, b.list[2].y);
            b.list += 3;
        }
        gs->addPath(spl, pen);
    }
}

/**
 * @brief Draw the given node onto a QGraphicsScene.
 * @param gs The scene.
 * @param n The node.
 */
void Graph::drawNode(QGraphicsScene *gs, Agnode_t *n)
{
    QString label = agget(n, const_cast<char *>("label"));
    QString shape = agget(n, const_cast<char *>("shape"));
    //std::cerr<<shape.toStdString()<<std::endl;

    QRectF nr = getNodeRectF(n);

    QColor strokeColor = getColor(n, "color");
    QColor fillColor = getColor(n, "fillcolor");
    QColor fontColor = getColor(n, "fontcolor");

    QGraphicsEllipseItem *el = new QGraphicsEllipseItem(0, 0, nr.width(), nr.height()); // cannot use QGraphicsEllipseItem(nr); does not set size on item correctly, problem to retrieve item onclick later
    el->setPos(nr.topLeft());
    el->setToolTip(label);
    el->setBrush(QBrush(fillColor));
    el->setPen(QPen(strokeColor));

    // Add user-type data as QVariant to GraphicsItem
    foreach(std::string o, userData) {
        void *d = getNodeRecord(n, const_cast<char *>(o.c_str()));
        if(d) {
            el->setData(userData.indexOf(o), QVariant::fromValue(d));
        }
    }

    gs->addItem(el);

    // Draw Widget?
    std::map<Agnode_t*, QGraphicsItem*>::const_iterator it;
    if((it = nodeWidgets.find(n)) != nodeWidgets.end()) {
        // Yes, widget
        QGraphicsItem *w = it->second;

        // center in node shape
        QPoint newCenter = nr.center().toPoint();
#ifndef _WINDOWS_
        newCenter.setX(newCenter.x() - w->boundingRect().width() / 2);
        newCenter.setY(newCenter.y() - w->boundingRect().height() / 2);
#endif
#ifdef _WINDOWS_
        // magic windows correction to avoid node-edge-gap
        w->setScale(1/0.75);
        newCenter.setX(newCenter.x() - w->boundingRect().width() / 0.75 / 2);
        newCenter.setY(newCenter.y() - w->boundingRect().height() / 0.75 / 2);
#endif
        w->setPos(newCenter);
        gs->addItem(w);

    } else {
        // No widget, text only

        QGraphicsTextItem *li = new QGraphicsTextItem();
        // TODO font color and formatting li->setFont();
        li->setHtml(label);
        QFont labelFont = li->font();
        // TODO font attributes bold, ...
        if(int fontsize = QString(agget(n, const_cast<char *>("fontsize"))).toInt())
            labelFont.setPointSize(fontsize);
        li->setFont(labelFont);
        li->setPos(nr.center().x() - 0.5 * li->boundingRect().width(), nr.center().y() - 0.5 * li->boundingRect().height());
        QBrush fgBrush;
        fgBrush.setColor(fontColor); // todo: check
        gs->setForegroundBrush(fgBrush);
        gs->addItem(li);
    }
}

/**
 * @brief See Graph::setEdgeAttribute(std::string nn1, std::string nn2, std::string en, std::string a, std::string v)
 */
void Graph::setEdgeAttribute(QString nn1, QString nn2, QString en, QString a, QString v)
{
    setEdgeAttribute(nn1.toStdString(), nn2.toStdString(), en.toStdString(), a.toStdString(), v.toStdString());
}

/**
 * @brief Set GraphViz edge-attribute a to value v on edge between nodes nn1 and nn2 with name en.
 * @param nn1
 * @param nn2
 * @param en
 * @param a
 * @param v
 */
void Graph::setEdgeAttribute(std::string nn1, std::string nn2, std::string en, std::string a, std::string v)
{
    // TODO keep map of edge objects and names?
    Agnode_t *n1 = agnode(g, const_cast<char *>(nn1.c_str()), FALSE);
    if(!n1)
        throw GraphVizQTException("Invalid node:", nn1);
    Agnode_t *n2 = agnode(g, const_cast<char *>(nn2.c_str()), FALSE);
    if(!n2)
        throw GraphVizQTException("Invalid node:", nn2);

    Agedge_t *e = agedge(g, n1, n2, const_cast<char *>(en.c_str()), FALSE);
    agsafeset(e, const_cast<char *>(a.c_str()), const_cast<char *>(v.c_str()), const_cast<char *>(defaultNodeAttributes[a].c_str())); // TODO default?
}


/**
 * @brief See Graph::setNodeAttribute(std::string nn, std::string a, std::string v)
 */

void Graph::setNodeAttribute(QString nn, QString a, QString v)
{
    setNodeAttribute(nn.toStdString(), a.toStdString(), v.toStdString());
}

/**
 * @brief Set GraphViz node-attribute a to value v for node named nn.
 * @param nn
 * @param a
 * @param v
 */
void Graph::setNodeAttribute(std::string nn, std::string a, std::string v)
{
    Agnode_t *n = agnode(g, const_cast<char *>(nn.c_str()), FALSE);
    agsafeset(n, const_cast<char *>(a.c_str()), const_cast<char *>(v.c_str()), const_cast<char *>(defaultNodeAttributes[a].c_str())); // TODO default?
}

/**
 * @brief Set GraphViz node-attribute a to value v for node n.
 * @param n
 * @param a
 * @param v
 */
void Graph::setNodeAttribute(Agnode_t *n, std::string a, std::string v)
{
    agsafeset(n, const_cast<char *>(a.c_str()), const_cast<char *>(v.c_str()), const_cast<char *>(defaultNodeAttributes[a].c_str())); // TODO default?
}

/**
 * @brief User-data-type to index mapping. Used for QVariants when adding user-type data to GraphicsItems.
 * @param type Typename.
 * @return The index for QGraphicsItem::data(int).
 */
int Graph::getDataIndex(std::string type)
{
    return userData.indexOf(type);
}

QString Graph::getLayoutEngine() const
{
    return layoutEngine;
}

void Graph::setLayoutEngine(QString e)
{
    if(e == layoutEngine)
        return;

    layoutEngine = e;
    layouted = false;
}

/**
 * @brief Get color for attribute attr of the given GraphViz object (node, edge, ...) o and convert to QColor.
 * @param o Object.
 * @param attr Attribute
 * @return Color as QColor.
 */
QColor Graph::getColor(void *o, QString attr)
{
    QString cs  = agget(o,  const_cast<char *>(attr.toStdString().c_str()));
    if(cs.size()) {
        if(cs.startsWith('#')) {
            // #RRGGBB
            return QColor::fromRgb(cs.right(6).toUInt(0, 16));
        }
        return QColor(cs);
    }

    // default values // todo check default attribute maps
    if(attr == EDGE_ATTR_COLOR) return QColor("black");
    if(attr == EDGE_ATTR_LABELFONTCOLOR) return QColor("black");
    if(attr == NODE_ATTR_COLOR) return QColor("black");
    if(attr == NODE_ATTR_FILLCOLOR) return QColor("lightgray");
    if(attr == NODE_ATTR_FONTCOLOR) return QColor("black");

    return QColor("black");
}

/**
 * @brief Initialize internal variables.
 */
void Graph::init()
{
    //gvAddLibrary(gvc, &gvplugin_dot_layout_LTX_library);
    //gvAddLibrary(gvc, &gvplugin_neato_layout_LTX_library);
    //gvAddLibrary(gvc, &gvplugin_core_LTX_library);

    agseterr(AGWARN);
    layouted = false;
    layoutEngine = GRAPHVIZQT_H_DEFAULT_LAYOUT_ENGINE;

    // Determine screen resolution
    QDesktopWidget dw;
    physicalDpiX = dw.physicalDpiX();
    physicalDpiY = dw.physicalDpiX();

    defaultEdgeAttributes[EDGE_ATTR_COLOR] = "black";
    defaultEdgeAttributes[EDGE_ATTR_PENWIDTH] = "1  ";
    defaultNodeAttributes[NODE_ATTR_COLOR] = "black";
    defaultNodeAttributes["fontsize"] = "12"; // TODO get defaults from graphviz?
}

/**
 * @brief Get GraphViz node object for given node name.
 * @param node Node identifier.
 * @return GraphViz node object.
 */
Agnode_t *Graph::getNode(std::string node)
{
    return agnode(g, const_cast<char *>(node.c_str()), FALSE);
}

/** @brief See Graph::getNode(std::string node). */
Agnode_t *Graph::getNode(QString node)
{
    return agnode(g, const_cast<char *>(node.toStdString().c_str()), FALSE);
}

GraphVizQTException::GraphVizQTException()
{
    std::cerr<<"GraphVizQTException"<<std::endl;
}

GraphVizQTException::GraphVizQTException(std::string e)
{
    std::cerr<<"GraphVizQTException: "<<e<<std::endl;
}

GraphVizQTException::GraphVizQTException(std::string e1, std::string e2)
{
    std::cerr<<"GraphVizQTException: "<<e1<<e2<<std::endl;
}

}
