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

#ifndef GRAPHVIZQT_H
#define GRAPHVIZQT_H

#include<map>
#include<iostream>
#include<exception>

#include<gvc.h>

#include<QString>
#include<QRectF>
#include<QGraphicsScene>
#include<QGraphicsItem>
#include<QList>

#include "gvplugin.h"

namespace graphvizqt {

class GraphVizQTException;

class GraphVizQT
{
public:
    GraphVizQT();
};

#define GRAPHVIZQT_H_DEFAULT_LAYOUT_ENGINE "dot"

/**
 * @brief The Graph class. QT-Graphviz interface. Draws a graph layouted by GraphViz onto a QGraphicsScene.
 */

class Graph : public QObject
{
    Q_OBJECT

    //TODO use map<id, ...> ; *edgePainter ; *nodePainter

public:
    Graph();
    Graph(std::string dot);
    Graph(QString dot);

    ~Graph();

    void *addNode(QString name, const char *recordType = "", size_t recordTypeSize = 0, QStyle *style = 0);
    void *addNode(std::string name, const char *recordType = "", size_t recordTypeSize = 0, QStyle *style = 0);
    void setNodeGraphicsItem(std::string node, QGraphicsItem *w);
    void addEdge(QString n1, QString n2, QString e, QStyle *style = 0);
    void addEdge(std::string n1, std::string n2, std::string e, QStyle *style = 0);
    void *getNodeRecord(QString nn, char *recordType = 0);
    void *getNodeRecord(std::string nn, char *recordType = 0);
    void *getNodeRecord(node_t *n, char *recordType = 0);

    void newGraph();
    int nNodes();
    int nEdges();
    // int nSubgraphs(); /** Number of subgraphs */
    void draw(QGraphicsScene *gs);
    int layout();
    int render(const QString format, const QString filename);
    void readFromDotFile(QString filename);
    void setGraphAttribute(QString a, QString v);
    void setGraphAttribute(std::string a, std::string v);
    void setNodeAttribute(QString nn, QString a, QString v);
    void setEdgeAttribute(QString nn1, QString nn2, QString en, QString a, QString v);
    void setNodeAttribute(std::string nn, std::string a, std::string v);
    void setNodeAttribute(Agnode_t *n, std::string a, std::string v);
    void setEdgeAttribute(std::string nn1, std::string nn2, std::string en, std::string a, std::string v);
    int getDataIndex(std::string type);

    QString getLayoutEngine() const;

public slots:
    void setLayoutEngine(QString e);

private:
    QRectF getNodeRectF(Agnode_t * n);
    void drawEdge(QGraphicsScene* gs,  Agedge_t *e);
    void drawNode(QGraphicsScene* gs, Agnode_t *n);
    Agnode_t *getNode(std::string node);
    Agnode_t *getNode(QString node);
    //void setStyle(Agedge_t *n, QStyle *style = 0);
    //void setStyle(Agnode_t *e, QStyle *style = 0);
    QColor getColor(void *o, QString a);
    void init();

    static const float EDGE_STYLE_BOLD_SCALING;

    // make not const to avoid const_cast<char *> for all graphviz functions
    static char EDGE_ATTR_LABEL[];
    static char EDGE_ATTR_COLOR[];
    static char EDGE_ATTR_STYLE[];
    static char EDGE_ATTR_LABELFONTSIZE[];
    static char EDGE_ATTR_LABELFONTCOLOR[];
    static char EDGE_ATTR_TOOLTIP[];
    static char EDGE_ATTR_PENWIDTH[];

    static char NODE_ATTR_COLOR[];
    static char NODE_ATTR_FILLCOLOR[];
    static char NODE_ATTR_FONTCOLOR[];
    static char NODE_ATTR_WIDTH[];
    static char NODE_ATTR_HEIGHT[];

    QList<std::string> userData; /**< Index for Qt setData(int, QVariant) functions. To bind user-types to nodes in graphicsscene. */
    std::map<std::string, std::string> defaultEdgeAttributes; // TODO make static
    std::map<std::string, std::string> defaultNodeAttributes;
    std::map<std::string, std::string> defaultGraphAttributes;

    std::map<Agnode_t*, QGraphicsItem*> nodeWidgets; /**< QWidgets that were added to nodes */

    GVC_t *gvc;             /**< GraphViz context. */
    Agraph_t *g;            /**< The GraphViz graph. */
    bool layouted;          /**< Is the graph layouted and can be plotted?. */
    QString layoutEngine;
    double physicalDpiX;    /**< Screen resolution (DPI) for converting GraphViz coordinates and distances. */
    double physicalDpiY;    /**< Screen resolution (DPI) for converting GraphViz coordinates and distances. */

};


class GraphVizQTException : public std::exception {
public:
    GraphVizQTException();
    GraphVizQTException(std::string e);
    GraphVizQTException(std::string e1, std::string e2);
};


}
#endif // GRAPHVIZQT_H
