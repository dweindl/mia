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

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

#include <QHBoxLayout>
#include <QtConcurrent>
#include <QFuture>
#include <QTreeView>

#ifdef MIA_WITH_HDF5
    #include <H5Cpp.h>
    #include "hdfwriter.h"
#endif

#include "spectrum.h"
#include "labeledcompound.h"
#include "libraryhit.h"

#include "miamainwindow.h"
#include "src/misc.h"
#include "src/settings.h"
#include "src/miaexception.h"
#include "nodecompoundtreemodel.h"
#include "experimentwizard.h"
#include "miaguiconstants.h"
#include "midplot.h"
#include "configdialog.h"

#ifdef MIA_WITH_METABOBASE
    #include "mddb_pgsql.h"
#endif

namespace mia {

MIAMainWindow::MIAMainWindow(QWidget *parent) :
    QMainWindow(parent)
{
    init();
    setupGUI();
}

MIAMainWindow::~MIAMainWindow()
{
    delete view;
    delete scene;
    delete LabelingNetworkSet::distCalc;
    delete g;

#ifdef MIA_WITH_NETCDF_IMPORT
    if(netCDFImportDialog)
        delete netCDFImportDialog;
#endif

    delete networkSet;

    if(excludeLib) delete excludeLib;
}

/**
 * @brief Start compound detection. (static function for QtConcurrent::run())
 * @param ds
 */
void MIAMainWindow::findLabeledCompounds(mia::NetworkLayer *ds)
{
    std::cout<<ds->getSettings().experiment<<": Starting compound detection."<<std::endl;

    qDebug("%s", ds->getSettings().toString().c_str());

    ds->findLabeledCompounds();

    std::cout<<ds->getSettings().experiment<<": Found "<<ds->cmpLab.size()<<" labeled compounds and "
            <<ds->cmpUnlab.size()<< " unlabeled compounds"<<std::endl<<std::endl;
}

void MIAMainWindow::recreateGraph()
{
    emit(progressText("Creating graph."));

    networkSet->setExcludeM0(qsettings.value("alignment_m0", 0).toInt());
    updateCompoundList();
    setupExperimentOverlayGraph();
}

/*
void MIAMainWindow::repaintGraph()
{
    // Create graph from selected MIDs in compound panel
    scene->clear();
    g->draw(scene);
}
*/

void MIAMainWindow::labelDetectionFinished()
{
    static bool waiting = false;

    if(waiting)
        return; // there is already another call waiting

    waiting = true;

    // all threads finished?
    foreach(LabelIdentificatorQThread *t, labidThreads) {
        t->wait();
    }
    while (!labidThreads.isEmpty())
        delete labidThreads.takeFirst();

    matchCompoundsAcrossExperiments(); // TODO separate thread?
    matchCompoundsAgainstLibrary();
    experimentListWidget->updateExperimentList(networkSet->getDatasets(), experimentColors);
    updateCompoundList();
    networkSet->createDistanceMatrices();
    recreateGraph();

    if(progressDialog) {
        delete progressDialog;
        progressDialog = 0;
    }

    QApplication::processEvents(); // first process finished() of other threads

    waiting = false;
}


void MIAMainWindow::compoundClicked(QModelIndex mi)
{
    // find corresponding widget to center graphics view on it
    int cmpID = compoundTreeView->model()->data(mi, Qt::UserRole).toInt();
    NodeWidget *w = nodeWidgets[cmpID];
    if(w)
        view->view()->centerOn(w);
}

/**
 * @brief Setup some class variables
 */
void MIAMainWindow::init()
{
    // overlay colors
    experimentColors.push_back(QColor(0xCD, 0x16, 0x16)); // red
    experimentColors.push_back(QColor(0x22, 0x8A, 0x27)); // green
    experimentColors.push_back(QColor(0x1B, 0x54, 0xD8)); // blue
    experimentColors.push_back(QColor(0xEB, 0xC8, 0x16)); // yellow
    experimentColors.push_back(QColor(0x81, 0x00, 0x81)); // violett
    experimentColors.push_back(QColor(0xB1, 0x32, 0x44)); // darkred
    experimentColors.push_back(QColor(0x89, 0xa3, 0x5c)); // some green
    experimentColors.push_back(QColor(0x00, 0xff, 0xff)); // cyan
    experimentColors.push_back(QColor(0x80, 0x80, 0x00)); // dark yellow
    experimentColors.push_back(QColor(0x00, 0x00, 0x80)); // dark blue

    experimentColors.push_back(QColor(0xCD, 0x16, 0x16)); // red
    experimentColors.push_back(QColor(0x22, 0x8A, 0x27)); // green
    experimentColors.push_back(QColor(0x1B, 0x54, 0xD8)); // blue
    experimentColors.push_back(QColor(0xEB, 0xC8, 0x16)); // yellow
    experimentColors.push_back(QColor(0x81, 0x00, 0x81)); // violett
    experimentColors.push_back(QColor(0xB1, 0x32, 0x44)); // darkred
    experimentColors.push_back(QColor(0x89, 0xa3, 0x5c)); // some green
    experimentColors.push_back(QColor(0x00, 0xff, 0xff)); // cyan
    experimentColors.push_back(QColor(0x80, 0x80, 0x00)); // dark yellow
    experimentColors.push_back(QColor(0x00, 0x00, 0x80)); // dark blue
// TODO: autogenerate... rainbow()?

    g = new graphvizqt::Graph();

    LabelingNetworkSet::distCalc = new MIDDistanceCalculator(NW_GAP_PENALTY);
    graphSizeWarningLimit = 200;
    excludeLib = 0;
    progressDialog = 0;
    networkSet = new LabelingNetworkSet();
#ifdef MIA_WITH_NETCDF_IMPORT
    netCDFImportDialog = 0;
#endif
#ifdef MIA_WITH_METABOBASE
    keggMapper = new KEGGReactionMapper();
    keggMapper->setupGraph();
#endif
}

/**
 * @brief Match spectra from the different experiments
 * Populate @nodes maps
 */
void MIAMainWindow::matchCompoundsAcrossExperiments()
{
    emit(progressText("Matching compounds across datasets."));

    qDeleteAll(nodeWidgets);
    nodeWidgets.clear();

    double mylibScoreCutoff = qsettings.value("cmp_matching_score_cutoff", CMP_MATCHING_SCORE_CUTOFF).toDouble(); // TODO use from Settings
    bool useLargestCommonIon = qsettings.value("nw_use_common_largest_ion", NW_USE_LARGEST_COMMON_ION).toBool();
    networkSet->matchCompoundsAcrossExperiments(mylibScoreCutoff, useLargestCommonIon);
}

void MIAMainWindow::startLabelDetection(NetworkLayer *ds)
{
    if(!progressDialog) {
        progressDialog = new QProgressDialog(this);
        progressDialog->setModal(true);
        progressDialog->setMinimumDuration(0);
        progressDialog->setMaximum(0);
        connect(this, SIGNAL(progressText(QString)), progressDialog, SLOT(setLabelText(QString)));
        progressDialog->show();
    }

    std::cout<<"Starting '"<<ds->getSettings().experiment<<"'\n";
    qDebug("%s", ds->getSettings().toString().c_str());
    emit(progressText(QString("Detecting labeled compounds for <%1>.").arg(QString::fromStdString(ds->getSettings().experiment))));

    LabelIdentificatorQThread *labidThread = new LabelIdentificatorQThread(ds, this);
    labidThread->start(QThread::LowPriority);
    labidThreads.push_back(labidThread);

    connect(labidThread, SIGNAL(progressMessage(QString)), progressDialog, SLOT(setLabelText(QString)));
//        connect(labidThread, SIGNAL(progressMax(int)),progressDialog,SLOT(setMaximum(int)));
//        connect(labidThread, SIGNAL(progress(int)),progressDialog,SLOT(setValue(int)));
    connect(labidThread, SIGNAL(finished()), this, SLOT(labelDetectionFinished()));
}

void MIAMainWindow::generateNodeWidgets()
{
    const QMap<int, NodeCompound*> nodes = networkSet->getNodeCompounds();
    for(int n = 0; n < nodes.size(); ++n) {
        NodeCompound *nc = nodes[n];

        // colors for the different experiments
        QList<QColor> colors;
        std::vector<std::string> exps = nc->getExperiments();
        QList<NetworkLayer *> datasets = networkSet->getDatasets();
        for(int ds = 0; ds < datasets.size(); ++ds) {
            if (std::find(exps.begin(), exps.end(), datasets[ds]->getSettings().experiment) != exps.end())
                colors.push_back(experimentColors[ds]);
        }

        NodeWidget *w = new NodeWidget(nc, colors);
        nodeWidgets[n] = w;
    }
}


void MIAMainWindow::matchCompoundsAgainstLibrary(QString libFile, bool overwriteNames)
{
    emit(progressText("Matching compounds against library."));

    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    networkSet->matchCompoundsAgainstLibrary(libFile, overwriteNames);
    QApplication::restoreOverrideCursor();
}

void MIAMainWindow::setupCompoundList()
{
    compoundListDockWidget = new QDockWidget(tr("Compounds"), this);
    QWidget *w = new QWidget(compoundListDockWidget);
    QVBoxLayout *layout = new QVBoxLayout(w);
    QToolBar *compoundToolBar = new QToolBar(w);

    // toolbar
    QAction *actExpand = new QAction(QIcon(":/gui/icons/zoom-in.png"), tr("&ExpandAll"), this);
    actExpand->setToolTip("Expand list");
    compoundToolBar->addAction(actExpand);
    QAction *actCollapse = new QAction(QIcon(":/gui/icons/zoom-out.png"), tr("&CollapseAll"), this);
    actCollapse->setToolTip("Collapse list");
    compoundToolBar->addAction(actCollapse);
/*
    QAction *actRefreshGraph = new QAction(QIcon::fromTheme("system-run"), tr("&Recreate graph"), this);
    actRefreshGraph->setToolTip("Recreate graph with selected MIDs.");
    compoundToolBar->addAction(actRefreshGraph);
*/
    layout->addWidget(compoundToolBar);

    // compound tree itself
    compoundTreeView = new QTreeView;
    compoundTreeView->setSelectionMode(QAbstractItemView::ExtendedSelection);

    connect(compoundTreeView, SIGNAL(clicked(QModelIndex)), this, SLOT(compoundClicked(QModelIndex)));
    connect(actExpand, SIGNAL(triggered()), compoundTreeView, SLOT(expandAll()));
    connect(actCollapse, SIGNAL(triggered()), compoundTreeView, SLOT(collapseAll()));
//    connect(actRefreshGraph, SIGNAL(triggered()), this, SLOT(repaintGraph()));
    layout->addWidget(compoundTreeView);

    w->setLayout(layout);
    compoundListDockWidget->setWidget(w);
}

void MIAMainWindow::setupGUI()
{
    setupToolbar();

    setupGraphOptionPanel();
    setupExperimentList();
    setupCompoundList();

    addDockWidget(Qt::RightDockWidgetArea, graphOptionsDockWidget);
    addDockWidget(Qt::RightDockWidgetArea, experimentListDockWidget);
    tabifyDockWidget(graphOptionsDockWidget, experimentListDockWidget);
    addDockWidget(Qt::LeftDockWidgetArea, compoundListDockWidget);

    // main canvas
    scene = new QGraphicsScene();
    view = new NWView(scene, g, this);

    // put all together
    setCentralWidget(view);

    // add test graph
    // graphvizqt::Graph g(std::string("graph{A--B;A--C;A--D;A--B;B--C; C--E; E--F; A[label=\"AAA\"]};"));
    // g->readFromDotFile("networkoverlaytest.dot");
}

void MIAMainWindow::setupToolbar()
{
    QToolBar *mainTB = new QToolBar(this);

#ifdef MIA_WITH_NETCDF_IMPORT
    QAction *actImportData = new QAction(QIcon(":/gui/icons/document-import.png"), tr("&Import data"), mainTB);
    actImportData->setToolTip("Import netCDF data or redetect compounds");
    mainTB->addAction(actImportData);
    connect(actImportData, SIGNAL(triggered()), this, SLOT(showDataImportDialog()));
#endif

    QAction *actOpenXML = new QAction(QIcon(":/gui/icons/document-open-xml.png"), tr("Load experiment from XML config"), mainTB);
    actOpenXML->setShortcut(QKeySequence::Open);
    actOpenXML->setToolTip("Load experiment from XML config");
    mainTB->addAction(actOpenXML);

    /* // save and open binary mia-data
    QAction *actOpen = new QAction(QIcon::fromTheme("document-open"), tr("&Open data"), mainTB);
    actOpen->setToolTip("Open data");
    mainTB->addAction(actOpen);

    QAction *actSave = new QAction(QIcon::fromTheme("document-save"), tr("&Save"), mainTB);
    actSave->setToolTip("Save data");
    mainTB->addAction(actSave);
 */

#ifdef MIA_WITH_HDF5
    QAction *actExportHDF = new QAction(QIcon(":/gui/icons/document-save-hdf5.png"), tr("&Save HDF5"), mainTB);
    actExportHDF->setToolTip("Export distances to HDF5 file");
    mainTB->addAction(actExportHDF);
    connect(actExportHDF, SIGNAL(triggered()), this, SLOT(saveHDF()));
#endif

    QAction *actExportMIDs = new QAction(QIcon(":/gui/icons/document-save-csv.png"), tr("&Save CSV"), mainTB);
    actExportMIDs->setToolTip("Export MIDs as comma separated values");
    mainTB->addAction(actExportMIDs);

    QAction *actSelectLibrary = new QAction(QIcon(":/gui/icons/edit-find.png"), tr("&Identify"), mainTB);
    actSelectLibrary->setToolTip("Select compound library for identification");
    mainTB->addAction(actSelectLibrary);

    QAction *actAbout = new QAction(QIcon(":/gui/icons/help-about.png"), tr("Abo&ut"), mainTB);
    actAbout->setToolTip(QString("About ").append(QApplication::applicationName()));
    mainTB->addAction(actAbout);

    QAction *actConfig = new QAction(QIcon(":/gui/icons/preferences-other.png"), tr("&Config"), mainTB);
    actConfig->setToolTip("Show config dialog");
    mainTB->addAction(actConfig);

    addToolBar(mainTB);

    //connect(actSave, SIGNAL(triggered()), this, SLOT(saveFile()));
    connect(actExportMIDs, SIGNAL(triggered()), this, SLOT(exportMIDs()));
    connect(actOpenXML, SIGNAL(triggered()), this, SLOT(openXMLFile()));
    //connect(actOpen, SIGNAL(triggered()), this, SLOT(openFile()));
    connect(actAbout, SIGNAL(triggered()), this, SLOT(showInfoDialog()));
    connect(actConfig, SIGNAL(triggered()), this, SLOT(showConfigDialog()));
    connect(actSelectLibrary, SIGNAL(triggered()), this, SLOT(selectCompoundLibraryAndIdentify()));
}

void MIAMainWindow::setupExperimentList()
{
    // tracer list
    experimentListDockWidget = new QDockWidget(tr("Datasets"), this);
    experimentListWidget = new ExperimentListWidget(this);

    experimentListDockWidget->setWidget(experimentListWidget);

    connect(experimentListWidget, SIGNAL(experimentSettingChanged(int,Settings)), this, SLOT(experimentSettingChanged(int,Settings)));
    connect(experimentListWidget, SIGNAL(experimentAdded(Settings)), this, SLOT(addExperiment(Settings)));
    connect(experimentListWidget, SIGNAL(experimentUsageChanged()), this, SLOT(experimentSelectionChanged()));
    connect(experimentListWidget, SIGNAL(experimentRemoved(NetworkLayer*)), this, SLOT(experimentRemoved(NetworkLayer*)));

}

/**
 * @brief MIAMainWindow::setupTracerOverlayGraph Create tracer overlay graph from distance matrices
 */
void MIAMainWindow::setupExperimentOverlayGraph()
{
    std::cout<<"Setup overlay graph"<<std::endl;

    // collect some settings
    double variationCutoff = variationSlider->value() / 10000.0;
    bool hideLessVarying = hideLessVaryingNodes->checkState() == Qt::Checked;
    bool hideFoundInLessExperiments = this->hideFoundInLessExperiments->checkState() == Qt::Checked;
    bool showUnconnectedNodes = qsettings.value("graph_show_unconnected_nodes", true).toBool();
    int excludeIfFoundInLessExperiments = qsettings.value("exclude_if_found_in_less_experiments", 1).toInt();

    if(!networkSet->getDatasets().size())
        return; // nothing to do

    // count edges
    int e = networkSet->getNumberOfEdges(variationCutoff, excludeIfFoundInLessExperiments);

    // check number of edges
    if(e > graphSizeWarningLimit) {
        if(! showGraphSizeWarning(e))
            return;
    }

    // reclaim NodeWidgets before resetting scene, otherwise they will be destroyed -> delete later manually
    if(nodeWidgets.size()) {
         QList<QGraphicsItem *> items = scene->items();
         foreach (QGraphicsItem* i, items) {
             if(dynamic_cast<NodeWidget *>(i))
                  scene->removeItem(i);
         }
    }

    g->newGraph();

    g->setGraphAttribute(std::string("overlap"), std::string("false"));
    g->setGraphAttribute(std::string("splines"), std::string("true"));

    // add nodes TODO only first time, then just delete and re-add edges
    if(!nodeWidgets.size()) { // TODO chekc if size = same, crash after reload -> implement unload()
        // generate nodeWidgets TODO: check if stuff changed and need to recreate nodes
        generateNodeWidgets();
    }

    // add nodes
    std::map<int, NodeCompound *> visNodes = networkSet->getNodesInGraph(showUnconnectedNodes, hideLessVarying, variationCutoff, hideFoundInLessExperiments, excludeIfFoundInLessExperiments);
    for(std::map<int, NodeCompound *>::iterator it = visNodes.begin();
        it != visNodes.end(); ++it) {

        NodeCompound *nc = it->second;
        int n = it->first;
        std::string nodeName = nc->getCompoundName();

        mynode *newNode = reinterpret_cast<mynode *>(g->addNode(nodeName, "mynode", sizeof(mynode)));
        newNode->nc = nc;

        // 70 RSD
        // 700*sd
        //factor = factor < 0 ? 0 : factor;
        g->setNodeAttribute(nodeName, "fillcolor",
                            QColor::fromRgb(0x41, 0x5E, 0x79).lighter(100 + 1200 * nc->getMaxIsotopomerSD()).name().toStdString());
        g->setNodeGraphicsItem(nodeName, nodeWidgets[n]);
    }

    addEdgesToGraph(excludeIfFoundInLessExperiments, variationCutoff);
    std::cout<<"Added "<<e<<" edges"<<std::endl;

    scene->clear();
    g->draw(scene);
    view->resetView();

    std::cout<<"Setup overlay graph done."<<std::endl;
}

void MIAMainWindow::setupGraphOptionPanel()
{
        graphOptionsDockWidget = new QDockWidget(tr("Graph options"), this);
        QWidget *nwWidget = new QWidget(graphOptionsDockWidget);
        QGridLayout *nwGrid = new QGridLayout(nwWidget);

        // layout engines
        QGroupBox *layoutGroupBox = new QGroupBox("Layout engine", nwWidget);
        QVBoxLayout *vl = new QVBoxLayout(layoutGroupBox);

        QListWidget *engineList = new QListWidget(nwWidget);
        QStringList listItemsEngine;
        listItemsEngine<<"dot"<<"neato"<<"twopi"<<"circo"<<"fdp"<<"sfdp";
        engineList->addItems(listItemsEngine);
        engineList->setCurrentRow(0);
        vl->addWidget(engineList);
        connect(engineList, SIGNAL(currentTextChanged(QString)), g, SLOT(setLayoutEngine(QString)));

        nwGrid->addWidget(layoutGroupBox);

        // distance measure
        QGroupBox *distanceGroupBox = new QGroupBox("Distance calculation", nwWidget);
        vl = new QVBoxLayout(distanceGroupBox);

        QLabel *metricLabel = new QLabel("Metric:", nwWidget);
        vl->addWidget(metricLabel);
        QListWidget *distanceList = new QListWidget(nwWidget);
        QStringList listItems;
        listItems<<"Canberra"<<"Euclidean"<<"Manhattan"<<"Cosine"<<"Custom";
        distanceList->addItems(listItems);
        distanceList->setCurrentRow(0);
        distanceList->sortItems();    
        vl->addWidget(distanceList);
        connect(distanceList, SIGNAL(currentTextChanged(QString)), this, SLOT(distanceMeasureChanged(QString)));

        // normalization measure
        QLabel *normalizerLabel = new QLabel("Normalization:", nwWidget);
        vl->addWidget(normalizerLabel);
        QListWidget *normalizerList = new QListWidget(nwWidget);
        QStringList listItemsNorm;
        listItemsNorm<<"NONE"<<"SUM"<<"PROD"<<"MAX"<<"MIN";
        normalizerList->addItems(listItemsNorm);
        normalizerList->setCurrentRow(1);
        vl->addWidget(normalizerList);
        connect(normalizerList, SIGNAL(currentTextChanged(QString)), this, SLOT(distanceMeasureNormChanged(QString)));

        // distance cutoff
        QLabel* distanceLabel = new QLabel("Distance cutoff:", nwWidget);
        vl->addWidget(distanceLabel);

        cutOffSlider = new QSlider(Qt::Horizontal, nwWidget);
        cutOffSlider->setRange(0, 1000);
        cutOffSlider->setSingleStep(1);
        vl->addWidget(cutOffSlider);
        connect(cutOffSlider, SIGNAL(sliderMoved(int)), this, SLOT(cutOffSliderChanged(int)));

        cutOffLabel = new QLabel("0 %", nwWidget);
        vl->addWidget(cutOffLabel);

#ifdef MIAMAINWINDOW_H_ENABLE_ZSCORE
        useZScore = new QCheckBox("Use ZScore", nwWidget);
        connect(useZScore, SIGNAL(clicked()), this, SLOT(useZScoreChanged()));
        vl->addWidget(useZScore);
#endif
        nwGrid->addWidget(distanceGroupBox);

        // experiment count
        QGroupBox *experimentCountGroupBox = new QGroupBox("Minimum experiments for compounds", nwWidget);
        QHBoxLayout *hl = new QHBoxLayout(experimentCountGroupBox);

        excludeIfFoundInLessExperiments = new QSpinBox(nwWidget);
        excludeIfFoundInLessExperiments->setToolTip("Exclude node from network if present in less than n experiments");
        excludeIfFoundInLessExperiments->setMinimum(1);
        excludeIfFoundInLessExperiments->setValue(qsettings.value("exclude_if_found_in_less_experiments", 1).toInt());
        hl->addWidget(excludeIfFoundInLessExperiments);
        connect(excludeIfFoundInLessExperiments, SIGNAL(valueChanged(int)), this, SLOT(excludeIfFoundInLessExperimentsChanged(int)));

        hideFoundInLessExperiments = new QCheckBox("Hide others", nwWidget);
        connect(hideFoundInLessExperiments, SIGNAL(clicked()), this, SLOT(hideFoundInLessExperimentsChanged()));
        hl->addWidget(hideFoundInLessExperiments);

        nwGrid->addWidget(experimentCountGroupBox);

        // variation
        QGroupBox *variationGroupBox = new QGroupBox("Variation cutoff", nwWidget);
        vl = new QVBoxLayout(variationGroupBox);

        variationSlider = new QSlider(Qt::Horizontal, nwWidget);
        variationSlider->setRange(0, 10000);
        variationSlider->setSingleStep(1);
        variationSlider->setValue(0);
        vl->addWidget(variationSlider);
        connect(variationSlider, SIGNAL(sliderMoved(int)), this, SLOT(variationSliderChanged(int)));

        variationLabel = new QLabel(QString::number(variationSlider->value()), nwWidget);
        vl->addWidget(variationLabel);

        hideLessVaryingNodes = new QCheckBox("Hide others", nwWidget);
        connect(hideLessVaryingNodes, SIGNAL(clicked()), this, SLOT(hideLessVaryingChanged()));
        vl->addWidget(hideLessVaryingNodes);

        nwGrid->addWidget(variationGroupBox);

        graphOptionsDockWidget->setWidget(nwWidget);
}

void MIAMainWindow::addEdgesToGraph(int excludeIfFoundInLessExperiments, double variationCutoff)
{
    // collect distance info for edge line scaling
    double overallMin, overallMax;
    networkSet->getMinMaxDistances(overallMin, overallMax);

    // add edges
    std::vector<LabelingDatasetEdge *> edges = networkSet->getEdges(excludeIfFoundInLessExperiments, variationCutoff);

    for(int i = 0; i < edges.size(); ++i) {
        LabelingDatasetEdge *e = edges[i];

        // label / name
        std::string edgeLabel = "";
        if(qsettings.value("graph_edge_label", "1").toBool()) {
            edgeLabel = QString::number(e->distance).toStdString();
        }

        g->addEdge(e->node1->getCompoundName(), e->node2->getCompoundName(), edgeLabel);

        // weight makes some layouts crash!
        // g->setEdgeAttribute(e->compoundName1, e->compoundName2, edgeLabel,
        //                     "weight", QString::number((int)(1.0 / e->distance)).toStdString());


        // color
        g->setEdgeAttribute(e->node1->getCompoundName(), e->node2->getCompoundName(), edgeLabel,
                            "color", experimentColors[e->datasetIndex % experimentColors.size()].name().toStdString());

        // pen width
        double minPenWidth = 3;
        double maxPenWidth = 30;
        double penWidth = minPenWidth;

        penWidth = maxPenWidth - (maxPenWidth - minPenWidth) * (e->distance - overallMin) / (overallMax - overallMin);

        g->setEdgeAttribute(e->node1->getCompoundName(), e->node2->getCompoundName(), edgeLabel,
                            "penwidth", QString::number(penWidth).toStdString());


#ifdef MIA_WITH_METABOBASE
        //int keggCon = keggMapper->findReactions(nodes[i]->getFeature("PRECURSOR_KEGG_ID"), nodes[j]->getFeature("PRECURSOR_KEGG_ID"), true);
        std::stringstream keggInfo;
        int keggCon = keggMapper->startFindPath(e->node1->getFeature("PRECURSOR_KEGG_ID"), e->node2->getFeature("PRECURSOR_KEGG_ID"), 2, keggInfo);
        if(keggCon) {
            g->setEdgeAttribute(e->node1->getCompoundName(), e->node2->getCompoundName(), edgeLabel,
                                "style", "bold");
        } else {
            g->setEdgeAttribute(e->node1->getCompoundName(), e->node2->getCompoundName(), edgeLabel,
                                "style", "dotted");
        }
        g->setEdgeAttribute(e->node1->getCompoundName(), e->node2->getCompoundName(), edgeLabel,
                            "edgetooltip", keggInfo.str());
#endif

        delete e;
    }
    // g->setGraphAttribute(std::string("splines"), "false");
}

bool MIAMainWindow::showGraphSizeWarning(int edges)
{
    QMessageBox msg;
    std::stringstream str;
    str<<"Graph contains "<<edges<<" edges. This will take a while to lay out. Continue anyway?";
    msg.setWindowTitle("Time...");
    msg.setIcon(QMessageBox::Question);
    msg.setText(QString::fromStdString(str.str()));
    msg.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msg.setDefaultButton(QMessageBox::No);
    msg.setWindowModality(Qt::ApplicationModal);

    int ret = msg.exec();

    return ret == QMessageBox::Yes;
}

void MIAMainWindow::exportMIDs()
{
    QString filename = QFileDialog::getSaveFileName(this, "Save MIDs as CSV", "", "CSV files (*.csv);;All files (*)");

    if(filename.isNull())
        return;

    QFile f(filename);
    f.open(QIODevice::WriteOnly);
    QTextStream out(&f);
    networkSet->exportMIDs(out);
}

void MIAMainWindow::updateCompoundList()
{
    QSortFilterProxyModel *sortProxy = new QSortFilterProxyModel(this);
    sortProxy->setSourceModel(new NodeCompoundTreeModel(networkSet->getNodeCompounds()));
    sortProxy->setDynamicSortFilter(true);

    compoundTreeView->setModel(sortProxy);
    compoundTreeView->setAlternatingRowColors(true);
    compoundTreeView->setSortingEnabled(true);
    compoundTreeView->sortByColumn(1, Qt::AscendingOrder);
    compoundTreeView->resizeColumnToContents(0);
    compoundTreeView->resizeColumnToContents(1);
    compoundTreeView->resizeColumnToContents(2);
    compoundTreeView->resizeColumnToContents(3);
}


void MIAMainWindow::showInfoDialog()
{
    QMessageBox::about(this, QString("About ").append(QApplication::applicationName()),
                       QString("<b>%1</b><br>"
                               "Version %2<br><br>"
                               "2013-2016<br>"
                               "&nbsp;&nbsp;&nbsp;&nbsp;Daniel Weindl &lt;<a href='mailto:sci@danielweindl.de'>sci@danielweindl.de</a>&gt;<br><br>"
                               "<b>Config</b><br>"
                               "<a href='%3'>%3</a><br><br>"
                               "<b>Documentation</b><br>"
                               "<a href='file:%4'>%4</a>"
                               "<a href='%5'>%5</a>"
                               ).arg(QApplication::applicationName(),
                                     QApplication::applicationVersion(),
                                     QSettings().fileName(),
                                 #ifndef _WINDOWS_
                                     "/usr/share/doc/mia/mia-doc.pdf",
                                 #endif
                                 #ifdef _WINDOWS_
                                     QApplication::applicationDirPath().append("/../doc/mia-doc.pdf"),
                                 #endif
                                     "http://massisotopolomeanalyzer.lu/download/mia-doc.pdf"
                                    ));

}

void MIAMainWindow::cutOffSliderChanged(int i)
{
#ifdef MIAMAINWINDOW_H_ENABLE_ZSCORE
    if(useZScore->checkState() == Qt::Checked) {
        for(int ds = 0; ds < datasets.size(); ++ds) {
            Settings s = datasets[ds]->getSettings();

            double newCutoff = i / 10.0;
            std::cout<<"Set cutoff to "<<newCutoff<<" "<<i<<std::endl;

            s.mid_distance_cutoff = newCutoff;
            datasets[ds]->setSettings(s);

            cutOffLabel->setText(QString("d < ") + QString::number(newCutoff));
        }
    } else {
#endif
        double percent = 100.0 * (i - cutOffSlider->minimum())
                / (cutOffSlider->maximum() - cutOffSlider->minimum());
        // set slider maximum to 70%, because only lower range of interest
        percent *= DISTANCE_CUTOFF_SLIDER_IMPLICIT_MAXIMUM;
        networkSet->setRelativeDistanceCutoff(percent);

        cutOffLabel->setText(QString::number(percent) + " %");

#ifdef MIAMAINWINDOW_H_ENABLE_ZSCORE
    }
#endif

    setupExperimentOverlayGraph();
}

void MIAMainWindow::variationSliderChanged(int i)
{
    variationLabel->setText(QString::number(i / 10000.0));
    setupExperimentOverlayGraph();
    updateCompoundList();
}

void MIAMainWindow::distanceMeasureChanged(QString cur)
{
    if(cur == "Euclidean")
        LabelingNetworkSet::distCalc->distanceMeasure = MIDDistanceCalculator::D_EUCLIDEAN;
    else if(cur == "Canberra")
        LabelingNetworkSet::distCalc->distanceMeasure = MIDDistanceCalculator::D_CANBERRA;
    else if(cur == "Cosine")
        LabelingNetworkSet::distCalc->distanceMeasure = MIDDistanceCalculator::D_COSINE;
    else if(cur == "Manhattan")
        LabelingNetworkSet::distCalc->distanceMeasure = MIDDistanceCalculator::D_MANHATTAN;
    else if(cur == "Custom")
        LabelingNetworkSet::distCalc->distanceMeasure = MIDDistanceCalculator::D_CUSTOM;

    // recalculate all distances
    networkSet->createDistanceMatrices();

    // reset distance slider and cutoff
    cutOffSlider->setSliderPosition(0);
    networkSet->setDistanceCutoff(0);

    setupExperimentOverlayGraph();
}

void MIAMainWindow::distanceMeasureNormChanged(QString cur)
{
    if(cur == "MAX")
        LabelingNetworkSet::distCalc->distanceNormalization = MIDDistanceCalculator::DN_MAX;
    else if(cur == "MIN")
        LabelingNetworkSet::distCalc->distanceNormalization = MIDDistanceCalculator::DN_MIN;
    else if(cur == "PROD")
        LabelingNetworkSet::distCalc->distanceNormalization = MIDDistanceCalculator::DN_PROD;
    else if(cur == "SUM")
        LabelingNetworkSet::distCalc->distanceNormalization = MIDDistanceCalculator::DN_SUM;
    else if(cur == "NONE")
        LabelingNetworkSet::distCalc->distanceNormalization = MIDDistanceCalculator::DN_ONE;

    // recalculate all distances
    networkSet->createDistanceMatrices();

    // reset slider pos and cutoff
    cutOffSlider->setSliderPosition(0);
    networkSet->setDistanceCutoff(0);

    setupExperimentOverlayGraph();
}

void MIAMainWindow::selectCompoundLibraryAndIdentify()
{
    // keep current identification if no new match?
    QMessageBox mb;
    mb.setWindowTitle("Identify compounds...");
    mb.setText("Keep current names if no new hit is found?");
    mb.setIcon(QMessageBox::Question);
    mb.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
    mb.setDefaultButton(QMessageBox::No);
    int keep = mb.exec();

    if(keep == QMessageBox::Cancel)
        return;

    // Select library
    QStringList libFiles = QFileDialog::getOpenFileNames(this, "Select compound library",
                                                   QSettings().value("cmp_id_library", QString::fromStdString(CMP_ID_LIBRARY)).toString(),
                                                         "Library files (*.lbr);;All files (*)");

    QStringList::Iterator libIt = libFiles.begin();
    while(libIt != libFiles.end()) {
        QString libFile = *libIt;

        if(libFile.isNull())
            return;

        QSettings().setValue("cmp_id_library", libFile);

        if(keep == QMessageBox::Yes) {
            matchCompoundsAgainstLibrary(libFile, false);
        } else if (QMessageBox::No) {
            matchCompoundsAgainstLibrary(libFile, true);
        }

        ++libIt;
    }
    qDeleteAll(nodeWidgets);
    nodeWidgets.clear(); // need to recreate with new labels TODO: keep ref to compound in the widget?

    setupExperimentOverlayGraph();
    experimentListWidget->updateExperimentList(networkSet->getDatasets(), experimentColors);
    updateCompoundList();
}

void MIAMainWindow::showConfigDialog()
{
    ConfigDialog* conf = new ConfigDialog;

    bool oldUseCommonIon = qsettings.value("nw_use_common_largest_ion", NW_USE_LARGEST_COMMON_ION).toBool();

    //conf->show();
    conf->exec();

    // update excludelib?
    if(excludeLib) delete excludeLib;
    QString libFile = qsettings.value("exclude_library_file", "").toString();
    if(libFile.length())
        excludeLib = gcms::LibrarySearch<int,float>::fromDisk(libFile.toStdString().c_str());

    bool newUseCommonIon = qsettings.value("nw_use_common_largest_ion", NW_USE_LARGEST_COMMON_ION).toBool();
    if(newUseCommonIon != oldUseCommonIon) {
        networkSet->setUseLargestCommonIon(newUseCommonIon);
        qDeleteAll(nodeWidgets);
        nodeWidgets.clear(); // ... and recreate widgets to show newly selected MID
    }

    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
    //matchCompoundsAcrossExperiments();
    //matchCompoundsAgainstLibrary();

    networkSet->setExcludeM0(qsettings.value("alignment_m0", 0).toInt());
    networkSet->createDistanceMatrices();
    setupExperimentOverlayGraph();
    experimentListWidget->updateExperimentList(networkSet->getDatasets(), experimentColors);
    updateCompoundList();

    QApplication::restoreOverrideCursor();
}

#ifdef MIAMAINWINDOW_H_ENABLE_ZSCORE
void MIAMainWindow::useZScoreChanged()
{
    if(useZScore->checkState() == Qt::Checked) {
        cutOffSlider->setMinimum(-30);
        cutOffSlider->setMaximum(0);
        cutOffSlider->setValue(-10);
    } else {
        cutOffSlider->setRange(0, 1000);
        cutOffSlider->setSingleStep(1);
    }

    LabelingNetworkSet::createDistanceMatrices();
}
#endif

void MIAMainWindow::hideLessVaryingChanged()
{
    setupExperimentOverlayGraph();
}

void MIAMainWindow::excludeIfFoundInLessExperimentsChanged(int i)
{
    qsettings.setValue("exclude_if_found_in_less_experiments", i);
    setupExperimentOverlayGraph();
}

void MIAMainWindow::hideFoundInLessExperimentsChanged()
{
    setupExperimentOverlayGraph();
}

void MIAMainWindow::experimentSettingChanged(int idx, Settings s)
{
    NetworkLayer *ds = networkSet->getDataset(idx);
    ds->setSettings(s);
    startLabelDetection(ds);
    // TODO rerun label detection if neccessary -> Settings::settingsChanged
    // No, TODO: delete and recreate ds with the new settings
}

void MIAMainWindow::addExperiment(Settings s)
{
    NetworkLayer *ds = new NetworkLayer(s);
    networkSet->addDataset(ds);
    startLabelDetection(ds);
}

void MIAMainWindow::experimentSelectionChanged()
{
    setupExperimentOverlayGraph();
}

void MIAMainWindow::experimentRemoved(NetworkLayer *ds)
{
    networkSet->removeDataset(ds);
    delete ds;
    labelDetectionFinished();
}

void MIAMainWindow::closeEvent(QCloseEvent *event)
{
    int q = QMessageBox::question(this, "Close MIA?", "Are your sure you want to close MIA?", QMessageBox::Ok, QMessageBox::Cancel);

#if MIA_DEBUG_LEVELg > 1
    // show all widgets
    foreach(QWidget* w, QApplication::topLevelWidgets()) {
        std::cout<<"-"<<w->metaObject()->className()<<std::endl;
    }
#endif

    if(q == QMessageBox::Ok) {
        event->accept();
    } else {
        event->ignore();
    }
}

#ifdef MIA_WITH_NETCDF_IMPORT
void MIAMainWindow::showDataImportDialog()
{
    if(!netCDFImportDialog) {
        netCDFImportDialog = new NetCDFImportDialog();
    }
    netCDFImportDialog->show();
}
#endif


/*
void MIAMainWindow::openFile()
{
   QString filename = QFileDialog::getOpenFileName(this, "Open file", "", "Data files (*.dat);;All files (*)");

   if(filename.isNull())
       return;

   std::cout<<"Loading...";
   QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

   if(!progressDialog) {
       progressDialog = new QProgressDialog(this);
       progressDialog->setModal(true);
       progressDialog->setMinimumDuration(0);
       progressDialog->setMaximum(0);
       connect(this, SIGNAL(progressText(QString)), progressDialog, SLOT(setLabelText(QString)));
       progressDialog->show();
   }

   QFile file(filename);
   file.open(QIODevice::ReadOnly);
   QDataStream in(&file);
   try {
       emit(progressText("Reading file."));
       in >> layers >> nodes;
       qDeleteAll(nodeWidgets);
       nodeWidgets.clear();
   } catch(DeserializationException const &e) {
       file.close();
       QApplication::restoreOverrideCursor();
       QMessageBox mb;
       mb.setWindowTitle("Error opening file...");
       mb.setText("Maybe you used data saved with an older program version?");
       mb.setDetailedText(e.what());
       mb.setIcon(QMessageBox::Critical);
       mb.setStandardButtons(QMessageBox::Ok);
   }

   file.close();

   // refresh everything
   labelDetectionFinished();

   setWindowTitle(QApplication::applicationName() + " - " + filename);

   QApplication::restoreOverrideCursor();
   std::cout<<"done."<<std::endl;
}
*/

void MIAMainWindow::openXMLFile()
{
    QString filename = QFileDialog::getOpenFileName(this, "Open file", qsettings.value("last_open_xml_path", "").toString(),
                                                    "XML files (*.xml);;All files (*)");

    if(!filename.length())
        return;

    std::cout<<"Loading...";
    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

    // TODO allow multiple xml loads
    networkSet->removeAllDatasets();

    try {
        QList<NetworkLayer*> datasets = QList<NetworkLayer*>::fromVector(QVector<NetworkLayer*>::fromStdVector(NetworkLayer::fromXMLFile(filename.toStdString())));

        if(datasets.size() > experimentColors.size()) {
            std::cerr<<"Define more colors"<<std::endl;
            return;
        }

        for(int tracerID = 0; tracerID < datasets.size(); ++tracerID){
            mia::NetworkLayer *ds = datasets[tracerID];
            networkSet->addDataset(ds);
            startLabelDetection(ds);
        }

        setWindowTitle(QApplication::applicationName() + " - " + filename);
        QApplication::restoreOverrideCursor();

    } catch (rapidxml::parse_error const &e) {
        QApplication::restoreOverrideCursor();
        std::cerr<<"Error reading XML file"<<std::endl;

        QMessageBox mb;
        mb.setWindowTitle("Error reading file...");
        mb.setText("Error reading XML file. Check syntax.");
        mb.setIcon(QMessageBox::Critical);
        mb.setStandardButtons(QMessageBox::Ok);
        mb.exec();
    } catch (...) {
        QApplication::restoreOverrideCursor();

        QMessageBox mb;
        mb.setWindowTitle("Error...");
        mb.setText("Error loading data.");
        mb.setIcon(QMessageBox::Critical);
        mb.setStandardButtons(QMessageBox::Ok);
        mb.exec();
    }

    qsettings.setValue("last_open_xml_path", filename);
}

/*
void MIAMainWindow::saveFile()
{
    QString filename = QFileDialog::getSaveFileName(this, "Save file", "", "DAT file (*.dat);;All files (*)");

    if(filename.isNull())
        return;

    std::cout<<"Saving...";
    QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));

    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QDataStream out(&file);

    out << datasets;
    out << nodes;
    //TODO << g;

    file.close();

    QApplication::restoreOverrideCursor();
    std::cout<<"done."<<std::endl;
}
*/

#ifdef MIA_WITH_HDF5
/**
 * @brief Export distance matrices to HDF5
 */
void MIAMainWindow::saveHDF()
{
    QString filename = QFileDialog::getSaveFileName(this, "Save HDF5", "", "HDF5 file (*.hdf *,hdf5);;All files (*)");

    if(filename.isNull())
        return;

    // open
    const H5std_string fileNameH5 = filename.toStdString();
    H5::H5File h5f(fileNameH5, H5F_ACC_TRUNC);


    // create group
    const H5std_string groupName = "Dists";
    H5::Group distsGroup = h5f.createGroup(groupName);

    // write experiment labels
    std::vector<std::string> stringVec;
    for(std::map<std::string, std::vector<std::vector<double> > >::iterator mapIt = distMats.begin();
        mapIt != distMats.end(); ++mapIt) {
        stringVec.push_back(mapIt->first);
    }
    HDFStringWriter::writeVector(distsGroup, "Experiments", stringVec);

    // write metabolite names
    stringVec.clear();
    foreach (NodeCompound* c, nodes.values()) {
        stringVec.push_back(c->getCompoundName());
    }
    HDFStringWriter::writeVector(distsGroup, "Metabolites", stringVec);

    // write distance matrices
    int count = 0;
    for(std::map<std::string, std::vector<std::vector<double> > >::iterator mapIt = distMats.begin();
        mapIt != distMats.end(); ++mapIt) {

        std::vector<std::vector<double> > dist = mapIt->second;

        hsize_t dim[2];
        dim[0] = dist.size();
        dim[1] = dist[0].size();
        H5::DataSpace dataspace(2, dim);
        std::stringstream name;
        name <<"/Dists/Matrix" << count++; // TODO make array, not Matrix1..n
        H5::DataSet dataset = h5f.createDataSet(name.str(), H5::PredType::NATIVE_DOUBLE, dataspace);

        double data[dim[0]][dim[1]];

        for(int i = 0; i < dim[0]; ++i)
            for(int j = 0; j < dim[1]; ++j)
                data[i][j] = dist[i][j];
        dataset.write(data, H5::PredType::NATIVE_DOUBLE);
    }

    h5f.close();
}
#endif

}
