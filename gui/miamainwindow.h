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

#ifndef MIAMAINWINDOW_H
#define MIAMAINWINDOW_H

#include <QtGui>

#include "rapidxml/rapidxml.hpp"

#include "compound.h"

#include "src/config.h"
#include "src/networklayer.h"
#include "src/nodecompound.h"
#include "labelingnetworkset.h"

#include "nwview.h"
#include "nodewidget.h"
#include "graphvizqt.h"
#include "labelidentificatorqthread.h"
#include "experimentlistwidget.h"
#include "graphvizmia.h"

// Let user choose z-score normalization for distance cutoff
// needs to be checked
#undef MIAMAINWINDOW_H_ENABLE_ZSCORE

#ifdef MIA_WITH_NETCDF_IMPORT
    #include "netcdfimportdialog.h"
#endif

#ifdef MIA_WITH_METABOBASE
    #include "kegg_reaction_test/keggreactionmapper.h"
#endif

namespace mia {

class MIAMainWindow;

class MIAMainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MIAMainWindow(QWidget *parent = 0);
    ~MIAMainWindow();

    static void findLabeledCompounds(NetworkLayer *ds);

signals:
    // TODO use!
    void graphChanged();
    void progressText(QString str);
    
public slots:
    void recreateGraph();
//    void repaintGraph();
    void labelDetectionFinished();
    void compoundClicked(QModelIndex mi);
    //void saveFile();
#ifdef MIA_WITH_HDF5
    void saveHDF();
#endif
    void exportMIDs();
    //void openFile();
    void openXMLFile();
    void updateCompoundList();
    void showInfoDialog();
    void cutOffSliderChanged(int i);
    void variationSliderChanged(int i);
    void distanceMeasureChanged(QString cur);
    void distanceMeasureNormChanged(QString cur);
    void selectCompoundLibraryAndIdentify();
    void showConfigDialog();
#ifdef MIAMAINWINDOW_H_ENABLE_ZSCORE
    void useZScoreChanged();
#endif
    void hideLessVaryingChanged();
    void excludeIfFoundInLessExperimentsChanged(int i);
    void hideFoundInLessExperimentsChanged();
    void experimentSettingChanged(int idx, Settings s);
    void addExperiment(Settings s);
    void experimentSelectionChanged();
    void experimentRemoved(NetworkLayer *ds);
    void closeEvent(QCloseEvent *event);
#ifdef MIA_WITH_NETCDF_IMPORT
    void showDataImportDialog();
#endif

private:
    QList<QColor> experimentColors; /** Colors for the different experiment layers */
    QSettings qsettings; /** Application user settings */

    NWView *view;
    QMap<int, NodeWidget*> nodeWidgets; /** All the different compounds found in any experiment, index is the ID-feature of the compound */
    graphvizqt::Graph *g; /** The network graph */
    int graphSizeWarningLimit;
    QList<LabelIdentificatorQThread*> labidThreads;
#ifdef MIA_WITH_METABOBASE
    KEGGReactionMapper *keggMapper;
#endif

#ifdef MIA_WITH_NETCDF_IMPORT
    NetCDFImportDialog *netCDFImportDialog;
#endif
    gcms::LibrarySearch<int,float> *excludeLib;

    // GUI stuff
    QDockWidget *compoundListDockWidget;
    QDockWidget *experimentListDockWidget;
    QDockWidget *graphOptionsDockWidget;
//    QMenu *fileMenu;
//    QToolBar *fileToolBar;
    QGraphicsScene *scene;
    QTreeView *compoundTreeView;
    QSlider* cutOffSlider;
    QLabel* cutOffLabel;
#ifdef MIAMAINWINDOW_H_ENABLE_ZSCORE
    QCheckBox* useZScore;
#endif
    QSlider* variationSlider;
    QLabel* variationLabel;
    QCheckBox* hideLessVaryingNodes;
    QSpinBox *excludeIfFoundInLessExperiments;
    QCheckBox* hideFoundInLessExperiments;
    QProgressDialog* progressDialog;
    ExperimentListWidget *experimentListWidget;
    LabelingNetworkSet *networkSet;

    void init();
    void matchCompoundsAcrossExperiments();
    void startLabelDetection(NetworkLayer *ds); // TODO move to labelingNetworkSet, needs signals & slots first
    void generateNodeWidgets();
    void matchCompoundsAgainstLibrary(QString libFile = "", bool overwriteNames = true);

    void setupGUI();
    void setupToolbar();
    void setupCompoundList(); /** GUI stuff */
    void setupExperimentList(); /** GUI stuff */
    void setupExperimentOverlayGraph(); /** Do multi-tracer overlay */
    void setupGraphOptionPanel();
    void addEdgesToGraph(int excludeIfFoundInLessExperiments, double variationCutoff);

    bool showGraphSizeWarning(int edges);
};
}
#endif // MIAMAINWINDOW_H
