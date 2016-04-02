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

#ifndef EXPERIMENTLISTWIDGET_H
#define EXPERIMENTLISTWIDGET_H

#include <QWidget>
#include <QSignalMapper>
#include "experimentwizard.h"
#include "../src/settings.h"
#include "../src/networklayer.h"

Q_DECLARE_METATYPE(mia::NetworkLayer*) // for use with UserRole

namespace mia {

class ExperimentListWidget : public QWidget
{
    Q_OBJECT
public:
    explicit ExperimentListWidget(QWidget *parent = 0);
    void updateExperimentList(QList<NetworkLayer*> layers, const QList<QColor> &experimentColors);
    QList<NetworkLayer*> getVisibleExperiments(QList<NetworkLayer*> datasets);

signals:
    void experimentSettingChanged(int idx, Settings s);
    void experimentAdded(Settings s);
    void experimentUsageChanged();
    void experimentRemoved(NetworkLayer *layer);

public slots:
    void useAllExperimentsChanged(int state);
    void experimentClicked(int row, int col);
    void experimentDoubleClicked(int row, int col);
    void showAddExperimentWizard();
    void experimentWizardFinished(int idx = -1);
    void experimentSelectionStateChanged(int idx);
    void removeExperiment();
    void saveXML();

private:
    QTableWidget *experimentTable;
    QPushButton *experimentAddButton;
    QPushButton *experimentRemoveButton;
    QPushButton *saveXMLButton;
    QSignalMapper *signalMapperEditSettings;
    ExperimentWizard *tw;
    QSignalMapper *experimentSelectionMapper;

    static const int nameColIdx = 1;
};

}
#endif // EXPERIMENTLISTWIDGET_H
