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

#include "experimentlistwidget.h"
#include <QtGui>

namespace mia {

ExperimentListWidget::ExperimentListWidget(QWidget *parent) :
    QWidget(parent)
{
    experimentAddButton = new QPushButton(QIcon(":/gui/icons/list-add.png"), tr("&Add"));
    connect(experimentAddButton, SIGNAL(clicked()), this, SLOT(showAddExperimentWizard()));

    experimentRemoveButton = new QPushButton(QIcon(":/gui/icons/list-remove.png"), tr("&Remove"));
    connect(experimentRemoveButton, SIGNAL(clicked()), this, SLOT(removeExperiment()));

    saveXMLButton = new QPushButton(QIcon(":/gui/icons/document-save.png"), tr("&Save XML"));
    connect(saveXMLButton, SIGNAL(clicked()), this, SLOT(saveXML()));

    QCheckBox *useAll = new QCheckBox("Use all/none for similarity analysis");
    useAll->setCheckState(Qt::Checked);
    connect(useAll, SIGNAL(stateChanged(int)), this, SLOT(useAllExperimentsChanged(int)));

    experimentTable = new QTableWidget();
    connect(experimentTable, SIGNAL(cellClicked(int,int)), this, SLOT(experimentClicked(int, int)));
    connect(experimentTable, SIGNAL(cellDoubleClicked(int,int)), this, SLOT(experimentDoubleClicked(int, int)));

    QGridLayout *experimentGrid = new QGridLayout();
    experimentGrid->addWidget(experimentAddButton);
    experimentGrid->addWidget(experimentRemoveButton);
    experimentGrid->addWidget(saveXMLButton);
    experimentGrid->addWidget(useAll);
    experimentGrid->addWidget(experimentTable);

    setLayout(experimentGrid);

    experimentSelectionMapper = new QSignalMapper(this);
    connect(experimentSelectionMapper, SIGNAL(mapped(int)), this, SLOT(experimentSelectionStateChanged(int)));
}

void ExperimentListWidget::experimentWizardFinished(int idx)
{
    if(tw->result() == QDialog::Rejected)
        return;

    Settings s = tw->getSettings(); // TODO add library fileselectiondialog button

    if(idx >= 0) {
        // edit dialog, settings were changes
        emit(experimentSettingChanged(idx, s));
    } else {
        // add new experiment
        if(tw->result() == QDialog::Rejected)
            return;

        emit(experimentAdded(s));
    }

    delete tw;
}

void ExperimentListWidget::useAllExperimentsChanged(int state)
{
    for(int row = 0; row < experimentTable->rowCount(); ++row) {
        QCheckBox *chk = static_cast<QCheckBox *>(experimentTable->cellWidget(row, 0));
        QVariant v = experimentTable->item(row, nameColIdx)->data(Qt::UserRole);
        NetworkLayer *layer = v.value<NetworkLayer *>();
        layer->setVisible(state == Qt::Checked);
        chk->setCheckState((Qt::CheckState)state);
    }
    emit(experimentUsageChanged());
}

void ExperimentListWidget::showAddExperimentWizard()
{
    tw = new ExperimentWizard(this); // QString("Dataset %1").arg(experimentTable->rowCount() + 1),
    tw->show();
    connect(tw, SIGNAL(finished(int)), this, SLOT(experimentWizardFinished()));
}

void ExperimentListWidget::experimentDoubleClicked(int row, int col)
{
    QVariant v = experimentTable->item(row, nameColIdx)->data(Qt::UserRole);
    NetworkLayer *layer = v.value<NetworkLayer *>();
    tw = new ExperimentWizard(layer->getSettings(), this);
    tw->show();
    signalMapperEditSettings = new QSignalMapper(this);
    connect(tw, SIGNAL(finished(int)), signalMapperEditSettings, SLOT(map()));
    signalMapperEditSettings->setMapping(tw, row);
    connect(signalMapperEditSettings, SIGNAL(mapped(int)), this, SLOT(experimentWizardFinished(int)));
}

void ExperimentListWidget::experimentClicked(int row, int col)
{
    std::cerr<<row<<col<<std::endl;
}

void ExperimentListWidget::updateExperimentList(QList<NetworkLayer*> layers, const QList<QColor> &experimentColors)
{
    experimentTable->setRowCount(layers.size());
    experimentTable->setColumnCount(4);
    experimentTable->setHorizontalHeaderLabels(QStringList()<<"Use"<<"Dataset"<<"Lab"<<"Unlab");
    QTableWidgetItem *item;
    int i = 0;
    foreach(NetworkLayer* layer, layers) {
        item = new QTableWidgetItem(QString::fromStdString(layer->getSettings().experiment));
        item->setTextColor(experimentColors[i % experimentColors.size()]); // TODO wrong color sequence. add colors to layer
        item->setData(Qt::UserRole, QVariant::fromValue(layer));

        int column = 0;

        QCheckBox *chk = new QCheckBox(experimentTable);
        chk->setChecked(layer->isVisible());
        connect(chk, SIGNAL(clicked()), experimentSelectionMapper, SLOT(map()));
        experimentSelectionMapper->setMapping(chk, i);
        experimentTable->setCellWidget(i, column++, chk);

        experimentTable->setItem(i, column++, item);
        experimentTable->setItem(i, column++, new QTableWidgetItem(QString::number(layer->cmpLab.size())));
        experimentTable->setItem(i, column++, new QTableWidgetItem(QString::number(layer->cmpUnlab.size())));
        ++i;
    }
    experimentTable->resizeColumnsToContents();
}

/**
 * @brief An experiment was selected or deselected
 * @param idx index in "layers"
 */
void ExperimentListWidget::experimentSelectionStateChanged(int idx)
{
    QCheckBox *chk = static_cast<QCheckBox *>(experimentTable->cellWidget(idx, 0));
    QVariant v = experimentTable->item(idx, nameColIdx)->data(Qt::UserRole);
    NetworkLayer *layer = v.value<NetworkLayer *>();
    layer->setVisible(chk->isChecked());

    emit(experimentUsageChanged());
}

void ExperimentListWidget::removeExperiment()
{
    QList<QTableWidgetItem*> items = experimentTable->selectedItems();
    if(!items.size())
        return;

    foreach(QTableWidgetItem* item, items) {
        QVariant v = experimentTable->item(item->row(), nameColIdx)->data(Qt::UserRole);
        NetworkLayer *layer = v.value<NetworkLayer *>();
        emit(experimentRemoved(layer));
    }
}

void ExperimentListWidget::saveXML()
{
    // write current datasets to XML file for easy reload
    QString filename = QFileDialog::getSaveFileName(this, "Save current dataset as XML", "", "XML files (*.xml);;All files (*)");

    if(filename.isNull())
        return;

    QFile f(filename);
    f.open(QIODevice::WriteOnly);
    QTextStream out(&f);

    out<<"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    out<<"<experiments>\n";

    for(int row = 0; row < experimentTable->rowCount(); ++row) {
        QVariant v = experimentTable->item(row, nameColIdx)->data(Qt::UserRole);
        NetworkLayer *layer = v.value<NetworkLayer *>();
        out<<layer->getSettings().toXML().c_str()<<"\n";
    }
    out<<"</experiments>\n";
}

}
