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

#include "configdialog.h"

#include <QtWidgets>
#include "nodewidget.h"
#include "miaguiconstants.h"

namespace mia {

ConfigDialog::ConfigDialog(QWidget *parent) :
    QDialog(parent)
{
    contentsWidget = new QListWidget;
    contentsWidget->setViewMode(QListView::IconMode);
    contentsWidget->setIconSize(QSize(48, 48));
    contentsWidget->setMaximumWidth(80);
    contentsWidget->setSpacing(2);

    pagesWidget = new QStackedWidget;

    pagePaths = new ConfigurationPagePaths;
    pagesWidget->addWidget(pagePaths);
    pageMisc = new ConfigurationPageMisc;
    pagesWidget->addWidget(pageMisc);

    QPushButton *closeButton = new QPushButton(tr("Close"));
    // TODO QPushButton *defaultButton = new QPushButton(tr("Defaults"));

    createTOCIcons();
    contentsWidget->setCurrentRow(0);

    connect(closeButton, SIGNAL(clicked()), this, SLOT(accept()));
    // TODO connect(defaultButton, SIGNAL(clicked()), this, SLOT(close()));

    QHBoxLayout *horizontalLayout = new QHBoxLayout;
    horizontalLayout->addWidget(contentsWidget);
    horizontalLayout->addWidget(pagesWidget, 1);

    QHBoxLayout *buttonsLayout = new QHBoxLayout;
    buttonsLayout->addStretch(1);
    buttonsLayout->addWidget(closeButton);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addLayout(horizontalLayout);
    mainLayout->addStretch(1);
    mainLayout->addSpacing(12);
    mainLayout->addLayout(buttonsLayout);
    setLayout(mainLayout);

    setWindowTitle(QApplication::applicationName() + tr(" Config"));

    connect(this, SIGNAL(accepted()), this, SLOT(saveConfig()));
}

void ConfigDialog::changePage(QListWidgetItem *current, QListWidgetItem *previous)
{
    if (!current)
        current = previous;

    pagesWidget->setCurrentIndex(contentsWidget->row(current));
}

void ConfigDialog::saveConfig()
{
    pagePaths->save();
    pageMisc->save();
}

void ConfigDialog::createTOCIcons()
{
    QListWidgetItem *pathsButton = new QListWidgetItem(contentsWidget);
    pathsButton->setIcon(QIcon::fromTheme("folder-open")); // emblem-system
    pathsButton->setText(tr("Paths"));
    pathsButton->setTextAlignment(Qt::AlignHCenter);
    pathsButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);

    connect(contentsWidget,
            SIGNAL(currentItemChanged(QListWidgetItem*, QListWidgetItem*)),
            this, SLOT(changePage(QListWidgetItem*, QListWidgetItem*)));

    QListWidgetItem *miscButton = new QListWidgetItem(contentsWidget);
    miscButton->setIcon(QIcon::fromTheme("configure"));
    miscButton->setText(tr("Misc"));
    miscButton->setTextAlignment(Qt::AlignHCenter);
    miscButton->setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);

    connect(contentsWidget,
            SIGNAL(currentItemChanged(QListWidgetItem*, QListWidgetItem*)),
            this, SLOT(changePage(QListWidgetItem*, QListWidgetItem*)));

}

ConfigurationPagePaths::ConfigurationPagePaths(QWidget *parent)
{
    QSettings s;

    QVBoxLayout *mainLayout = new QVBoxLayout;

    // XML
    QLabel *labelXMLPath = new QLabel("XML directory:");
    textXMLPath = new QLineEdit(s.value("open_xml_path", QDir::homePath()).toString());
    QPushButton *buttonXMLPath = new QPushButton(QIcon::fromTheme("folder-open"), "Select");
    connect(buttonXMLPath, SIGNAL(clicked()), this, SLOT(selectXMLPath()));

    mainLayout->addWidget(labelXMLPath);
    mainLayout->addWidget(textXMLPath);
    mainLayout->addWidget(buttonXMLPath);

    // Lib
    QLabel *labelLibPath = new QLabel("Library directory:");
    textLibPath = new QLineEdit(s.value("open_library_path", QDir::homePath()).toString());
    mainLayout->addWidget(labelLibPath);
    mainLayout->addWidget(textLibPath);

    // Data
    QLabel *labelDataPath = new QLabel("Binary data directory:");
    textDataPath = new QLineEdit(s.value("open_binary_path", QDir::homePath()).toString());
    mainLayout->addWidget(labelDataPath);
    mainLayout->addWidget(textDataPath);

    QLabel *labelChromatogramPath = new QLabel("Chromatogram directory:");
    textChromatogramPath = new QLineEdit(s.value("chromatogram_path", QDir::homePath()).toString());
    mainLayout->addWidget(labelChromatogramPath);
    mainLayout->addWidget(textChromatogramPath);

    // ExcludeLib
    QLabel *labelExcludeLib = new QLabel("Exclude library:");
    textExcludeLib = new QLineEdit(s.value("exclude_library_file", "").toString());
    mainLayout->addWidget(labelExcludeLib);
    mainLayout->addWidget(textExcludeLib);

    setLayout(mainLayout);
}

void ConfigurationPagePaths::save()
{
    QSettings s;
    s.setValue("open_xml_path", textXMLPath->text());
    s.setValue("open_library_path", textLibPath->text());
    s.setValue("open_binary_path", textDataPath->text());
    s.setValue("chromatogram_path", textChromatogramPath->text());
    s.setValue("exclude_library_file", textExcludeLib->text());
}

void ConfigurationPagePaths::selectXMLPath()
{
    QString dir = QFileDialog::getExistingDirectory(this, "Select directory...", textXMLPath->text());

    if(dir.isNull())
        return;

    textXMLPath->setText(dir);
}

ConfigurationPageMisc::ConfigurationPageMisc(QWidget *parent)
{
    initNodeLayoutGroupBox();
    initM0OptionsGroupBox();

    QVBoxLayout *vl = new QVBoxLayout(this);

    showUnconnectedNodes = new QCheckBox("Show unconnected nodes", this);
    showUnconnectedNodes->setChecked(s.value("graph_show_unconnected_nodes", true).toBool());
    vl->addWidget(showUnconnectedNodes);

    useLargestCommonIon = new QCheckBox("Use largest common ion", this);
    useLargestCommonIon->setChecked(s.value("nw_use_common_largest_ion", NW_USE_LARGEST_COMMON_ION).toBool());
    vl->addWidget(useLargestCommonIon);

    vl->addWidget(nodeLayoutGroupBox);

    vl->addWidget(m0OptionsGroupBox);

    setLayout(vl);
}

void ConfigurationPageMisc::save()
{
    s.setValue("nodewidget_layout_direction", nodeLayoutGroup->checkedId());

    s.setValue("alignment_m0", m0OptionsGroup->checkedId());

    s.setValue("graph_show_unconnected_nodes", showUnconnectedNodes->isChecked());

    s.setValue("nw_use_common_largest_ion", useLargestCommonIon->isChecked());
}

void ConfigurationPageMisc::initNodeLayoutGroupBox()
{
    int nodeLayout = s.value("nodewidget_layout_direction", NodeWidget::LAYOUT_VERTICAL).toInt();

    nodeLayoutGroup = new QButtonGroup(this);

    nodeLayoutVertical = new QRadioButton("Vertical", this);
    nodeLayoutHorizontal = new QRadioButton("Horizontal", this);
    nodeLayoutSquare = new QRadioButton("Square", this);

    switch (nodeLayout) {
    case NodeWidget::LAYOUT_VERTICAL:
        nodeLayoutVertical->setChecked(true);
        break;
    case NodeWidget::LAYOUT_HORIZONTAL:
        nodeLayoutHorizontal->setChecked(true);
        break;
    case NodeWidget::LAYOUT_SQUARE:
        nodeLayoutSquare->setChecked(true);
        break;
    }

    nodeLayoutGroup->addButton(nodeLayoutVertical, NodeWidget::LAYOUT_VERTICAL);
    nodeLayoutGroup->addButton(nodeLayoutHorizontal, NodeWidget::LAYOUT_HORIZONTAL);
    nodeLayoutGroup->addButton(nodeLayoutSquare, NodeWidget::LAYOUT_SQUARE);

    nodeLayoutGroupBox = new QGroupBox("MID plot layout", this);
    QVBoxLayout *vl = new QVBoxLayout(nodeLayoutGroupBox);
    vl->addWidget(nodeLayoutVertical);
    vl->addWidget(nodeLayoutHorizontal);
    vl->addWidget(nodeLayoutSquare);
    nodeLayoutGroupBox->setLayout(vl);
}

void ConfigurationPageMisc::initM0OptionsGroupBox()
{
    int m0 = s.value("alignment_m0", 0).toInt();

    m0OptionsGroup = new QButtonGroup(this);

    m0Include = new QRadioButton("Include", this);
    m0Omit = new QRadioButton("Exclude M0", this);
    m0OmitScaleBasepeak = new QRadioButton("Exclude M0, normalize to basepeak", this);
    m0OmitScaleSum = new QRadioButton("Exclude M0, normalize to sum", this);

    switch (m0) {
    case 0:
        m0Include->setChecked(true);
        break;
    case 1:
        m0Omit->setChecked(true);
        break;
    case 2:
        m0OmitScaleBasepeak->setChecked(true);
        break;
    case 3:
        m0OmitScaleSum->setChecked(true);
        break;
    default:
        m0Include->setChecked(true);
    }

    m0OptionsGroup->addButton(m0Include, 0);
    m0OptionsGroup->addButton(m0Omit, 1);
    m0OptionsGroup->addButton(m0OmitScaleBasepeak, 2);
    m0OptionsGroup->addButton(m0OmitScaleSum, 3);

    m0OptionsGroupBox = new QGroupBox("Distance calculation M0 options", this);
    QVBoxLayout *vl = new QVBoxLayout(m0OptionsGroupBox);
    vl->addWidget(m0Include);
    vl->addWidget(m0Omit);
    vl->addWidget(m0OmitScaleBasepeak);
    vl->addWidget(m0OmitScaleSum);

    m0OptionsGroupBox->setLayout(vl);
}


}
