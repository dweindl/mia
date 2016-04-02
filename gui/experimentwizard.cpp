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

#include "experimentwizard.h"
#include "miaguiconstants.h"
#include <iostream>

namespace mia {

ExperimentWizard::ExperimentWizard(QWidget *parent) : QWizard(parent)
{
    addPage(new ExperimentWizardFiles);
    addPage(new ExperimentWizardSettings);
    setWindowTitle(tr("Add tracer..."));
    setDefaultProperty("QDoubleSpinBox", "value", "valueChanged");
}

ExperimentWizard::ExperimentWizard(QString title, QWidget *parent) : QWizard(parent)
{
    addPage(new ExperimentWizardFiles(title));
    addPage(new ExperimentWizardSettings);
    setWindowTitle(tr("Add tracer..."));
    setDefaultProperty("QDoubleSpinBox", "value", "valueChanged");
}

ExperimentWizard::ExperimentWizard(Settings s, QWidget *parent) : QWizard(parent)
{
    addPage(new ExperimentWizardFiles(s));
    addPage(new ExperimentWizardSettings(s));
    setWindowTitle(tr("Edit tracer..."));
}

void ExperimentWizard::accept()
{
    QDialog::accept();
}

Settings ExperimentWizard::getSettings()
{
    Settings s = Settings();
    s.experiment = field("tracer_name").toString().toStdString();
    ExperimentWizardFiles* tw = dynamic_cast<ExperimentWizardFiles*>(page(0));
    foreach (QString i, tw->getLabeledFiles())
        s.labFiles.push_back(i.toStdString());
    foreach (QString i, tw->getUnlabeledFiles())
        s.unlabFiles.push_back(i.toStdString());
    s.cmp_id_ri_tol = field("cmp_id_ri_tol").toDouble();
    s.cmp_id_library = field("cmp_id_library").toString().toStdString();
    // s.cmp_id_mass_filter
    s.cmp_id_score_cutoff = field("cmp_id_score_cutoff").toDouble();
    //s.cmp_id_use_ri
    s.labels_max_hits = field("cmp_id_maxhits").toInt();
    // TODO for different tracer atoms: s.lid_correction_ratio
    s.lid_maximal_frag_dev= field("lid_maximal_frag_dev").toDouble();
    s.lid_min_frag_num = field("lid_min_frag_num").toInt();
    s.lid_req_label_amount = field("lid_req_label_amount").toDouble();
    s.lid_req_r2 = field("lid_req_r2").toDouble();
    s.lid_min_m0 = field("lid_min_m0").toDouble();
    s.lid_max_mass_isotopomer = field("lid_max_mass_isotopomer").toInt();
    s.nw_gap_penalty = field("nw_gap_penalty").toDouble();
    s.mid_distance_cutoff = field("mid_distance_cutoff").toDouble();

    return s;
}

ExperimentWizardFiles::ExperimentWizardFiles(QWidget *parent)
{
    tracerNameText = new QLineEdit();
    registerField("tracer_name*", tracerNameText);
    labFilesList = new QListWidget();
    labFilesList->setSelectionMode(QListWidget::ExtendedSelection);
    unlabFilesList = new QListWidget();
    unlabFilesList->setSelectionMode(QListWidget::ExtendedSelection);
    init();
}

ExperimentWizardFiles::ExperimentWizardFiles(QString title, QWidget *parent)
{
    tracerNameText = new QLineEdit(title);
    registerField("tracer_name*", tracerNameText);
    labFilesList = new QListWidget();
    labFilesList->setSelectionMode(QListWidget::ExtendedSelection);
    unlabFilesList = new QListWidget();
    unlabFilesList->setSelectionMode(QListWidget::ExtendedSelection);
}

ExperimentWizardFiles::ExperimentWizardFiles(Settings s, QWidget *parent)
{
    tracerNameText = new QLineEdit(QString::fromStdString(s.experiment));
    registerField("tracer_name", tracerNameText);
    labFilesList = new QListWidget();

    QStringList files;
    for(std::vector<std::string>::iterator it = s.labFiles.begin(); it != s.labFiles.end(); ++it)
        files.push_back(QString::fromStdString(*it));
    labFilesList->addItems(files);

    unlabFilesList = new QListWidget();
    files.clear();
    for(std::vector<std::string>::iterator it = s.unlabFiles.begin(); it != s.unlabFiles.end(); ++it)
        files.push_back(QString::fromStdString(*it));
    unlabFilesList->addItems(files);

    init();
}

QList<QString> ExperimentWizardFiles::getLabeledFiles()
{
    QList<QString> files;
    for(int i = 0; i < labFilesList->count(); ++i) {
        QListWidgetItem *item = labFilesList->item(i);
        files.append(item->text());
    }
    return files;
}

QList<QString> ExperimentWizardFiles::getUnlabeledFiles()
{
    QList<QString> files;
    for(int i = 0; i < unlabFilesList->count(); ++i) {
        QListWidgetItem *item = unlabFilesList->item(i);
        files.append(item->text());
    }
    return files;
}

void ExperimentWizardFiles::init()
{
    setTitle(tr("Select files"));

    tracerNameLabel = new QLabel(tr("Tracer name:"));
    labFilesLabel = new QLabel(tr("Labeled chromatograms:"));
    unlabFilesLabel = new QLabel(tr("Unlabeled chromatograms:"));

    addLabFilesButton = new QPushButton(QIcon(":/gui/icons/list-add.png"), tr("&Add"));
    removeLabFilesButton = new QPushButton(QIcon(":/gui/icons/list-remove.png"), tr("&Remove"));
    addUnlabFilesButton = new QPushButton(QIcon(":/gui/icons/list-add.png"), tr("&Add"));
    removeUnlabFilesButton = new QPushButton(QIcon(":/gui/icons/list-remove.png"), tr("&Remove"));

    registerField("lab_files", labFilesList);
    registerField("unlab_files", unlabFilesList);

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(tracerNameLabel);
    layout->addWidget(tracerNameText);

    layout->addWidget(labFilesLabel);
    layout->addWidget(labFilesList);
    layout->addWidget(addLabFilesButton);
    layout->addWidget(removeLabFilesButton);

    layout->addWidget(unlabFilesLabel);
    layout->addWidget(unlabFilesList);
    layout->addWidget(addUnlabFilesButton);
    layout->addWidget(removeUnlabFilesButton);

    setLayout(layout);

    connect(addLabFilesButton, SIGNAL(clicked()), this, SLOT(addLabeledChromatogram()));
    connect(removeLabFilesButton, SIGNAL(clicked()), this, SLOT(removeLabeledChromatogram()));
    connect(addUnlabFilesButton, SIGNAL(clicked()), this, SLOT(addUnlabeledChromatogram()));
    connect(removeUnlabFilesButton, SIGNAL(clicked()), this, SLOT(removeUnlabeledChromatogram()));
}

void ExperimentWizardFiles::addLabeledChromatogram()
{
    QSettings settings;
    QStringList files = QFileDialog::getOpenFileNames(this, "Select labeled chromatograms", settings.value("dir_add_chromatogram", "").toString(), "MetaboliteDetector files (*.cmp)");

    if(files.size()) {
        labFilesList->addItems(files);
        QDir dir(files[0]);
        settings.setValue("dir_add_chromatogram", dir.absolutePath());
    }
}

void ExperimentWizardFiles::addUnlabeledChromatogram()
{
    QSettings settings;
    QStringList files = QFileDialog::getOpenFileNames(this, "Select labeled chromatograms", settings.value("dir_add_chromatograms", "").toString(), "MetaboliteDetector files (*.cmp)");

    if(files.size()) {
        unlabFilesList->addItems(files);
        QDir dir(files[0]);
        settings.setValue("dir_add_chromatogram", dir.absolutePath());
    }
}

void ExperimentWizardFiles::removeLabeledChromatogram()
{
    QList<QListWidgetItem*> items = labFilesList->selectedItems();
    foreach(QListWidgetItem* item, items) {
        delete item;
    }
    labFilesList->update();
}

void ExperimentWizardFiles::removeUnlabeledChromatogram()
{
    QList<QListWidgetItem*> items = unlabFilesList->selectedItems();
    foreach(QListWidgetItem* item, items) {
        delete item;
    }
    unlabFilesList->update();
}

ExperimentWizardSettings::ExperimentWizardSettings(QWidget *parent)
{

    QSettings settings;

    cmp_id_ri_tol_spinbox = new QSpinBox();
    cmp_id_ri_tol_spinbox->setValue(settings.value("cmp_id_ri_tol", CMP_ID_RI_TOL).toInt());
    cmp_id_score_cutoff_dblspinbox = new QDoubleSpinBox();
    cmp_id_score_cutoff_dblspinbox->setValue(settings.value("cmp_id_score_cutoff", CMP_ID_SCORE_CUTOFF).toDouble());
    cmp_id_maxhits_spinbox = new QSpinBox();
    cmp_id_maxhits_spinbox->setValue(settings.value("cmp_id_maxhits", LABELS_MAX_HITS).toInt());
//    cmp_id_library_text = new QLineEdit(settings.value("cmp_id_library", QString::fromStdString(CMP_ID_LIBRARY)).toString());

    lid_maximal_frag_dev_dblspinbox = new QDoubleSpinBox();
    lid_maximal_frag_dev_dblspinbox->setValue(settings.value("lid_maximal_frag_dev", LID_MAXIMAL_FRAG_DEV).toDouble());
    lid_min_frag_num_spinbox = new QSpinBox();
    lid_min_frag_num_spinbox->setValue(settings.value("lid_min_frag_num", LID_MIN_FRAG_NUM).toInt());
    lid_req_label_amount_dblspinbox = new QDoubleSpinBox();
    lid_req_label_amount_dblspinbox->setValue(settings.value("lid_req_label_amount", LID_REQ_LABEL_AMOUNT).toDouble());
    lid_req_r2_dblspinbox = new QDoubleSpinBox();
    lid_req_r2_dblspinbox->setValue(settings.value("lid_req_r2", LID_REQ_R2).toDouble());
    lid_max_mass_isotopomer_spinbox = new QSpinBox();
    lid_max_mass_isotopomer_spinbox->setValue(settings.value("lid_max_mass_isotopomer", LID_MAX_MASS_ISOTOPOMER).toInt());
    lid_min_m0_dblspinbox = new QDoubleSpinBox();
    lid_min_m0_dblspinbox->setValue(settings.value("lid_min_m0", LID_MIN_M0).toDouble());
    nw_gap_penalty_dblspinbox = new QDoubleSpinBox();
    nw_gap_penalty_dblspinbox->setValue(settings.value("nw_gap_penalty", NW_GAP_PENALTY).toDouble());
    mid_distance_cutoff_dblspinbox = new QDoubleSpinBox();
    mid_distance_cutoff_dblspinbox->setValue(settings.value("mid_distance_cutoff", MID_DISTANCE_CUTOFF).toDouble());

    init();
}

ExperimentWizardSettings::ExperimentWizardSettings(Settings s, QWidget *parent)
{
    cmp_id_ri_tol_spinbox = new QSpinBox();
    cmp_id_ri_tol_spinbox->setValue(s.cmp_id_ri_tol);
    cmp_id_score_cutoff_dblspinbox = new QDoubleSpinBox();
    cmp_id_score_cutoff_dblspinbox->setValue(s.cmp_id_score_cutoff);
    cmp_id_maxhits_spinbox = new QSpinBox();
    cmp_id_maxhits_spinbox->setValue(s.labels_max_hits);
//    cmp_id_library_text = new QLineEdit(QString::fromStdString(s.cmp_id_library));

    lid_maximal_frag_dev_dblspinbox = new QDoubleSpinBox();
    lid_maximal_frag_dev_dblspinbox->setValue(s.lid_maximal_frag_dev);
    lid_min_frag_num_spinbox = new QSpinBox();
    lid_min_frag_num_spinbox->setValue(s.lid_min_frag_num);
    lid_req_label_amount_dblspinbox = new QDoubleSpinBox();
    lid_req_label_amount_dblspinbox->setValue(s.lid_req_label_amount);
    lid_req_r2_dblspinbox = new QDoubleSpinBox();
    lid_req_r2_dblspinbox->setValue(s.lid_req_r2);
    lid_max_mass_isotopomer_spinbox = new QSpinBox();
    lid_max_mass_isotopomer_spinbox->setValue(s.lid_max_mass_isotopomer);
    lid_min_m0_dblspinbox = new QDoubleSpinBox();
    lid_min_m0_dblspinbox->setValue(s.lid_min_m0);
    nw_gap_penalty_dblspinbox = new QDoubleSpinBox();
    nw_gap_penalty_dblspinbox->setValue(s.nw_gap_penalty);
    mid_distance_cutoff_dblspinbox = new QDoubleSpinBox();
    mid_distance_cutoff_dblspinbox->setValue(s.mid_distance_cutoff);

    init();
}

void ExperimentWizardSettings::init()
{
    setTitle(tr("Settings..."));

    cmp_id_ri_tol_spinbox->setToolTip(tr("Retention index tolerance for library matching."));
    cmp_id_maxhits_spinbox->setToolTip(tr("Maximum number of library hits to show (-1 = show all)"));

    cmp_id_ri_tol_label = new QLabel(tr("RI &Tolerance:"));
    cmp_id_ri_tol_label->setBuddy(cmp_id_ri_tol_spinbox);
    cmp_id_score_cutoff_label = new QLabel(tr("Compound identification cutoff score:"));
    cmp_id_maxhits_label = new QLabel(tr("Show top n hits:"));
//    cmp_id_library_label = new QLabel(tr("Library for compound identification:"));
    lid_maximal_frag_dev_label = new QLabel(tr("Maximum fragment deviation:"));
    lid_min_frag_num_label = new QLabel(tr("Minimum number of labeled fragments:"));
    lid_req_label_amount_label = new QLabel(tr("Required amount of isotopic enrichment:"));
    lid_req_r2_label = new QLabel(tr("Minimum R^2:"));
    lid_max_mass_isotopomer_label = new QLabel(tr("Ignore compounds with M_n | with n > ...:"));
    lid_min_m0_label = new QLabel(tr("Minimum M0 abd.:"));
    nw_gap_penalty_label = new QLabel(tr("Gap penalty for Needleman-Wunsch-Scoring:"));
    mid_distance_cutoff_label = new QLabel(tr("Distance cutoff for graph edges:"));

    cmp_id_maxhits_spinbox->setMinimum(-1);
    cmp_id_maxhits_spinbox->setMaximum(5);
    cmp_id_score_cutoff_dblspinbox->setMinimum(0);
    cmp_id_score_cutoff_dblspinbox->setMaximum(1);
    cmp_id_score_cutoff_dblspinbox->setSingleStep(0.01);
    cmp_id_ri_tol_spinbox->setMinimum(0);
    cmp_id_ri_tol_spinbox->setMaximum(1E6);
    cmp_id_ri_tol_spinbox->setSingleStep(10);
    lid_maximal_frag_dev_dblspinbox->setMinimum(0);
    lid_maximal_frag_dev_dblspinbox->setMaximum(1);
    lid_maximal_frag_dev_dblspinbox->setSingleStep(0.01);
    lid_min_frag_num_spinbox->setMinimum(1);
    lid_min_frag_num_spinbox->setMaximum(10);
    lid_req_label_amount_dblspinbox->setMinimum(0);
    lid_req_label_amount_dblspinbox->setMaximum(1);
    lid_req_label_amount_dblspinbox->setSingleStep(0.01);
    lid_req_r2_dblspinbox->setMinimum(0);
    lid_req_r2_dblspinbox->setMaximum(1);
    lid_req_r2_dblspinbox->setSingleStep(0.01);
    lid_max_mass_isotopomer_spinbox->setMinimum(1);
    lid_max_mass_isotopomer_spinbox->setMaximum(100);
    lid_min_m0_dblspinbox->setMinimum(0);
    lid_min_m0_dblspinbox->setMaximum(1);
    lid_min_m0_dblspinbox->setSingleStep(0.01);
    nw_gap_penalty_dblspinbox->setMinimum(0);
    nw_gap_penalty_dblspinbox->setMaximum(1);
    nw_gap_penalty_dblspinbox->setSingleStep(0.01);
    mid_distance_cutoff_dblspinbox->setMinimum(0);
    mid_distance_cutoff_dblspinbox->setMaximum(1);
    mid_distance_cutoff_dblspinbox->setSingleStep(0.01);

    registerField("cmp_id_ri_tol", cmp_id_ri_tol_spinbox);
    registerField("cmp_id_score_cutoff", cmp_id_score_cutoff_dblspinbox, "value");
    registerField("cmp_id_maxhits", cmp_id_maxhits_spinbox);
//    registerField("cmp_id_library", cmp_id_library_text);
    registerField("lid_maximal_frag_dev", lid_maximal_frag_dev_dblspinbox, "value");
    registerField("lid_min_frag_num", lid_min_frag_num_spinbox);
    registerField("lid_req_label_amount", lid_req_label_amount_dblspinbox, "value");
    registerField("lid_req_r2", lid_req_r2_dblspinbox, "value");
    registerField("lid_min_m0", lid_min_m0_dblspinbox, "value");
    registerField("lid_max_mass_isotopomer", lid_max_mass_isotopomer_spinbox);
    registerField("nw_gap_penalty", nw_gap_penalty_dblspinbox, "value");
    registerField("mid_distance_cutoff", mid_distance_cutoff_dblspinbox, "value");

    QGridLayout *layout = new QGridLayout;

    QGridLayout *cmpIDLayout = new QGridLayout;
    cmpIDLayout->addWidget(cmp_id_ri_tol_label, 0, 0);
    cmpIDLayout->addWidget(cmp_id_ri_tol_spinbox, 0, 1);
    cmpIDLayout->addWidget(cmp_id_score_cutoff_label, 1, 0);
    cmpIDLayout->addWidget(cmp_id_score_cutoff_dblspinbox, 1, 1);
    cmpIDLayout->addWidget(cmp_id_maxhits_label, 2, 0);
    cmpIDLayout->addWidget(cmp_id_maxhits_spinbox, 2, 1);
//    cmpIDLayout->addWidget(cmp_id_library_label);
//    cmpIDLayout->addWidget(cmp_id_library_text);
    QGroupBox *cmpIDGroupBox = new QGroupBox("Compound identification", this);
    cmpIDGroupBox->setLayout(cmpIDLayout);
    layout->addWidget(cmpIDGroupBox, 0, 0, 1, 2);

    QGridLayout *lidLayout = new QGridLayout;
    int row = 0;
    lidLayout->addWidget(lid_maximal_frag_dev_label, row, 0);
    lidLayout->addWidget(lid_maximal_frag_dev_dblspinbox, row++, 1);
    lidLayout->addWidget(lid_min_frag_num_label, row, 0);
    lidLayout->addWidget(lid_min_frag_num_spinbox, row++, 1);
    lidLayout->addWidget(lid_req_label_amount_label, row, 0);
    lidLayout->addWidget(lid_req_label_amount_dblspinbox, row++, 1);
    lidLayout->addWidget(lid_req_r2_label, row, 0);
    lidLayout->addWidget(lid_req_r2_dblspinbox, row++, 1);
    lidLayout->addWidget(lid_min_m0_label, row, 0);
    lidLayout->addWidget(lid_min_m0_dblspinbox, row++, 1);
    lidLayout->addWidget(lid_max_mass_isotopomer_label, row, 0);
    lidLayout->addWidget(lid_max_mass_isotopomer_spinbox, row++, 1);

    QGroupBox *lidGroupBox = new QGroupBox("Label detection", this);
    lidGroupBox->setLayout(lidLayout);
    layout->addWidget(lidGroupBox, 1, 0, 1, 2);

    layout->addWidget(nw_gap_penalty_label, 2, 0);
    layout->addWidget(nw_gap_penalty_dblspinbox, 2, 1);
    layout->addWidget(mid_distance_cutoff_label, 3, 0);
    layout->addWidget(mid_distance_cutoff_dblspinbox, 3, 1);
    setLayout(layout);
}

/*
void ExperimentWizardSettings::selectLibrary()
{
    QSettings settings;
    QString file = QFileDialog::getOpenFileName(this, "Select compound library", settings.value("dir_compound_library", "").toString(), "MetaboliteDetector libraries (*.lbr)");
    cmp_id_library_text->setText(file);
    QDir dir(file);
    settings.setValue("dir_compound_library", dir.absolutePath());
}
*/
}
