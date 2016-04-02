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

#ifndef ExperimentWizard_H
#define ExperimentWizard_H

#include <QWizard>
#include <QtWidgets>
#include "src/settings.h"

namespace mia {

class ExperimentWizard;
class ExperimentWizardFiles;
class ExperimentWizardSettings;

/**
 * @brief The ExperimentWizard class shows a wizard to add a new label-experiment.
 */
class ExperimentWizard : public QWizard
{
    Q_OBJECT
public:
    explicit ExperimentWizard(QWidget *parent = 0);
    explicit ExperimentWizard(QString title, QWidget *parent = 0);
    explicit ExperimentWizard(Settings s, QWidget *parent = 0);

    void accept();
    Settings getSettings();
};

/**
 * @brief The file-selection wizard page.
 */
class ExperimentWizardFiles : public QWizardPage {
    Q_OBJECT

public:
    ExperimentWizardFiles(QWidget *parent = 0);
    ExperimentWizardFiles(QString title, QWidget *parent = 0);
    ExperimentWizardFiles(Settings s, QWidget *parent = 0);
    QList<QString> getLabeledFiles();
    QList<QString> getUnlabeledFiles();

private:
    void init();
    QLabel *tracerNameLabel;
    QLabel *labFilesLabel;
    QLabel *unlabFilesLabel;
    QLineEdit *tracerNameText;
    QListWidget *labFilesList;
    QListWidget *unlabFilesList;
    QPushButton *addLabFilesButton;
    QPushButton *addUnlabFilesButton;
    QPushButton *removeLabFilesButton;
    QPushButton *removeUnlabFilesButton;

public slots:
    void addLabeledChromatogram();
    void addUnlabeledChromatogram();
    void removeLabeledChromatogram();
    void removeUnlabeledChromatogram();
};

/**
 * @brief The settings wizard page.
 */
class ExperimentWizardSettings: public QWizardPage {
    Q_OBJECT

public:
    ExperimentWizardSettings(QWidget *parent = 0);
    ExperimentWizardSettings(Settings s, QWidget *parent = 0);

public slots:
//    void selectLibrary();

private:
    void init();

    QLabel *cmp_id_ri_tol_label;
    QSpinBox *cmp_id_ri_tol_spinbox;
    QLabel *cmp_id_score_cutoff_label;
    QDoubleSpinBox *cmp_id_score_cutoff_dblspinbox;
    QLabel *cmp_id_maxhits_label;
    QSpinBox *cmp_id_maxhits_spinbox;
//    QLabel *cmp_id_library_label;
//    QLineEdit *cmp_id_library_text;
    QLabel *lid_maximal_frag_dev_label;
    QDoubleSpinBox *lid_maximal_frag_dev_dblspinbox;
    QLabel *lid_min_frag_num_label;
    QSpinBox *lid_min_frag_num_spinbox;
    QLabel *lid_req_label_amount_label;
    QDoubleSpinBox *lid_req_label_amount_dblspinbox;
    QLabel *lid_req_r2_label;
    QDoubleSpinBox *lid_req_r2_dblspinbox;
    QLabel *lid_max_mass_isotopomer_label;
    QSpinBox *lid_max_mass_isotopomer_spinbox;
    QLabel *lid_min_m0_label;
    QDoubleSpinBox *lid_min_m0_dblspinbox;
    QLabel *nw_gap_penalty_label;
    QDoubleSpinBox *nw_gap_penalty_dblspinbox;
    QLabel *mid_distance_cutoff_label;
    QDoubleSpinBox *mid_distance_cutoff_dblspinbox;

    /* TODO
    bool lid_filter_by_conf_interval; // NTFD: hardcoded
    int lid_min_signal_to_noise; // NTFD: hardcoded
    double lid_required_spec_freq;  // NTFD: hardcoded

    double lid_sensitivity;
    double lid_maximal_frag_dev;
    double lid_correction_ratio;
    //bool nw_exclude_m0; /** Gap penalty for needleman wunsch scoring */

};
}
#endif // ExperimentWizard_H
