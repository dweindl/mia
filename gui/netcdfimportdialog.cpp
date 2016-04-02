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

#include "netcdfimportdialog.h"
#include "compounddetector.h"
#include "defaultpeakdetector.h"
#include "cdfimporter.h"

namespace mia {

NetCDFImportDialog::NetCDFImportDialog(QWidget *parent) : QWidget(parent)
{
    setWindowTitle("MIA - Data import");

    listFiles = new QListWidget();
    listFiles->setSelectionMode(QListWidget::ExtendedSelection);

    QPushButton *btAddFiles = new QPushButton(QIcon(":/gui/icons/list-add.png"), tr("&Add"), this);
    QPushButton *btRemoveFiles = new QPushButton(QIcon(":/gui/icons/list-remove.png"), tr("&Remove"), this);

    connect(btAddFiles, SIGNAL(clicked()), this, SLOT(addFileClicked()));
    connect(btRemoveFiles, SIGNAL(clicked()), this, SLOT(removeFileClicked()));

    QLabel *lbDeconvolutionWidth = new QLabel("Deconvolution width", this);
    lbDeconvolutionWidth->setToolTip("Deconvolution width in scans");

    sbDeconvolutionWidth = new QDoubleSpinBox(this);
    sbDeconvolutionWidth->setMinimum(0);
    sbDeconvolutionWidth->setValue(5);
    sbDeconvolutionWidth->setMaximum(100);

    QLabel *lbPeakThreshold = new QLabel("Peak threshold", this);
    lbPeakThreshold->setToolTip("Peak threshold");

    sbPeakThreshold = new QDoubleSpinBox(this);
    sbPeakThreshold->setMinimum(0);
    sbPeakThreshold->setValue(10);
    sbPeakThreshold->setMaximum(100);

    QLabel *lbMinPeakHeight = new QLabel("Minimum peak height", this);
    lbPeakThreshold->setToolTip("Minimum peak height");

    sbMinPeakHeight = new QDoubleSpinBox(this);
    sbMinPeakHeight->setMinimum(0);
    sbMinPeakHeight->setValue(10);
    sbMinPeakHeight->setMaximum(100);

    QLabel *lbMinPeaks = new QLabel("Minimum number of peaks", this);
    lbMinPeaks->setToolTip("Minimum number of peaks");

    sbMinPeaks = new QSpinBox(this);
    sbMinPeaks->setMinimum(1);
    sbMinPeaks->setValue(25);
    sbMinPeaks->setMaximum(1000);

    QPushButton *btAccept = new QPushButton("Start", this);
    QPushButton *btCancel = new QPushButton("Cancel", this);

    connect(btAccept, SIGNAL(clicked()), this, SLOT(okayClicked()));
    connect(btCancel, SIGNAL(clicked()), this, SLOT(cancelClicked()));


    QGridLayout *layout = new QGridLayout(this);
    int row = 0;
    layout->addWidget(lbDeconvolutionWidth, row, 0);
    layout->addWidget(sbDeconvolutionWidth, row, 1);
    ++row;
    layout->addWidget(lbPeakThreshold, row, 0);
    layout->addWidget(sbPeakThreshold, row, 1);
    ++row;
    layout->addWidget(lbMinPeakHeight, row, 0);
    layout->addWidget(sbMinPeakHeight, row, 1);
    ++row;
    layout->addWidget(lbMinPeaks, row, 0);
    layout->addWidget(sbMinPeaks, row, 1);
    ++row;
    layout->addWidget(listFiles, row, 0, 2, 1);
    layout->addWidget(btAddFiles, row, 1);
    ++row;
    layout->addWidget(btRemoveFiles, row, 1);
    layout->setRowStretch(row, 1);

    layout->addWidget(btAccept, 0, 2);
    layout->addWidget(btCancel, 1, 2);
}

void NetCDFImportDialog::addFileClicked()
{
    QSettings settings;
    QStringList files = QFileDialog::getOpenFileNames(this,
                                                      "Select chromatograms",
                                                      settings.value("dir_add_chromatogram", "").toString(),
                                                      "MetaboliteDetector files (*.cmp);;netCDF files (*.cdf);;All files (*.*)"
                                                      );

    if(files.size()) {
        listFiles->addItems(files);
        QDir dir(files[0]);
        settings.setValue("dir_add_chromatogram", dir.absolutePath());
    }
}

void NetCDFImportDialog::removeFileClicked()
{
    QList<QListWidgetItem*> items = listFiles->selectedItems();
    foreach(QListWidgetItem* item, items) {
        delete item;
    }
    listFiles->update();
}

void NetCDFImportDialog::okayClicked()
{
    QProgressDialog progress("Compound detection...\n(This might take a few minutes)",
                             "Abort",
                             0, 2 * listFiles->count(),
                             this);
    progress.setWindowModality(Qt::WindowModal);
    progress.show();
    QCoreApplication::processEvents();

    // peak detection / deconvolution settings
    gcms::GCMSSettings::AN_DECONVOLUTION_WIDTH = sbDeconvolutionWidth->value();
    gcms::GCMSSettings::AN_MIN_PEAK_HEIGHT = sbMinPeakHeight->value();
    gcms::GCMSSettings::AN_PEAK_THRESHOLD_BEGIN = sbPeakThreshold->value();
    gcms::GCMSSettings::AN_PEAK_THRESHOLD_END = -1.0 * gcms::GCMSSettings::AN_PEAK_THRESHOLD_BEGIN;


    gcms::GCMSSettings::AN_MIN_PEAK_NUMBER = sbMinPeaks->value();

    for(int i = 0; i < listFiles->count(); ++i) {

        QListWidgetItem *item = listFiles->item(i);
        // files.append(item->text());
        QString file = item->text();
        QString fileBase = file.left(file.length() - 4); // strip file ending
        QString ext = file.right(4);

        std::cout<<file.toStdString()<<"...\n";
        progress.setValue(2 * i);
        QCoreApplication::processEvents();

        if(ext.compare(".cdf", Qt::CaseInsensitive) == 0) {
            // import netCDF

            std::cout<<"\tnetCDF import..."<<file.toStdString()<<"\n";

            try
            {
                  gcms::CdfImporter imp(file.toStdString().c_str());
                  imp.importData(fileBase.toStdString().c_str());
            }
            catch ( gcms::GCMSScanException e )
            {
                  // emit ( errorMessage( tr ( e.getMessage().c_str() ) ) );
            }
            catch ( ... )
            {
                  // emit ( errorMessage( ( tr ( "Error during NetCDF import" ) ) ));
            }

            ext = ".bin";
        }

        progress.setValue(2 * i + 1);
        QCoreApplication::processEvents();

        if (progress.wasCanceled())
            break;

        if(ext.compare(".bin", Qt::CaseInsensitive) == 0 || ext.compare(".cmp", Qt::CaseInsensitive) == 0) {
            // detect compounds
            std::cout<<"\tcompound detection..."<<fileBase.toStdString()<<"\n";
            gcms::DefaultPeakDetector<int,float>* detector=new gcms::DefaultPeakDetector<int,float>();
            CompoundDetector cd(fileBase, detector, true);
            cd.run();
        }

        std::cout<<"done\n";

        if (progress.wasCanceled())
            break;
    }

    if(progress.wasCanceled()) {
        QMessageBox::information(this, "MIA", "Compound detection finished.", QMessageBox::Ok);
    } else {
        QMessageBox::critical(this, "MIA", "Compound detection canceled.", QMessageBox::Ok);
    }

    close();
}

void NetCDFImportDialog::cancelClicked()
{
    close();
}

}
