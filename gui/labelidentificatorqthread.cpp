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

#include "labelidentificatorqthread.h"
#include "QApplication"
#include "QMessageBox"
#include "gcmsscanexception.h"

namespace mia {

LabelIdentificatorQThread::LabelIdentificatorQThread(LabelingDataset *ds, QObject *parent) :
    QThread(parent), ds(ds)
{
}

void LabelIdentificatorQThread::setProgress(size_t value, size_t max)
{
    emit(progressMax(max));
    emit(progress(value));
}

void LabelIdentificatorQThread::setMessage(std::string message)
{
    emit(progressMessage(message.c_str()));
}

void LabelIdentificatorQThread::run()
{
    try {
        ds->findLabeledCompounds(this);
    } catch (gcms::GCMSScanException e) {
        QApplication::restoreOverrideCursor();

        QMessageBox mb;
        mb.setWindowTitle("Error...");
        mb.setText(QString::fromStdString("Error loading data:\n\n%1").arg(QString::fromStdString(e.getMessage())));
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
}

}
