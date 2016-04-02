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

#ifndef LABELIDENTIFICATORQTHREAD_H
#define LABELIDENTIFICATORQTHREAD_H

#include <QThread>
#include "labelidentificator.h"
#include "src/labelingdataset.h"

namespace mia {

class LabelIdentificatorQThread : public QThread, public labid::LabelIdentificatorProgressListener
{
    Q_OBJECT
public:
    explicit LabelIdentificatorQThread(LabelingDataset *ds, QObject *parent = 0);
    void setProgress(size_t value, size_t max);
    void setMessage(std::string message);
    void run();

signals:
    void progress(int value);
    void progressMessage(QString message);
    void progressMax(int max);
    void errorMessage(QString message);

public slots:

private:
    LabelingDataset *ds;
};


}

#endif // LABELIDENTIFICATORQTHREAD_H
