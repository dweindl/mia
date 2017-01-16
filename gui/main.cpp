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

#include <string>

#include <QApplication>
#include <QDataStream>

#include"miamainwindow.h"
#include"graphvizqt.h"
#include"miaguiconstants.h"

#ifdef QT_STATIC
#ifndef _WINDOWS_
    Q_IMPORT_PLUGIN(QXcbIntegrationPlugin)
#endif
#ifdef _WINDOWS_
    Q_IMPORT_PLUGIN(QWindowsIntegrationPlugin)
#endif
#endif

int main(int argc, char *argv[])
{
    //Q_INIT_RESOURCE(application);
    QApplication app(argc, argv);
    app.setApplicationVersion(MIA_VERSION);
    std::cout<<"MIA version "<<app.applicationVersion().toStdString()<<std::endl;

    app.setOrganizationName("LCSB");
    app.setOrganizationDomain("lcsb.uni.lu");
    app.setApplicationName("MIA - Mass Isotopolome Analyzer");
    app.setWindowIcon(QIcon(":/gui/icons/programmIcon16x16.png"));

    // Splashscreen
    QPixmap pixmap(":/gui/splash.png");
    QSplashScreen splash(pixmap);
    splash.show();
    app.processEvents();

    mia::MIAMainWindow mw;
    mw.resize(mw.maximumSize());
    mw.showMaximized();

    splash.close();

    return app.exec();
}
