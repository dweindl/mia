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

#ifndef CONFIGDIALOG_H
#define CONFIGDIALOG_H

#include <QSettings>
#include <QDialog>
#include <QtWidgets>
#include <QStackedWidget>

namespace mia {

class ConfigurationPagePaths;
class ConfigurationPageMisc;

//class ConfigurationPageCompoundIdentification;
//class ConfigurationPageColors;

class ConfigDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ConfigDialog(QWidget *parent = 0);
    
signals:
    
public slots:
    void changePage(QListWidgetItem *current, QListWidgetItem *previous);
    void saveConfig();

private:
    QListWidget *contentsWidget; // TOC list
    QStackedWidget *pagesWidget; // individual config pages

    ConfigurationPagePaths* pagePaths;
    ConfigurationPageMisc* pageMisc;

    void createTOCIcons();

};

class ConfigurationPagePaths : public QWidget
{
    Q_OBJECT

public:
    ConfigurationPagePaths(QWidget *parent = 0);
    void save();

public slots:
    void selectXMLPath();

private:

    QLineEdit *textXMLPath;
    QLineEdit *textLibPath;
    QLineEdit *textDataPath;
    QLineEdit *textExcludeLib;
    QLineEdit *textChromatogramPath;
};


class ConfigurationPageMisc : public QWidget
{
    Q_OBJECT

public:
    ConfigurationPageMisc(QWidget *parent = 0);
    void save();

public slots:

private:
    QButtonGroup *nodeLayoutGroup;
    QRadioButton *nodeLayoutVertical;
    QRadioButton *nodeLayoutHorizontal;
    QRadioButton *nodeLayoutSquare;
    QGroupBox *nodeLayoutGroupBox;

    QButtonGroup *m0OptionsGroup;
    QRadioButton *m0Include;
    QRadioButton *m0Omit;
    QRadioButton *m0OmitScaleBasepeak;
    QRadioButton *m0OmitScaleSum;
    QGroupBox *m0OptionsGroupBox;

    QCheckBox *showUnconnectedNodes;
    QCheckBox *useLargestCommonIon;

    void initNodeLayoutGroupBox();
    void initM0OptionsGroupBox();

    QSettings s;

};


}

#endif // CONFIGDIALOG_H
