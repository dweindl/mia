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

#ifndef MIDPLOT_H
#define MIDPLOT_H

#include <QWidget>
#include <QColor>

/**
 * @brief The MIDPlot class creates a bar plot for the given MID. Adapted from ntfd-gui.
 */
class MIDPlot : public QWidget
{
      Q_OBJECT
public:
      MIDPlot(const QVector<double>& mids, const QString& title, double r2, int ion, QWidget* parent);

      virtual void paintEvent ( QPaintEvent * event );
      virtual QSize sizeHint () const;

      void setBarColor(QColor c);
      void setCI(const QVector<double>& ci);

      QColor getBgColor() const;
      void setBgColor(const QColor &value);

      QColor getBorderColor() const;
      void setBorderColor(const QColor &value);

      void mouseDoubleClickEvent(QMouseEvent *event);

public slots:
      void copyToClipboard();

private:
      int margin;
      QVector<double> mids;
      QVector<double> cis;
      QString title;
      QString subtitle;
      double sum;
      double r2;
      int ion;
      int bar_width;
      int space;
      int height;
      QColor barColor;
      QColor bgColor;
      QColor borderColor;
      bool displayValues;
};

#endif // MIDPLOT_H
