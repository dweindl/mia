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

#include "midplot.h"

#include <iostream>
#include <cassert>
#include <cmath>

#include<QPainter>
#include<QApplication>
#include<QClipboard>

MIDPlot::MIDPlot(const QVector<double>& mids, const QString& title, double r2, int ion, QWidget* parent)
    :QWidget(parent),
      mids(mids),
      title(title),
      r2(r2),
      ion(ion),
      bar_width(10),
      space(5),
      height(220)
{
    bgColor = Qt::lightGray; // this->palette().window().color();
    borderColor = Qt::black;
    barColor = Qt::blue;
    margin = 10;

    // bar width calculation
    QFont f = font();
    f.setPointSize(12);
    QFontMetrics fm = QFontMetrics(f);
    bar_width = qMax(bar_width, fm.width("M10"));

    setToolTip(title);
    displayValues = 1;

    sum = 0;
    for(int i = 0; i < mids.size(); ++i)
        sum += fabs(mids[i]);

    subtitle = "Ion: " + QString::number(ion)
            + " R2 = " + QString::number(r2, 'f', 3)
            + " S = " + QString::number(sum, 'f', 2);
}

void MIDPlot::paintEvent ( QPaintEvent * event )
{

    // background box
    QPainter p(this);
    p.fillRect(rect(), bgColor);

    QPen oldPen = p.pen();
    QPen newPen = p.pen();
    newPen.setColor(borderColor);
    p.setPen(newPen);
    p.drawRect(rect());
    p.setPen(oldPen);


    QFont f(font());

    //paint title
    f.setPointSize(16);
    p.setFont(f);
    QRect rect_title1(margin, margin, 1.01 * (width() - 2 * margin), QFontMetrics(f).height());
    p.drawText(rect_title1, Qt::AlignCenter, title);

    // subtitle
    f.setPointSize(14);
    p.setFont(f);
    QRect rect_title2 = rect_title1;
    rect_title2.moveTop(rect_title2.bottom());
    p.drawText(rect_title2, Qt::AlignCenter, subtitle);
    p.translate(0, QFontMetrics(f).height() * 2.5);

    f.setPointSize(12);
    p.setFont(f);

    int n_mids = mids.size();
    int bw = bar_width;
    int bh = height - margin; //border
    int fh = QFontMetrics(f).height();
    bh -= fh;

    int x = margin + qMax(0.0, 0.5 * (QFontMetrics(f).width(title) - n_mids * bw - (n_mids-1)*space));
    //int y = height-5;
    for(int i = 0; i < n_mids; i++)
    {
        //paint bar
        double rc_mid = mids.at(i);
        rc_mid = qMin(1.0, rc_mid);
        rc_mid = qMax(0.0, rc_mid);
        int rc_bh = bh * rc_mid + 0.5;
        p.fillRect(QRect(x, height - rc_bh - fh, bw, rc_bh), barColor);
        p.drawRect(QRect(x, height - rc_bh - fh, bw, rc_bh));

        if(cis.size() && cis[i] > 0) { // -1 if NaN or similar
            int ciLen = bh * cis[i] + 0.5;
            QPen oldPen = p.pen();
            QPen pen = p.pen();
            pen.setWidth(2);
            p.setPen(pen);
            p.drawLine(x + 0.5 * bw, height - rc_bh - fh + ciLen, x + 0.5 * bw, height - rc_bh - fh - ciLen);
            p.setPen(oldPen);
        }

        if(displayValues) {
            if(int mid = mids.at(i) * 100) {
                QRect rectVal(x, height - rc_bh - 2*fh, bw, fh);
                QFont vf;
                vf.setPointSize(10);
                p.setFont(vf);
                p.drawText(rectVal, Qt::AlignHCenter, QString::number(mid));
                p.setFont(f);
            }
        }

        // "M_n"
        QRect rect_text(x, height - margin, bw, fh);
        p.drawText(rect_text, Qt::AlignHCenter, "M" + QString::number(i));
        x += bw+space;
    }
}

QSize MIDPlot::sizeHint () const
{
    int width = bar_width * mids.size() + space*(mids.size() - 1) + 2* margin;

    QFont f;
    f.setPointSize(16);
    int text_width1 = QFontMetrics(f).width(title) + space * 2;
    f.setPointSize(14);
    int text_width2 = QFontMetrics(f).width(subtitle) + margin * 2;
    int rc_width = qMax(width, text_width1);
    rc_width = qMax(rc_width, text_width2);
    int rc_height = height + QFontMetrics(f).height() * 3 + 5;
    return QSize(rc_width, rc_height);
}

void MIDPlot::setBarColor(QColor c)
{
    barColor = c;
}

void MIDPlot::setCI(const QVector<double> &ci)
{
    assert(ci.size() == mids.size());
    cis = ci;
}
QColor MIDPlot::getBgColor() const
{
    return bgColor;
}

void MIDPlot::setBgColor(const QColor &value)
{
    bgColor = value;
}
QColor MIDPlot::getBorderColor() const
{
    return borderColor;
}

void MIDPlot::setBorderColor(const QColor &value)
{
    borderColor = value;
}

void MIDPlot::mouseDoubleClickEvent(QMouseEvent *event)
{
    copyToClipboard();
}

void MIDPlot::copyToClipboard()
{
    QString colSep = "\t";
    QString rowSep = "\n";

    QClipboard *clipboard = QApplication::clipboard();

    QStringList sl;
    sl<<title;
    for(int i = 0; i < mids.size(); ++i) {
        sl<<QString::number(mids[i]) + colSep + QString::number(cis[i]);
    }

    clipboard->setText(sl.join(rowSep));
}



