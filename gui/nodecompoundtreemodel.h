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

#ifndef NODECOMPOUNDTREEMODEL_H
#define NODECOMPOUNDTREEMODEL_H

#include <QAbstractItemModel>
#include <QStyledItemDelegate>
#include <QStringList>
#include <QPainter>

#include <labeledcompound.h>

#include "src/nodecompound.h"

namespace mia {

class NodeCompoundTreeItem;
class NodeCompoundTreeModel;


class NodeCompoundTreeModel : public QAbstractItemModel
{
    Q_OBJECT
public:
    explicit NodeCompoundTreeModel(QObject *parent = 0);
    explicit NodeCompoundTreeModel(std::vector<NodeCompound *>, QObject *parent = 0);
    explicit NodeCompoundTreeModel(QMap<int, NodeCompound *>, QObject *parent = 0);

    ~NodeCompoundTreeModel();

    QVariant data(const QModelIndex &index, int role) const;
    bool setData(const QModelIndex &index, const QVariant &value,
                 int role = Qt::EditRole);
    Qt::ItemFlags flags(const QModelIndex &index) const;
    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const;
    QModelIndex index(int row, int column,
                      const QModelIndex &parent = QModelIndex()) const;
    QModelIndex parent(const QModelIndex &index) const;

    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    void setChildrenData(const QModelIndex &index, const QVariant &value, int role);

signals:
    

private:
    NodeCompoundTreeItem *rootItem;
    static const int CHECKBOX_USE_COLUMN = 999; // column number of the checkbox // TODO = 0 to enable checkboxes
};

class NodeCompoundTreeItem {
public:
    //explicit NodeCompoundTreeItem(const QList<QVariant> &data, NodeCompoundTreeItem *parent = 0);
    //explicit NodeCompoundTreeItem(QList<NodeCompound *nc>, NodeCompoundTreeItem *parent = 0); //
    explicit NodeCompoundTreeItem(); // for root
    explicit NodeCompoundTreeItem(NodeCompound *nc, NodeCompoundTreeItem *parent = 0); // compound level
    explicit NodeCompoundTreeItem(std::string tracer, labid::LabeledCompound *lc, NodeCompoundTreeItem *parent = 0); // compound->tracer level
    explicit NodeCompoundTreeItem(labid::LabeledCompound *lc, int _ionIdx, NodeCompoundTreeItem *parent = 0); // compound->tracer->ion level
    explicit NodeCompoundTreeItem(int n, double abundance, double ci, NodeCompoundTreeItem *parent = 0); // compound->tracer->ion level

    ~NodeCompoundTreeItem();

    void appendChild(NodeCompoundTreeItem *child);

    NodeCompoundTreeItem *child(int row);
    int childCount() const;
    int columnCount() const;
    QVariant data(int column, int role = Qt::DisplayRole) const;
    bool setData(int column, const QVariant &value, int role);
    int row() const;
    NodeCompoundTreeItem *parent();

private:
    enum ItemRole {
        Compound,
        LabeledCompound,
        LabeledIon,
        MassIsotopomer
    };

    ItemRole role;  // or better subclass for different roles?
    QList<NodeCompoundTreeItem*> childItems;
    QList<QVariant> itemData;
    NodeCompoundTreeItem *parentItem;

    NodeCompound* nc;
    labid::LabeledCompound *lc;
    int ionIdx;
    QVariant userRoleData;

    QStringList header;
    static const int CHECKBOX_USE_COLUMN = 999; // column number of the checkbox // TODO = 0 to enable checkboxes
};

}
#endif // NODECOMPOUNDTREEMODEL_H
