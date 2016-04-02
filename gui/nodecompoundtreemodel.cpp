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

#include "nodecompoundtreemodel.h"
#include "src/config.h"

#include <QCheckBox>
#include <QRadioButton>
#include <QApplication>

namespace mia {

NodeCompoundTreeModel::NodeCompoundTreeModel(QObject *parent) :
    QAbstractItemModel(parent)
{
}

NodeCompoundTreeModel::NodeCompoundTreeModel(std::vector<NodeCompound *> ncs, QObject *parent) : QAbstractItemModel(parent)
{
    //TODO
    rootItem = new NodeCompoundTreeItem;

    // setup model data
    for(std::vector<NodeCompound*>::iterator it = ncs.begin(); it != ncs.end(); ++it) {
        // create top-level items i.e. abstract compounds
        rootItem->appendChild(new NodeCompoundTreeItem(*it, rootItem));
    }
}

NodeCompoundTreeModel::NodeCompoundTreeModel(QMap<int, NodeCompound *> ncs, QObject *parent)
{
    //TODO
    rootItem = new NodeCompoundTreeItem;

    // setup model data
    foreach(NodeCompound *nc, ncs) {
        // create top-level items i.e. abstract compounds
        rootItem->appendChild(new NodeCompoundTreeItem(nc, rootItem));
    }
}

NodeCompoundTreeModel::~NodeCompoundTreeModel()
{
    delete rootItem;
}

QVariant NodeCompoundTreeModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    NodeCompoundTreeItem *item = static_cast<NodeCompoundTreeItem*>(index.internalPointer());
    return item->data(index.column(), role);
}

bool NodeCompoundTreeModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if(index.isValid()){
        switch (role) {
        case Qt::DisplayRole:
        case Qt::EditRole:
        case Qt::CheckStateRole:
            if(index.column() == CHECKBOX_USE_COLUMN) {
                // set checkmarks hierachically
                // TODO: set parents grey
                setChildrenData(index, value, role);
            }
            emit dataChanged(index, index);
            return static_cast<NodeCompoundTreeItem*>(index.internalPointer())->setData(index.column(), value, role);
            break;
        }
    }
    return false;
}

void NodeCompoundTreeModel::setChildrenData(const QModelIndex &index, const QVariant &value, int role)
{
    //if(index.column() == 0) {
    NodeCompoundTreeItem* item = static_cast<NodeCompoundTreeItem*>(index.internalPointer());
    for(int i = 0; i < item->childCount(); ++i) {
        setData(this->index(i, index.column(), index), value, role);
        setChildrenData(this->index(i, index.column(), index), value, role);
    }
    emit dataChanged(this->index(0, index.column(), index), this->index(item->childCount() - 1, index.column(), index));
}

Qt::ItemFlags NodeCompoundTreeModel::flags(const QModelIndex &index) const
{
    if (!index.isValid())
        return 0;

    if(index.column() == CHECKBOX_USE_COLUMN) // checkboxes in first column | move to item?
        return QAbstractItemModel::flags(index) | Qt::ItemIsUserCheckable;

    return QAbstractItemModel::flags(index) | Qt::ItemIsEditable; // Allows copy/paste
}

QVariant NodeCompoundTreeModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (orientation == Qt::Horizontal && role == Qt::DisplayRole)
        return rootItem->data(section);

    return QVariant();
}

int NodeCompoundTreeModel::rowCount(const QModelIndex &parent) const
{
    NodeCompoundTreeItem *parentItem;
    if (parent.column() > 0)
        return 0;

    if (!parent.isValid())
        parentItem = rootItem;
    else
        parentItem = static_cast<NodeCompoundTreeItem*>(parent.internalPointer());

    return parentItem->childCount();
}

int NodeCompoundTreeModel::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return static_cast<NodeCompoundTreeItem*>(parent.internalPointer())->columnCount();
    else
        return rootItem->columnCount();
}

QModelIndex NodeCompoundTreeModel::index(int row, int column, const QModelIndex &parent) const
{
    if (!hasIndex(row, column, parent))
        return QModelIndex();

    NodeCompoundTreeItem *parentItem;

    if (!parent.isValid())
        parentItem = rootItem;
    else
        parentItem = static_cast<NodeCompoundTreeItem*>(parent.internalPointer());

    NodeCompoundTreeItem *childItem = parentItem->child(row);
    if (childItem)
        return createIndex(row, column, childItem);
    else
        return QModelIndex();
}

QModelIndex NodeCompoundTreeModel::parent(const QModelIndex &index) const
{
    if (!index.isValid())
        return QModelIndex();

    NodeCompoundTreeItem *childItem = static_cast<NodeCompoundTreeItem*>(index.internalPointer());
    NodeCompoundTreeItem *parentItem = childItem->parent();

    if (parentItem == rootItem)
        return QModelIndex();

    return createIndex(parentItem->row(), 0, parentItem);
}

NodeCompoundTreeItem::NodeCompoundTreeItem()
{
    // Use root item for header data
    parentItem = 0;
    header = QStringList();
    header /*<< "Use"*/ << "Compound / Experiment"<< "Abd." << "CI";

    foreach(QString s, header){
        itemData << s;
    }
}

NodeCompoundTreeItem::NodeCompoundTreeItem(NodeCompound *_nc, NodeCompoundTreeItem *parent)
{
    // Abstract compound level
    role = Compound;
    parentItem = parent;
    nc = _nc;
    userRoleData = QString::fromStdString(nc->getFeature(COMPOUND_GROUPING_FEATURE));

    // create children
    std::vector<std::string> tracers = _nc->getExperiments();
    for(std::vector<std::string>::iterator it = tracers.begin(); it != tracers.end(); ++it) {
        labid::LabeledCompound *lc =_nc->getLabeledCompound(*it);
        if(lc > 0)
            appendChild(new NodeCompoundTreeItem(*it, lc, this));
    }

    // add data
    // itemData.push_back(Qt::Checked);
    itemData.push_back(QString::fromStdString(_nc->getCompoundName()));
    itemData.push_back(QString::number(childCount()).append(" experiments"));
}

NodeCompoundTreeItem::NodeCompoundTreeItem(std::string tracer, labid::LabeledCompound *_lc, NodeCompoundTreeItem *parent)
{
    // The tracer experiment
    role = LabeledCompound;
    parentItem = parent;
    lc = _lc;

    // create children
    int numIons = _lc->getLabeledIons().size();
    for(int i = 0; i < numIons; ++i) {
        appendChild(new NodeCompoundTreeItem(lc, i, this));
    }

    // add data
    // itemData.push_back(Qt::Checked);
    itemData.push_back(QString::fromStdString(tracer));
    itemData.push_back(QString::number(childCount()).append(" ions"));
}

NodeCompoundTreeItem::NodeCompoundTreeItem(labid::LabeledCompound *_lc, int _ionIdx, NodeCompoundTreeItem *parent)
{
    // Labeled ion from one compound
    role = LabeledIon;
    parentItem = parent;
    lc = _lc;
    ionIdx = _ionIdx;

    // create children
    const std::vector< double > abd = lc->getIsotopomers().at(ionIdx);
    const std::vector< double > ci = lc->getConfidenceIntervals().at(ionIdx);
    for(int i = 0; i < abd.size(); ++i) {
        appendChild(new NodeCompoundTreeItem(i, abd[i], ci[i], this));
    }

    // add data
    // itemData.push_back(Qt::Checked);
    itemData.push_back(QString("m/z %1").arg(lc->getLabeledIons().at(ionIdx)));
    itemData.push_back(QVariant(lc->getR2s().at(ionIdx)));
    itemData.push_back(QString("R2 = ").append(QString::number(childCount())));
}

NodeCompoundTreeItem::NodeCompoundTreeItem(int n, double abundance, double ci, NodeCompoundTreeItem *parent)
{
    // Single mass isotopomer
    role = MassIsotopomer;
    parentItem = parent;

    // add data
    // itemData.push_back(Qt::Checked);
    itemData.push_back(QString("M").append(QString::number(n)));
    itemData.push_back(QVariant(abundance));
    itemData.push_back(QVariant(ci));
}

NodeCompoundTreeItem::~NodeCompoundTreeItem()
{
    qDeleteAll(childItems);
}

void NodeCompoundTreeItem::appendChild(NodeCompoundTreeItem *child)
{
    childItems.append(child);
}

NodeCompoundTreeItem *NodeCompoundTreeItem::child(int row)
{
    return childItems.value(row);
}

int NodeCompoundTreeItem::childCount() const
{
    return childItems.count();
}

int NodeCompoundTreeItem::columnCount() const
{
    return itemData.count();
}

QVariant NodeCompoundTreeItem::data(int column, int role) const
{
    switch(role){
    case Qt::EditRole:
    case Qt::DisplayRole:
        if(column != CHECKBOX_USE_COLUMN) {
            return itemData.value(column);
        }
        break;
    case Qt::CheckStateRole:
        if(column == CHECKBOX_USE_COLUMN) {
            return itemData.value(column);
        }
        break;
    case Qt::TextAlignmentRole:
        /*if(column == 0) {
            return Qt::Checked;
        }*/
        break;
    case Qt::UserRole:
        return userRoleData;
    //case Qt::ForegroundRole:
//        if(nc && nc->)  return Qt::gray;
    }

    return QVariant();
}

bool NodeCompoundTreeItem::setData(int column, const QVariant &value, int role)
{
    itemData[column] = value; // TODO roles?
    return true;
}

int NodeCompoundTreeItem::row() const
{
    if(parentItem)
        return parentItem->childItems.indexOf(const_cast<NodeCompoundTreeItem*>(this));
    return 0;
}

NodeCompoundTreeItem *NodeCompoundTreeItem::parent()
{
    return parentItem;
}

}
