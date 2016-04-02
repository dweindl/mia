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

#include "networklayer.h"

namespace mia {

NetworkLayer::NetworkLayer() : LabelingDataset()
{
    visible = true;
}

NetworkLayer::NetworkLayer(LabelingDataset *ds) : LabelingDataset(*ds)
{
    visible = true;
}

NetworkLayer::NetworkLayer(Settings s): LabelingDataset(s)
{
    visible = true;
}

bool NetworkLayer::isVisible() const
{
    return visible;
}

void NetworkLayer::setVisible(bool visibility)
{
    visible = visibility;
}

std::vector<NetworkLayer *> NetworkLayer::fromXMLFile(std::string file)
{
    std::vector<NetworkLayer *> l;
    std::vector<LabelingDataset *> datasets = LabelingDataset::fromXMLFile(file);
    for(int i = 0; i < datasets.size(); ++i)
        l.push_back(new NetworkLayer(datasets[i]));
    return l;
}

}
