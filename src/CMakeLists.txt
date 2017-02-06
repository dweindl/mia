# MIA - Mass Isotopolome Analyzer
# Copyright (C) 2012-15 Daniel Weindl <daniel@danielweindl.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

find_package(Qt5Concurrent REQUIRED)

set(SRC_LIST
    ../alg/alglibinternal.cpp
    ../alg/alglibmisc.cpp
    ../alg/ap.cpp
    ../alg/linalg.cpp
    ../alg/specialfunctions.cpp
    ../alg/statistics.cpp
    config.h
    labelingdataset.cpp
    labelingnetworkset.cpp
    miaexception.cpp
    middistancecalculator.cpp
    misc.cpp
    networklayer.cpp
    nodecompound.cpp
    serializationqt.cpp
    settings.cpp
    utilities.cpp
)

if(MIA_WITH_NETCDF_IMPORT)
  set(SRC_LIST ${SRC_LIST} compounddetector.cpp)
endif()

include_directories(
 ${CMAKE_CURRENT_SOURCE_DIR}
 ${CMAKE_CURRENT_BINARY_DIR}
)

add_library(mia ${SRC_LIST})

qt5_use_modules(mia Core Concurrent)


if(MIA_WITH_METABOBASE)
    target_link_libraries(mia ${POSTGRESQL_LIBRARY} ${METABOBASE_LIBRARY})
endif()

set (Boost_LIBRARIES
    ${Boost_FILESYSTEM_LIBRARY_RELEASE}
    ${Boost_REGEX_LIBRARY_RELEASE}
    ${Boost_SYSTEM_LIBRARY_RELEASE}
    ${Boost_IOSTREAMS_LIBRARY_RELEASE}
)

target_link_libraries(mia
    ${LabId_LIBRARY}
    ${GCMS_LIBRARY}
    ${Boost_LIBRARIES}
    ${ZLIB_LIBRARIES}
    ${GSL_LIBRARIES}
)
