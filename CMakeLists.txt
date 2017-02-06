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

project(mia)

cmake_minimum_required(VERSION 2.8)

# Set version
execute_process(COMMAND sh -c "git describe --abbrev=4 --dirty=-dirty --always --tags | cut -c2- | tr -d '\n'" OUTPUT_VARIABLE GIT_VERSION)
message( STATUS "Building version ${GIT_VERSION}")
add_definitions(-DMIA_VERSION="${GIT_VERSION}")

MESSAGE( STATUS "CMAKE_SYSTEM: " ${CMAKE_SYSTEM} )
MESSAGE( STATUS "CMAKE_SYSTEM_NAME: " ${CMAKE_SYSTEM_NAME} )
MESSAGE( STATUS "CMAKE_SYSTEM_VERSION: " ${CMAKE_SYSTEM_VERSION} )
MESSAGE( STATUS "CMAKE_SYSTEM_PROCESSOR: " ${CMAKE_SYSTEM_PROCESSOR} )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-s -O3") # -O3
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")
set(CMAKE_AUTOMOC TRUE)

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    if(POLICY CMP0020)
        cmake_policy(SET CMP0020 NEW)
    endif()
endif()

option(STATIC_LINKING "Link statically?" ON)
if(STATIC_LINKING)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(Boost_USE_STATIC_LIBS ON)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DMIA_DEBUG_LEVEL=1)
endif()

option(MIA_WITH_METABOBASE "Use METABOBASE interface?" OFF)
if(MIA_WITH_METABOBASE)
    add_definitions(-DMIA_WITH_METABOBASE) # METABOBASE extension
    set(METABOBASE_LIBRARY "$ENV{HOME}/src/mddb/build/src/libMDDB.a")
    set(METABOBASE_INCLUDE_DIR "$ENV{HOME}/src/mddb/src/")
    find_package(PQ REQUIRED)
    include_directories(${METABOBASE_INCLUDE_DIR})
endif()

option(MIA_WITH_NETCDF_IMPORT "Add netCDF import?" ON)
if(MIA_WITH_NETCDF_IMPORT)
    add_definitions(-DMIA_WITH_NETCDF_IMPORT)
    find_package(NetCDF REQUIRED)
endif()

find_package(LabId REQUIRED)
find_package(GSL REQUIRED)
find_package(GCMS REQUIRED)
#find_package(GraphViz REQUIRED)
find_package(Boost COMPONENTS regex system filesystem REQUIRED)
find_package(ZLIB REQUIRED)

include_directories(
 ${CMAKE_CURRENT_SOURCE_DIR}
 ${CMAKE_CURRENT_BINARY_DIR}
 ${POSTGRESQL_INCLUDE_DIR}
 ${NetCDF_INCLUDE_DIR}
 ${GCMS_INCLUDE_DIR}
 ${GSL_INCLUDE_DIR}
 ${LabId_INCLUDE_DIR}
)

add_subdirectory(src) # Lib
add_subdirectory(gui) # GUI

#############################
# BEGIN CPACK configuration #
#############################
set(CPACK_PACKAGE_VERSION ${GIT_VERSION})
set(CPACK_PACKAGE_CONTACT "Daniel WEINDL <sci@danielweindl.de>")
set(CPACK_PACKAGE_NAME "mia")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${GIT_VERSION}")
set(CPACK_PACKAGE_EXECUTABLES "mia-gui" "MIA - Mass Isotopolome Analyzer")

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(CPACK_GENERATOR "NSIS")
    set(CPACK_NSIS_MUI_ICON "${CMAKE_SOURCE_DIR}/gui/icons/programmIcon16x16.ico")
    set(CPACK_NSIS_MUI_UNIICON ${CPACK_NSIS_MUI_ICON})
    set(CPACK_NSIS_HELP_LINK "http://massisotopolomeanalyzer.lu/")
    set(CPACK_NSIS_URL_INFO_ABOUT "http://massisotopolomeanalyzer.lu/")
    set(CPACK_NSIS_CONTACT "sci@danielweindl.de")
    set(CPACK_NSIS_CREATE_ICONS "CreateShortCut '\$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\MIA.lnk' '\$INSTDIR\\\\bin\\\\mia-gui.exe'")
    set(CPACK_NSIS_MUI_FINISHPAGE_RUN "mia-gui.exe")
    install(FILES "${CMAKE_SOURCE_DIR}/doc/${CPACK_PACKAGE_NAME}-doc.pdf" DESTINATION "doc\\\\")
    install(FILES "${CMAKE_SOURCE_DIR}/win32/msvcr100.dll" DESTINATION "bin\\\\")
else()
    set(CPACK_GENERATOR "DEB") # ;TGZ
    exec_program("dpkg --print-architecture" OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE)
    set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${GIT_VERSION}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6,libpq5,libgsl0ldbl,libx11-xcb1,libxkbcommon-x11-0,libjpeg8,libpng12-0") # objdump -p gui/mia-gui | grep NEEDED
    set(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6,libqt5core5a,libqt5gui5,libqt5widgets5,libnetcdf11,libgsl2,libboost-regex1.58.0,libgvc6,libqt5svg5")
    set(CPACK_DEBIAN_PACKAGE_SECTION "Science")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Mass Isotopolome Analyzer
    Non-targeted stable isotope labeling analysis tool")
    set(CPACK_DEBIAN_PACKAGE_SUGGESTS "metabolitedetector")
    install(FILES "${CMAKE_SOURCE_DIR}/doc/${CPACK_PACKAGE_NAME}-doc.pdf" DESTINATION "share/doc/${CPACK_PACKAGE_NAME}")
    install(FILES "${CMAKE_SOURCE_DIR}/${CPACK_PACKAGE_NAME}.xpm" DESTINATION "share/pixmaps")
    install(FILES "${CMAKE_SOURCE_DIR}/deb/${CPACK_PACKAGE_NAME}.desktop" DESTINATION "share/applications")
    install(FILES "${CMAKE_SOURCE_DIR}/deb/copyright" DESTINATION "share/doc/${CPACK_PACKAGE_NAME}")
    exec_program("cd ${CMAKE_SOURCE_DIR}/deb; gzip -k9f mia-gui.1; gzip -k9f changelog;")
    install(FILES "${CMAKE_SOURCE_DIR}/deb/mia-gui.1.gz" DESTINATION "share/man/man1/")
    install(FILES "${CMAKE_SOURCE_DIR}/deb/libstdc++.so.6" DESTINATION "share/${CPACK_PACKAGE_NAME}")
    install(FILES "${CMAKE_SOURCE_DIR}/deb/changelog.gz" DESTINATION "share/doc/${CPACK_PACKAGE_NAME}")
endif()
include(CPack)
###########################
# END CPACK configuration #
###########################
