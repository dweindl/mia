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

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(CMAKE_PREFIX_PATH /usr/local/Qt-5.3.1/)
    if(POLICY CMP0043)
        cmake_policy(SET CMP0043 OLD)
    endif()
else ()
endif()


find_package(Qt5Widgets REQUIRED)
find_package(MD REQUIRED)
# find_package(GraphViz REQUIRED)
set(GraphViz_INCLUDE_DIR "/usr/local/include/graphviz" CACHE STRING "GraphViz include directory")

set(CMAKE_AUTOMOC TRUE)

set(GUI_SRC_LIST
    configdialog.cpp
    experimentlistwidget.cpp
    experimentwizard.cpp
    graphvizqt.cpp
    labelidentificatorqthread.cpp
    main.cpp
    miaguiconstants.h
    miamainwindow.cpp
    midplot.cpp
    nodecompounddescriptionwidget.cpp
    nodecompoundtreemodel.cpp
    nodewidget.cpp
    nwview.cpp
)

option(MIA_WITH_HDF5 "Include HDF5 support?" OFF)
if(MIA_WITH_HDF5)
    find_package(HDF5 REQUIRED)
    add_definitions(-DMIA_WITH_HDF5)
    set(GUI_SRC_LIST ${GUI_SRC_LIST} hdfwriter.cpp)
endif()

if(MIA_WITH_METABOBASE)
    set(GUI_SRC_LIST ${GUI_SRC_LIST} ../kegg_reaction_test/keggreactionmapper.cpp)
endif()

if(MIA_WITH_NETCDF_IMPORT)
    set(GUI_SRC_LIST ${GUI_SRC_LIST} netcdfimportdialog.cpp)
endif()

qt5_add_resources(GUI_RESOURCES gui-resources.qrc)

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ../src
    ${QWT_INCLUDE_DIR}
    ${MD_INCLUDE_DIR}
    ${LABID_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIRS}
    ${GraphViz_INCLUDE_DIR}
)

add_definitions(-DWITH_CGRAPH)

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(RES_FILES "win32/mia.rc")
    set(CMAKE_RC_COMPILER_INIT windres)
    ENABLE_LANGUAGE(RC)
    SET(CMAKE_RC_COMPILE_OBJECT "<CMAKE_RC_COMPILER> <FLAGS> -O coff <DEFINES> -i <SOURCE> -o <OBJECT>")
endif()

add_executable(mia-gui ${RES_FILES} ${GUI_SRC_LIST} ${GUI_RESOURCES})

qt5_use_modules(mia-gui Widgets Core Concurrent PrintSupport Svg Network)

################## BEGIN LINKING OPTIONS ##################
set (GraphViz_LIBRARY_PATH "$ENV{HOME}/src/_libs/graphviz-2.30.1" CACHE STRING "GraphViz root dir")
set (GraphViz_LIBRARIES 
  ${GraphViz_LIBRARY_PATH}/plugin/core/.libs/libgvplugin_core_C.a
  ${GraphViz_LIBRARY_PATH}/plugin/dot_layout/.libs/libgvplugin_dot_layout_C.a
  ${GraphViz_LIBRARY_PATH}/plugin/neato_layout/.libs/libgvplugin_neato_layout_C.a
  ${GraphViz_LIBRARY_PATH}/lib/gvc/.libs/libgvc_C.a
  ${GraphViz_LIBRARY_PATH}/lib/cgraph/.libs/libcgraph_C.a
  ${GraphViz_LIBRARY_PATH}/lib/cdt/.libs/libcdt_C.a
  ${GraphViz_LIBRARY_PATH}/lib/pathplan/.libs/libpathplan_C.a
  ${GraphViz_LIBRARY_PATH}/lib/cgraph/.libs/libcgraph_C.a
  ${GraphViz_LIBRARY_PATH}/lib/common/.libs/libcommon_C.a
)

set (GraphViz_LIBRARIES 
  ${GraphViz_LIBRARY_PATH}/graphviz/libgvplugin_core.so
  ${GraphViz_LIBRARY_PATH}/graphviz/libgvplugin_dot_layout.so
  ${GraphViz_LIBRARY_PATH}/graphviz/libgvplugin_neato_layout.so
  ${GraphViz_LIBRARY_PATH}/libgvc.so
  ${GraphViz_LIBRARY_PATH}/libcgraph.so
  ${GraphViz_LIBRARY_PATH}/libcdt.so
  ${GraphViz_LIBRARY_PATH}/libpathplan.so

)


set (Boost_LIBRARIES
    ${Boost_FILESYSTEM_LIBRARY_RELEASE}
    ${Boost_REGEX_LIBRARY_RELEASE}
    ${Boost_SYSTEM_LIBRARY_RELEASE}
)

if(CMAKE_SYSTEM_NAME STREQUAL "Windows")

    set(CROSSCOMPILE_ROOT "$ENV{HOME}/src/crosscompile/")
    set(STATIC_QT_LIBRARIES ${CROSSCOMPILE_ROOT}/mxe/usr/i686-w64-mingw32.static/qt5/plugins/platforms/libqwindows.a  -lwinspool -L${CROSSCOMPILE_ROOT}/mxe/usr/i686-w64-mingw32.static/qt5/lib -lQt5PlatformSupport -lQt5Gui -lcomdlg32 -loleaut32 -limm32 -lglu32 -lopengl32 -ljpeg -lharfbuzz -lcairo -lgobject-2.0 -lfontconfig -lusp10 -lmsimg32 -lgdi32 -lpixman-1 -lffi -lexpat -lfreetype -lbz2 -lpng16 -lharfbuzz_too -lglib-2.0 -lwinmm -lshlwapi -lpcre -lintl -liconv -lpng -lQt5Core -lole32 -luuid -lws2_32 -ladvapi32 -lshell32 -luser32 -lkernel32 -lz -lsicuin -lsicuuc -lsicudt -lpcre16)

    target_link_libraries(mia-gui
        ${GraphViz_LIBRARIES}
        ${CROSSCOMPILE_ROOT}/win32/graphviz-2.38.0/libltdl/.libs/libltdlc.a -ldl
        #${CROSSCOMPILE_ROOT}/win32/graphviz-2.38.0/libltdl/.libs/dlopen.a
        #  /usr/lib/x86_64-linux-gnu/libexpat.a
        ${CMAKE_BINARY_DIR}/src/libmia.a
        ${MD_LIBRARY}
    )

    if(MIA_WITH_HDF5)
        target_link_libraries(mia-gui
            #${HDF5_hdf5_LIBRARY_RELEASE}
            #${HDF5_C_LIBRARIES}
            #${HDF5_CXX_LIBRARIES}
            #${HDF5_HL_LIBRARIES}
            ${CROSSCOMPILE_ROOT}/mxe/usr/i686-w64-mingw32.static/lib/libhdf5.a
            ${CROSSCOMPILE_ROOT}/mxe/usr/i686-w64-mingw32.static/lib/libhdf5_cpp.a
        )
    endif()

  target_link_libraries(mia-gui
    ${LabId_LIBRARY}
    ${GCMS_LIBRARY}
    ${Boost_LIBRARIES}
    # ${POSTGRESQL_LIBRARY}
    ${ZLIB_LIBRARIES}
    ${GSL_LIBRARIES}
    # ${QWT_LIBRARY}
    ${STATIC_QT_LIBRARIES}
    -lharfbuzz -lharfbuzz_too  -lexpat -lpcre -lpcre16 -lsicuuc -lsicudt -lssl -lcrypt32 -lcrypto -lregex -lws2_32 -lmsvcr100
    )
    if(MIA_WITH_NETCDF_IMPORT)
        target_link_libraries(mia-gui
            ${NetCDF_LIBRARY}
        )
    endif()

else()
    set (STATIC_QT_LIBRARIES
        /usr/local/Qt-5.3.1/plugins/platforms/libqxcb.a
        # from libqxcb.prl: -lX11-xcb -lXi -lSM -lICE -ldbus-1 -lxcb -L~/src/_libs/qt5/qtbase/src/plugins/platforms/xcb/xcb-static -lxcb-static -lxkbcommon-x11 -lxkbcommon -L/usr/local/Qt-5.3.1/lib -lQt5PlatformSupport -lfontconfig -lfreetype -lXrender -lXext -lX11 -ludev -lQt5DBus -lQt5Gui -ljpeg -lpng -lEGL -lQt5Core -lz -licui18n -licuuc -licudata -lm -ldl -pthread -lgthread-2.0 -lglib-2.0 -lrt -lGL -lpthread
        -Wl,-Bstatic -lX11-xcb -lXi -lSM -lICE -ldbus-1 -lxcb -L/home/dweidl/src/_libs/qt5/qtbase/src/plugins/platforms/xcb/xcb-static -lxcb-static -lxkbcommon-x11 -lxkbcommon -L/usr/local/Qt-5.3.1/lib -lQt5PlatformSupport -lfontconfig -lfreetype -lXrender -lXext -lX11  -lQt5DBus -lQt5Gui -ljpeg -lpng -lQt5Core -lpcre16 -lz -licui18n -licuuc -licudata -lm -pthread -lgthread-2.0 -lglib-2.0 -lrt  -Wl,-Bdynamic -ludev -lEGL -lGL -lXau -lXdmcp -lltdl -ldl
    )

    target_link_libraries(mia-gui -Wl,-rpath=/usr/share/mia/ -Wl,-rpath=/usr/lib/graphviz/
        ${GraphViz_LIBRARIES}
        -lexpat
        mia
    )

    if(MIA_WITH_METABOBASE)
        target_link_libraries(mia-gui
            ${CMAKE_SOURCE_DIR}/../mddb/build/src/libMDDB.a
            #/usr/lib/x86_64-linux-gnu/libpq.a
            ${POSTGRESQL_LIBRARY}
        )
    endif()

    if(MIA_WITH_NETCDF_IMPORT)
        target_link_libraries(mia-gui
            ${NetCDF_LIBRARY}
        )
    endif()

    if(MIA_WITH_HDF5)
        target_link_libraries(mia-gui
            #${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_DEFINITIONS}
            /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a
            /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a
        )
    endif()

    target_link_libraries(mia-gui
#        ${MD_LIBRARY}
        ${LabId_LIBRARY}
        ${GCMS_LIBRARY}
        ${Boost_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${GSL_LIBRARIES}
        ${QWT_LIBRARY}
	${NetCDF_LIBRARY}
#        ${STATIC_QT_LIBRARIES}
    )
endif()

################## END LINKING OPTIONS ##################

################## BEGIN INSTALL OPTIONS ##################
install(TARGETS mia-gui DESTINATION bin)
################## END INSTALL OPTIONS ##################
