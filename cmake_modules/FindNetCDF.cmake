FIND_PATH(NetCDF_INCLUDE_DIR netcdf.h /usr/include/ /usr/local/include/)

FIND_LIBRARY(NetCDF_LIBRARY NAMES netcdf PATHS /usr/lib /usr/local/lib) 

IF (NetCDF_INCLUDE_DIR AND NetCDF_LIBRARY)
   SET(NetCDF_FOUND TRUE)
ENDIF (NetCDF_INCLUDE_DIR AND NetCDF_LIBRARY)


IF (NetCDF_FOUND)
   IF (NOT NetCDF_FIND_QUIETLY)
      MESSAGE(STATUS "Found NetCDF: ${NetCDF_LIBRARY}")
   ENDIF (NOT NetCDF_FIND_QUIETLY)
ELSE (NetCDF_FOUND)
   IF (NetCDF_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find NetCDF")
   ENDIF (NetCDF_FIND_REQUIRED)
ENDIF (NetCDF_FOUND)