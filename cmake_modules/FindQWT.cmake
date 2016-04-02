# Find Qwt

# Once run this will define: 
# 
# QWT_FOUND       = system has QWT lib
#
# QWT_LIBRARY     = full path to the QWT library
#
# QWT_INCLUDE_DIR      = where to find headers 
#

FIND_PATH(QWT_INCLUDE_DIR qwt.h 
  PATHS /usr/local/qwt-6.1.0-rc3/include/ "$ENV{LIB_DIR}/include" "$ENV{LIB_DIR}/include/qwt"
  )

FIND_LIBRARY(QWT_LIBRARY qwt
              PATHS  /usr/local/qwt-6.1.0-rc3/lib/ )

IF (NOT QWT_LIBRARY)
   #try using ubuntu lib naming
  FIND_LIBRARY(QWT_LIBRARY 222 PATHS 
    /usr/local/qwt-5.2.2/lib
    )
ENDIF (NOT QWT_LIBRARY)

IF (QWT_INCLUDE_DIR AND QWT_LIBRARY)
  SET(QWT_FOUND TRUE)
ENDIF (QWT_INCLUDE_DIR AND QWT_LIBRARY)

IF (QWT_FOUND)
  IF (NOT QWT_FIND_QUIETLY)
    MESSAGE(STATUS "Found QWT: ${QWT_LIBRARY}")
    MESSAGE(STATUS "Found QWT: ${QWT_INCLUDE_DIR}")
  ENDIF (NOT QWT_FIND_QUIETLY)
ELSE (QWT_FOUND)
MESSAGE(FATAL_ERROR "QWT not found")
  IF (QWT_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find QWT")
  ENDIF (QWT_FIND_REQUIRED)
ENDIF (QWT_FOUND)
