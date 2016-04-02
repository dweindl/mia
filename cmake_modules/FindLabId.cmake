INCLUDE(LibFindMacros)

FIND_PATH(LabId_INCLUDE_DIR
  NAMES labeledcompound.h
  PATHS ~/src/labid/src
)
 
FIND_LIBRARY(LabId_LIBRARY
  NAMES liblabid.a
  PATHS ~/src/labid/build/src
)

set(LabId_PROCESS_INCLUDES LabId_INCLUDE_DIR )
set(LabId_PROCESS_LIBS LabId_LIBRARY )
libfind_process(LabId)

IF (LabId_INCLUDE_DIR AND LabId_LIBRARY)
   SET(LabId_FOUND TRUE)
ENDIF (LabId_INCLUDE_DIR AND LabId_LIBRARY)

IF (LabId_FOUND)
   IF (NOT LabId_FIND_QUIETLY)
      MESSAGE(STATUS "Found LabId: ${LabId_INCLUDE_DIR}")
      MESSAGE(STATUS "Found LabId: ${LabId_LIBRARY}")
   ENDIF (NOT LabId_FIND_QUIETLY)
ELSE (LabId_FOUND)
   IF (LabId_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find LabId")
   ENDIF (LabId_FIND_REQUIRED)
ENDIF (LabId_FOUND)
