INCLUDE(LibFindMacros)

FIND_PATH(InChi_INCLUDE_DIR
  NAMES inchi_api.h
  PATHS /home/dweindl/software/Inchi/1.04/INCHI-1-API/INCHI_API/inchi_dll/
)
 
FIND_LIBRARY(InChi_LIBRARY
   NAMES libinchi.so.1
   PATHS /home/dweindl/software/Inchi/1.04/INCHI-1-API/INCHI_API/gcc_so_makefile/result/
)

set(InChi_PROCESS_INCLUDES InChi_INCLUDE_DIR )
set(InChi_PROCESS_LIBS InChi_LIBRARY )
libfind_process(InChi)


IF (InChi_INCLUDE_DIR AND InChi_LIBRARY)
   SET(InChi_FOUND TRUE)
ENDIF (InChi_INCLUDE_DIR AND InChi_LIBRARY)

IF (InChi_FOUND)
   IF (NOT InChi_FIND_QUIETLY)
      MESSAGE(STATUS "Found InChi: ${InChi_INCLUDE_DIR}")
      MESSAGE(STATUS "Found InChi: ${InChi_LIBRARY}")
   ENDIF (NOT InChi_FIND_QUIETLY)
ELSE (InChi_FOUND)
   IF (InChi_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find InChi")
   ENDIF (InChi_FIND_REQUIRED)
ENDIF (InChi_FOUND)
