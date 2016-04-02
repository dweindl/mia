INCLUDE(LibFindMacros)

# Include dir
FIND_PATH(GCMS_INCLUDE_DIR
  NAMES compound.h
  PATHS /usr/include/gcms /usr/local/include/gcms
)
FIND_PATH(andims_INCLUDE_DIR
  NAMES ms10.h
  PATHS /usr/include/andi-ms /usr/local/include/andims
)

# Finally the library itself
FIND_LIBRARY(GCMS_LIBRARY
  NAMES gcmslib
  PATHS /usr/local/lib /usr/lib
)
FIND_LIBRARY(andims_LIBRARY
  NAMES andims
  PATHS /usr/lib /usr/local/lib
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(GCMS_PROCESS_INCLUDES GCMS_INCLUDE_DIR andims_INCLUDE_DIR)
set(GCMS_PROCESS_LIBS GCMS_LIBRARY andims_LIBRARY)
libfind_process(GCMS)

IF (GCMS_INCLUDE_DIR AND GCMS_LIBRARY)
   SET(GCMS_FOUND TRUE)
ENDIF (GCMS_INCLUDE_DIR AND GCMS_LIBRARY)


IF (GCMS_FOUND)
   IF (NOT GCMS_FIND_QUIETLY)
      MESSAGE(STATUS "Found GCMS: ${GCMS_LIBRARY}")
   ENDIF (NOT GCMS_FIND_QUIETLY)
ELSE (GCMS_FOUND)
   IF (GCMS_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find GCMS")
   ENDIF (GCMS_FIND_REQUIRED)
ENDIF (GCMS_FOUND)
