INCLUDE(LibFindMacros)

# Include dir
FIND_PATH(MD_INCLUDE_DIR
  NAMES chromatogramplot.h
  PATHS /usr/include/metabolitedetector /usr/local/include/metabolitedetector ~/src/metabolitedetector/src/
)

# Finally the library itself
FIND_LIBRARY(MD_LIBRARY
  NAMES metabolitedetector
  PATHS /usr/lib /usr/local/lib ~/src/metabolitedetector/build/src
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(MD_PROCESS_INCLUDES MD_INCLUDE_DIR)
set(MD_PROCESS_LIBS MD_LIBRARY)
libfind_process(MD)
