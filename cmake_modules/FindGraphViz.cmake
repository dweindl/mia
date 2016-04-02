include(LibFindMacros)

# Include dir
find_path(GraphViz_INCLUDE_DIR
  NAMES gvc.h
  PATHS /usr/local/include/graphviz/ /usr/include/graphviz/
)

# Library
#find_library(GraphViz_graph_LIBRARY
#  NAMES libgraph.so
#  PATHS /usr/local/lib /usr/lib
#)
find_library(GraphViz_gvc_LIBRARY
  NAMES libgvc.so
  PATHS /usr/local/lib /usr/lib
)
find_library(GraphViz_cgraph_LIBRARY
  NAMES libcgraph.so
  PATHS /usr/local/lib /usr/lib
)
find_library(GraphViz_cdt_LIBRARY
  NAMES libcdt.so
  PATHS /usr/local/lib /usr/lib
)


# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(GraphViz_PROCESS_INCLUDES GraphViz_INCLUDE_DIR)
set(GraphViz_PROCESS_LIBS  GraphViz_gvc_LIBRARY GraphViz_cdt_LIBRARY GraphViz_cgraph_LIBRARY)
libfind_process(GraphViz)
