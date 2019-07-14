project (relate_lib)

add_library(relateStatic STATIC cxxopts data sample anc mutations)
target_link_libraries(relateStatic gzstreamStatic)
add_library(relateShared SHARED cxxopts data sample anc mutations)
target_link_libraries(relateShared gzstreamShared)
