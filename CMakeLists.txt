cmake_minimum_required (VERSION 3.0)
project(denseCrf3d)


find_package( ITK REQUIRED)
add_subdirectory(external)
add_subdirectory(external/mia)
add_subdirectory(external/itkio)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)




list(APPEND SRC ${HDR})

INCLUDE( ${ITK_USE_FILE} )


include_directories( include )
# file(GLOB_RECURSE HDR "include/*.h")


add_subdirectory (src)
add_subdirectory (app)
