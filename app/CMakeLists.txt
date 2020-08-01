cmake_minimum_required(VERSION 3.0)

project(dense3DCrfInference)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
# include(BasicConfig)
cond_find_package(mia REQUIRED)
cond_find_package(itkio REQUIRED)

find_package(Boost COMPONENTS program_options REQUIRED)

include_directories( include )
file(GLOB_RECURSE SRC "src/*.cpp")
file(GLOB_RECURSE HDR "include/*.h")
list(APPEND SRC ${HDR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_BUILD_TYPE Debug)

add_executable(${PROJECT_NAME} ${SRC})
# set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_PDB_NAME_DEBUG "${PROJECT_NAME}${CMAKE_DEBUG_POSTFIX}")

target_link_libraries(${PROJECT_NAME} PUBLIC ${ITK_LIBRARIES})
target_link_libraries(${PROJECT_NAME} PUBLIC ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(${PROJECT_NAME} PUBLIC "densecrf")
target_link_libraries(${PROJECT_NAME} PUBLIC "mia")
target_link_libraries(${PROJECT_NAME} PUBLIC "itkio")
