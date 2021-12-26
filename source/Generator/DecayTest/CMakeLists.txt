#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(DecayTestProject)

#add_definitions(-std=c++11)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIR})

include_directories(${PROJECT_SOURCE_DIR}/../include) 

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
file(GLOB_RECURSE sources "${PROJECT_SOURCE_DIR}/../src/*.cc" )
file(GLOB_RECURSE headers "${PROJECT_SOURCE_DIR}/../include/*.hh")

add_executable(DecayTest main.cpp ${sources} ${headers})
target_link_libraries(DecayTest ${Geant4_LIBRARIES} )
target_link_libraries(DecayTest ${ROOT_LIBRARIES})



#
add_custom_target(DecayTestProject DEPENDS DecayTest)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#

