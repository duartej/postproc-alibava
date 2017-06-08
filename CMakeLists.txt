cmake_minimum_required(VERSION 2.8)

# Allowing to introduce the variable as option 
# reuquies to place it before project
set(CMAKE_INSTALL_PREFIX $ENV{HOME}/.local)

project(postproc-alibava)

# Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -ldl -std=c++11")

add_executable(fortythieves src/fortythieves.cc)

# ROOT dependencies
set(ROOTLIBSTOLINK Tree RIO)
find_package(ROOT REQUIRED ${ROOTLIBSTOLINK})
include_directories(${ROOT_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ROOT_LIBRARIES})

target_link_libraries(fortythieves ${LIBS})

install(TARGETS fortythieves
    RUNTIME DESTINATION bin)
#install(PROGRAMS utils/processhzroot.py 
#    DESTINATION bin
#    RENAME processhzroot)
#install(PROGRAMS utils/hzplots.py 
#    DESTINATION bin
#    RENAME hzplots)

