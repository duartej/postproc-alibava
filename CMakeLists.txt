cmake_minimum_required(VERSION 2.8)

# Allowing to introduce the variable as option 
# reuquies to place it before project
set(CMAKE_INSTALL_PREFIX $ENV{HOME}/.local)

project(postproc-alibava)

# Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -ansi -pedantic -W -Wall -Wshadow -fPIC -ldl -std=c++11")

# The include directory
include_directories("${PROJECT_SOURCE_DIR}/include")
# and the source dir to build some useful libraries
add_subdirectory(src)

add_executable(fortythieves bin/fortythieves.cc)

target_link_libraries(fortythieves IOManager)

install(TARGETS fortythieves
    RUNTIME DESTINATION bin)
#install(PROGRAMS utils/processhzroot.py 
#    DESTINATION bin
#    RENAME processhzroot)
#install(PROGRAMS utils/hzplots.py 
#    DESTINATION bin
#    RENAME hzplots)

