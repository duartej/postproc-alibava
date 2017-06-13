cmake_minimum_required(VERSION 2.8)

project(postproc-alibava)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "default install path" FORCE )
endif()
MESSAGE(STATUS "Installation path at: ${CMAKE_INSTALL_PREFIX}")

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
