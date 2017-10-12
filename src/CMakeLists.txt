# Some useful libraries
add_library(IOManagersAlibava SHARED IOManager.cxx IOAlibavaReader.cxx AlibavaPostProcessor.cxx AlibavaDiagnosis.cxx)
# The fortythieves post-processor library
add_library(AlibavaSensorAnalysis SHARED IOFortythieves.cxx)

# Make sure the compiler can find include files for the library
# when other libraries or executables link to IOManager
target_include_directories(IOManagersAlibava PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(AlibavaSensorAnalysis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# ROOT dependencies
set(ROOTLIBSTOLINK Tree RIO)
find_package(ROOT REQUIRED ${ROOTLIBSTOLINK})
include_directories(${ROOT_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ROOT_LIBRARIES})

target_link_libraries(IOManagersAlibava ${LIBS})
target_link_libraries(AlibavaSensorAnalysis ${LIBS})
# 
install(
    TARGETS IOManagersAlibava 
    RUNTIME DESTINATION lib
    LIBRARY DESTINATION lib
    COMPONENT library)

install(
    TARGETS AlibavaSensorAnalysis
    RUNTIME DESTINATION lib
    LIBRARY DESTINATION lib
    COMPONENT library)
