# ---- Create the mechanism to integrate the alibavaSkifftools 
#      into the cmake build sytem, therefore
find_program(PYTHON "python")


if(PYTHON)
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(SETUP_CFG   "${CMAKE_CURRENT_BINARY_DIR}/setup.cfg")
    set(DEPS        "${CMAKE_CURRENT_SOURCE_DIR}/python/__init__.py")
    set(OUTPUT      "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")

    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    add_custom_command(OUTPUT ${OUTPUT}
                       COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/setup.cfg ${SETUP_CFG}
                       COMMAND ${PYTHON} ${SETUP_PY} build
                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
                       DEPENDS ${DEPS})

    add_custom_target(target ALL DEPENDS ${OUTPUT})
    
    # -- Not installed in the user area
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(USER_PATH_OPT "--user")
    else()
        set(USER_PATH_OPT "--prefix ${CMAKE_INSTALL_PREFIX}")
    endif()

    install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install ${USER_PATH_OPT})")
endif()
