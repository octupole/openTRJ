cmake_minimum_required( VERSION 3.0)


# Project name
project(test NONE)
find_package(PythonInterp)

set(CMD_PIP "NONE")
if(PYTHONINTERP_FOUND)
    if(PYTHON_VERSION_STRING GREATER_EQUAL 3.4)
        find_path(PIP_PATH pip3)
        if(${PIP_PATH} STREQUAL "PIP_PATH-NOTFOUND")
            find_path(PIP_PATH pip)
            if(NOT ${PIP_PATH} STREQUAL "PIP_PATH-NOTFOUND")
                set(CMD_PIP_TEST ${PIP_PATH}/pip --version)
                execute_process(COMMAND ${CMD_PIP_TEST}
                RESULT_VARIABLE VAR_OUTPUT
                OUTPUT_VARIABLE CMD_OUTPUT)
                if(VAR_OUTPUT EQUAL 0)
                    string(REGEX MATCH "python *3" HAS_MATCHED ${CMD_OUTPUT} )
                    if(HAS_MATCHED)
                        set(CMD_PIP ${PIP_PATH}/pip)
                    endif()
                endif()
            endif()
        else()
            set(CMD_PIP ${PIP_PATH}/pip3)
        endif()
    else()
        message(WARNING "Your version of python is ${PYTHON_VERSION_STRING} < 3.4")
    endif()
else()
    message(WARNING "pip not found cannot install python scripts")
endif()
if(NOT CMD_PIP STREQUAL "NONE")
    set(CMD_PIP "${CMD_PIP} install -e ${CMAKE_SOURCE_DIR}/python --user")
    message(STATUS "Installing python file with: ${CMD_PIP}")
    install(CODE "execute_process(COMMAND ${CMD_PIP} )")
    install(CODE "execute_process(COMMAND ln -f -s ~/.local/bin/Rg ${CMAKE_INSTALL_PREFIX} )")
    install(CODE "execute_process(COMMAND ln -f -s ~/.local/bin/Voro ${CMAKE_INSTALL_PREFIX})")
    install(CODE "MESSAGE(STATUS \"Python3 scripts installed in ${CMAKE_INSTALL_PREFIX}/bin\")")
endif()


