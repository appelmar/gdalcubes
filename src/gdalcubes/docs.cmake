find_program (DOXYGEN NAMES doxygen)
if(NOT DOXYGEN)
    message(WARNING "could not find doxygen, skipping documentation build")
endif()

execute_process(COMMAND ${DOXYGEN} src/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})



