find_program (DOXYGEN NAMES doxygen)
if(NOT DOXYGEN)
    message(WARNING "could not find doxygen, skipping documentation build")
endif()
find_program (MKDOCS NAMES mkdocs)
if(NOT MKDOCS)
    message(WARNING "could not find mkdocs, skipping documentation build")
endif()

execute_process(COMMAND ${MKDOCS} build
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc)
execute_process(COMMAND ${DOXYGEN} ../src/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc)

file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/doc/site/*.html")
foreach(file ${files})
    file (READ ${file} INDEX_HTML_IN)
    string (REPLACE "cpp-api.html" "doxygen/index.html" INDEX_HTML_OUT ${INDEX_HTML_IN})
    file (WRITE ${file} ${INDEX_HTML_OUT})
endforeach()

