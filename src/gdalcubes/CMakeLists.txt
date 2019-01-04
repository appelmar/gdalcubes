cmake_minimum_required(VERSION 3.8)
project(gdalcubes LANGUAGES CXX C VERSION 0.0.1)

#set(COLLECTION_FORMAT_VERSION_MAJOR 0)
#set(COLLECTION_FORMAT_VERSION_MINOR 0)
#set(COLLECTION_FORMAT_VERSION_PATCH 0)

#string(TIMESTAMP build_date "%Y%m%d")
#string(TIMESTAMP build_time "%H%M%S")

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type specified, using Release")
    set(CMAKE_BUILD_TYPE Release)
endif ()



#if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug") # enable profiling and coverage reports (gcc specific, only for debug)
#    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg --coverage")
#    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg --coverage" )
#    #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg --coverage")
#    #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=undefined -D_FORTIFY_SOURCE=0")
#    #set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  -fsanitize=address -fsanitize=undefined -D_FORTIFY_SOURCE=0")
#endif ( )


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
#add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
#add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)


if (DEFINED USE_EXPRTK)
    add_definitions(-DUSE_EXPRTK)
endif()


if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCURL_STATICLIB")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -DCURL_STATICLIB")
endif ()
if (MINGW AND USE_EXPRTK)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj") # for apply_pixel object due to exprtk library
endif ()



# cppcheck
find_program (CPPCHECK NAMES cppcheck)
if(${CPPCHECK-NOTFOUND})
    message(WARNING "could not find cppcheck, skipping static code analysis.")
else ()
    add_custom_target(cppcheck ALL
            COMMAND ${CPPCHECK} --enable=warning,performance  --language=c++ src 2> cppcheck.log
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            COMMENT "Running Cppcheck" VERBATIM)
    set_target_properties(cppcheck PROPERTIES EXCLUDE_FROM_ALL true)
endif()


# code formatting
find_program (CLANG_FORMAT NAMES clang-format)
if(NOT CLANG_FORMAT)
    message(WARNING "could not find clang-format, skipping automated source code formatting.")
else ()
    file(GLOB_RECURSE ALL_SOURCE_FILES  src/*.cpp src/*.h)
    # clang-format -i -style=file $
    add_custom_target(clangformat ALL
                      COMMAND ${CLANG_FORMAT} -i -style=file ${ALL_SOURCE_FILES}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                      COMMENT "Auto-formatting source code" VERBATIM)
    set_target_properties(clangformat PROPERTIES EXCLUDE_FROM_ALL true)
endif()





# documentation
find_program (DOXYGEN NAMES doxygen)
if(NOT DOXYGEN)
    message(WARNING "could not find doxygen, skipping documentation build")
endif()
find_program (MKDOCS NAMES mkdocs)
if(NOT MKDOCS)
    message(WARNING "could not find mkdocs, skipping documentation build")
endif()
if (DOXYGEN AND MKDOCS)
    file(GLOB_RECURSE DOC_FILES "doc/mkdocs.yml" "doc/docs/*" )
    add_custom_target(docs ALL
            COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/docs.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            DEPENDS ${DOC_FILES}
            COMMENT "Building documentation" VERBATIM)
    set_target_properties(docs PROPERTIES EXCLUDE_FROM_ALL true)
endif()


# install collection formats
if ( WIN32 )
    install(DIRECTORY formats/ DESTINATION $ENV{AllUsersProfile}/gdalcubes/formats FILES_MATCHING PATTERN "*.json")
else()
    install(DIRECTORY formats/ DESTINATION /usr/lib/gdalcubes/formats FILES_MATCHING PATTERN "*.json")
endif()


add_subdirectory(src)

# cpack stuff
set(CPACK_PACKAGE_NAME "gdalcubes")
set(CPACK_PACKAGE_VENDOR "Marius Appel <marius.appel@uni-muenster.de>")
set(CPACK_PACKAGE_CONTACT "Marius Appel <marius.appel@uni-muenster.de>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "gdalcubes - Earth observation data cubes from GDAL image collections")
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE)
set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.md)

include(CPack)