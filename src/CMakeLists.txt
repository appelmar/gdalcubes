

execute_process(
        COMMAND git describe --tags --always --abbrev=8
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_DESC
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
        COMMAND git rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE
)


configure_file( ${CMAKE_CURRENT_LIST_DIR}/build_info.h.in ${CMAKE_CURRENT_LIST_DIR}/build_info.h )


## find GDAL
find_path(GDAL_INCLUDEDIR "gdal.h" HINTS $ENV{GDAL_INCLUDEDIR} PATH_SUFFIXES gdal  NO_DEFAULT_PATH)
find_path(GDAL_INCLUDEDIR "gdal.h" HINTS $ENV{GDAL_INCLUDEDIR} PATH_SUFFIXES gdal)
if (GDAL_INCLUDEDIR-NOTFOUND)
    message(FATAL_ERROR "GDAL headers not found")
endif ()

find_library(GDAL_LIBRARY NAMES gdal gdal_i HINTS $ENV{GDAL_LIBDIR} NO_DEFAULT_PATH)
find_library(GDAL_LIBRARY NAMES gdal gdal_i HINTS $ENV{GDAL_LIBDIR})
if (NOT GDAL_LIBRARY)
    message(FATAL_ERROR "GDAL library not found!")
else ()
    message(STATUS "GDAL library found at ${GDAL_LIBRARY}")
endif ()



#set(Boost_DEBUG ON)
## find boost
#if (WIN32)
    set(BOOST_USE_STATIC_LIBS ON)
#endif ()




## find SQLite
find_path(sqlite_INCLUDEDIR "sqlite3.h" HINTS $ENV{SQLITE_INCLUDEDIR} NO_DEFAULT_PATH)
find_path(sqlite_INCLUDEDIR "sqlite3.h" HINTS $ENV{SQLITE_INCLUDEDIR})
if (sqlite_INCLUDEDIR-NOTFOUND)
    message(FATAL_ERROR "SQLite headers not found")
endif ()

find_library(SQLITE_LIBRARY NAMES sqlite sqlite3 HINTS $ENV{SQLITE_LIBDIR} NO_DEFAULT_PATH)
find_library(SQLITE_LIBRARY NAMES sqlite sqlite3 HINTS $ENV{SQLITE_LIBDIR})
if (NOT SQLITE_LIBRARY)
    message(FATAL_ERROR "SQLite library not found!")
else ()
    message(STATUS "SQLite library found at ${SQLITE_LIBRARY}")
endif ()



find_package(Threads)



## find netcdf C library
find_path(NETCDF_INCLUDEDIR "netcdf.h" HINTS $ENV{NETCDF_INCLUDEDIR} NO_DEFAULT_PATH)
find_path(NETCDF_INCLUDEDIR "netcdf.h" HINTS $ENV{NETCDF_INCLUDEDIR})
if (NETCDF_INCLUDEDIR-NOTFOUND)
    message(FATAL_ERROR "netcdf C library headers not found")
endif ()

find_library(NETCDF_LIBRARY netcdf HINTS $ENV{NETCDF_LIBDIR} NO_DEFAULT_PATH)
find_library(NETCDF_LIBRARY netcdf HINTS $ENV{NETCDF_LIBDIR})
if (NOT NETCDF_LIBRARY)
    message(FATAL_ERROR "netcdf C library not found!")
else ()
    message(STATUS "netcdf C library found at ${NETCDF_LIBRARY}")
endif ()



# find libcurl
find_path(CURL_INCLUDEDIR "curl.h" HINTS $ENV{CURL_INCLUDEDIR} PATH_SUFFIXES curl NO_DEFAULT_PATH)
find_path(CURL_INCLUDEDIR "curl.h" HINTS $ENV{CURL_INCLUDEDIR} PATH_SUFFIXES curl)
if (CURL_INCLUDEDIR-NOTFOUND)
    message(FATAL_ERROR "libcurl header not found")
endif ()

find_library(CURL_LIBRARY NAMES curl curl-dualssl HINTS $ENV{CURL_LIBDIR} NO_DEFAULT_PATH)
find_library(CURL_LIBRARY NAMES curl curl-dualssl HINTS $ENV{CURL_LIBDIR})
if (NOT CURL_LIBRARY)
    message(FATAL_ERROR "libcurl not found!")
else ()
    message(STATUS "libcurl found at ${CURL_LIBRARY}")
endif ()





include_directories (${GDAL_INCLUDEDIR} ${CURL_INCLUDEDIR} ${sqlite_INCLUDEDIR} ${NETCDF_INCLUDEDIR})

get_filename_component(LDIR_1  ${NETCDF_LIBRARY} DIRECTORY)
get_filename_component(LDIR_2  ${SQLITE_LIBRARY} DIRECTORY)
get_filename_component(LDIR_3  ${GDAL_LIBRARY}  DIRECTORY)
get_filename_component(LDIR_4  ${CURL_LIBRARY} DIRECTORY)
link_directories(${LDIR_1} ${LDIR_2} ${LDIR_3} ${LDIR_4})


# set(ADD_LIBS "-lgdal -lsqlite3 -lspatialite -lproj -lgeos_c -lgeos -ljson-c -lnetcdf -lmariadbclient -lpq -lintl -lwebp -lcurl -lssh2 -lssl -lcrypto -lkea -lhdf5_cpp -lhdf5_hl -lhdf5 -lexpat -lfreexl -lcfitsio -lmfhdf -ldf -lxdr -lopenjp2 -ljasper -lpng16 -ljpeg -ltiff -lgeotiff -lgif -lxml2 -llzma -lszip -lz -lodbc32 -lodbccp32 -liconv -lpsapi -lws2_32 -lcrypt32 -lwldap32 -lsecur32 -lgdi32 -lssh2 -lssl -lcrypto")




file(GLOB SOURCE_FILES
    "*.cpp" "*.h"
)

if (NOT DEFINED USE_EXPRTK)
    list(APPEND SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/external/tinyexpr/tinyexpr.c)
endif()

add_subdirectory("external/tiny-process-library")



list(REMOVE_ITEM SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/example.cpp)
list(REMOVE_ITEM SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/gdalcubes.cpp)
list(REMOVE_ITEM SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/server.cpp)



add_library(libgdalcubes_static STATIC  ${SOURCE_FILES})
#set_property(TARGET libgdalcubes_static PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(libgdalcubes_static PROPERTIES OUTPUT_NAME "gdalcubes")
target_link_libraries(libgdalcubes_static ${GDAL_LIBRARY} ${SQLITE_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${NETCDF_LIBRARY} ${CURL_LIBRARY} tiny-process-library)


add_library(libgdalcubes_shared SHARED  ${SOURCE_FILES})
#set_property(TARGET libgdalcubes_shared PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(libgdalcubes_shared PROPERTIES OUTPUT_NAME "gdalcubes")
target_link_libraries(libgdalcubes_shared ${GDAL_LIBRARY} ${SQLITE_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${NETCDF_LIBRARY} ${CURL_LIBRARY} tiny-process-library)


install(TARGETS libgdalcubes_static libgdalcubes_shared RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib/static)


# install header files
install(DIRECTORY . DESTINATION include/gdalcubes
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")


file(GLOB TEST_FILES
        "test/*.cpp" "*.h")

add_executable(gdalcubes_test ${TEST_FILES})
target_link_libraries (gdalcubes_test libgdalcubes_static)


find_package(Boost 1.65 COMPONENTS program_options system) # system is required for error codes
if (Boost_FOUND)
    message(STATUS "Found Boost libraries ${Boost_LIBRARIES}")


    add_executable(gdalcubes ${CMAKE_CURRENT_SOURCE_DIR}/gdalcubes.cpp)
    target_link_libraries (gdalcubes libgdalcubes_static ${Boost_LIBRARIES})
    install(TARGETS gdalcubes RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib/static)
    target_include_directories(gdalcubes PRIVATE ${Boost_INCLUDE_DIRS})


    add_executable(gdalcubes_example ${CMAKE_CURRENT_SOURCE_DIR}/example.cpp)
    target_link_libraries (gdalcubes_example libgdalcubes_static ${Boost_LIBRARIES})
    target_include_directories(gdalcubes_example PRIVATE ${Boost_INCLUDE_DIRS})


    find_library(CPPRESTSDK_LIB cpprest)
    find_library(OPENSSL_LIB ssl)
    find_library(CRYPTO_LIB crypto)
    if (NOT CPPRESTSDK_LIB OR NOT OPENSSL_LIB OR NOT CRYPTO_LIB)
        message(WARNING "Skipping build of gdalcubes_worker due to missing libraries.")
    else ()
        message(STATUS "cpprestsdk library found at ${CPPRESTSDK_LIB}")
        message(STATUS "openssl library found at ${OPENSSL_LIB}")
        message(STATUS "crypto library found at ${CRYPTO_LIB}")
        add_executable(gdalcubes_server ${CMAKE_CURRENT_SOURCE_DIR}/server.cpp)
        target_link_libraries (gdalcubes_server libgdalcubes_static ${CPPRESTSDK_LIB} ${OPENSSL_LIB} ${CRYPTO_LIB} ${Boost_LIBRARIES})
        target_include_directories(gdalcubes_server PRIVATE ${Boost_INCLUDE_DIRS})
        install(TARGETS gdalcubes_server RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib/static)
    endif ()


else ()
    message(WARNING "Could not find Boost libraries, skipping build of gdalcubes_server and command line interface")
endif ()









