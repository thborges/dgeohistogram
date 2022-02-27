cmake_minimum_required(VERSION 3.0.0)
cmake_policy(SET CMP0054 NEW)

project(dgeohistogram)

find_library(GEOS_C NAMES geos_c 
    PATHS 
    /usr/lib64 
    /usr/lib 
    /usr/local/lib )

IF(GEOS_C)
    MESSAGE ("GEOS C lib: ${GEOS_C}")
    get_filename_component(GEOS_DIR ${GEOS_C} DIRECTORY)
    include_directories(${GEOS_DIR}/../include)
    link_directories(${GEOS_DIR})
    MESSAGE ("GEOS C include: ${GEOS_DIR}/../include")
    SET(GEOS_C_FOUND TRUE)
ENDIF(GEOS_C)

find_library(GDAL NAMES gdal 
    NO_DEFAULT_PATH
    PATHS 
    /usr/local/lib
    /usr/lib64 
    /usr/lib )


IF(GDAL)
    MESSAGE ("GDAL lib: ${GDAL}")
    get_filename_component(GDAL_DIR ${GDAL} DIRECTORY)
    include_directories(${GDAL_DIR}/../include)
    link_directories(${GDAL_DIR})
    MESSAGE ("GDAL include: ${GDAL_DIR}/../include")
    SET(GDAL_FOUND TRUE)
ENDIF(GDAL)

IF (NOT GEOS_C_FOUND OR NOT GDAL_FOUND)
MESSAGE (SEND_ERROR "Some dependencies were not found! Geos: ${GEOS_C_FOUND}, GDAL: ${GDAL_FOUND}")
ENDIF (NOT GEOS_C_FOUND OR NOT GDAL_FOUND)

set(CMAKE_CXX_STANDARD 17)

include_directories(include rtree)

file(GLOB wquery_SRC CONFIGURE_DEPENDS
    "experiments/wquery.cpp"
    "src/*.cpp"
    "rtree/*.c"
)

add_executable(wquery ${wquery_SRC})
target_link_libraries(wquery gdal geos geos_c)

set_target_properties(wquery PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    VS_DEBUGGER_COMMAND_ARGUMENTS "datasets/hidrografia.shp")