######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.0)
project(sv_align VERSION 0.0.0)

# Compiler options.
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_CXX_FLAGS "-fexceptions -fno-new-ttp-matching -fsanitize=address")
set(CMAKE_CXX_FLAGS "-fexceptions -fno-new-ttp-matching")

#add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization
add_definitions(-O3 -Wall)              # Much optimization

# Definitions needed to eliminate runtime dependency
# on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test executables.
set(SOURCES
        src/BinaryIO.cpp
        src/FastaReaderLite.cpp
        src/GFAReader.cpp
        src/VCFReader.cpp
        )


project(sv_align)
add_library(sv_align SHARED ${SOURCES})

# To make sure the library is named runlength_analysis.so,
# get rid of the "lib" prefix.
set_target_properties(sv_align PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(sv_align PROPERTIES DEFINE_SYMBOL "")

############################################
# ---------------------------------------- #
# -------- Generating executables -------- #
# ---------------------------------------- #
############################################

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# stupid "experimental" filesystem library for gcc <8.0
link_libraries(stdc++fs)

find_package( Boost 1.60 COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )

find_package(VGio REQUIRED)


# -------- TESTS --------

set(FILENAME_PREFIX test_VCFReader)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN")
target_link_libraries(${FILENAME_PREFIX} sv_align Threads::Threads ${Boost_LIBRARIES} stdc++fs VGio::VGio)

set(FILENAME_PREFIX test_FastaReaderLite)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN")
target_link_libraries(${FILENAME_PREFIX} sv_align Threads::Threads ${Boost_LIBRARIES} stdc++fs VGio::VGio)

set(FILENAME_PREFIX test_GFAReader)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN")
target_link_libraries(${FILENAME_PREFIX} sv_align Threads::Threads ${Boost_LIBRARIES} stdc++fs VGio::VGio)

set(FILENAME_PREFIX test_VGio)
add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN")
target_link_libraries(${FILENAME_PREFIX} sv_align Threads::Threads ${Boost_LIBRARIES} stdc++fs VGio::VGio)

# -------- EXECUTABLES --------

set(FILENAME_PREFIX generate_haploblocks_from_vcf)
add_executable(${FILENAME_PREFIX} src/executables/${FILENAME_PREFIX}.cpp)
set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN")
target_link_libraries(${FILENAME_PREFIX} sv_align Threads::Threads ${Boost_LIBRARIES} stdc++fs)

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        sv_align
        DESTINATION ${BINARY_INSTALL_DIR})

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

