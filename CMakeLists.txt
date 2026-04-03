# CMake minimum version and project name
cmake_minimum_required(VERSION 3.14)
project(columba)

# Enforce compiler restrictions
if(WIN32 AND CMAKE_COMPILER_IS_MINGW)
    message(FATAL_ERROR "MinGW-w64 is not supported. Please use a different compiler.")
endif()

# Set default output directories early
if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
endif()
if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
endif()
if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif()

# ---- Configuration Options ----

# Option to enable run-length compression (default: OFF)
option(RUN_LENGTH_COMPRESSION "Enable run-length compression" OFF)

# Option to enable move structures for the phi operations (default: OFF)
option(PHI_MOVE "Enable phi move structures" OFF)
if(PHI_MOVE AND NOT RUN_LENGTH_COMPRESSION)
    message(WARNING "PHI_MOVE is not supported with Vanilla (non-RLC) mode and will be ignored.")
endif()

# Option to enable developer mode (default: OFF)
option(DEVELOPER_MODE "Enable developer mode" OFF)
if(DEVELOPER_MODE)
    add_definitions(-DDEVELOPER_MODE)
    message(STATUS "Developer mode enabled.")
else()
    message(STATUS "Developer mode disabled.")
endif()

# Compiler flags for warnings and C++ standard
set(COMMON_FLAGS "-Wall -pedantic -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -std=gnu++11")
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mpopcnt" HAS_POPCNT)
if(HAS_POPCNT)
    set(COMMON_FLAGS "${COMMON_FLAGS} -mpopcnt")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")

# Set default 32-bit or 64-bit type based on RUN_LENGTH_COMPRESSION
if(RUN_LENGTH_COMPRESSION)
    set(THIRTY_TWO_DEFAULT OFF)
    add_definitions(-DRUN_LENGTH_COMPRESSION)
    message(STATUS "Run-length compression enabled: using the bidirectional move structure (b-move).")

    if(PHI_MOVE)
        add_definitions(-DPHI_MOVE)
        message(STATUS "Phi move structures enabled.")
    else()
        message(STATUS "Phi move structures disabled.")
    endif()
else()
    set(THIRTY_TWO_DEFAULT ON)
    message(STATUS "Run-length compression disabled: using the bidirectional FM-index.")
endif()

option(THIRTY_TWO "Use 32-bit types" ${THIRTY_TWO_DEFAULT})
if(THIRTY_TWO)
    add_definitions(-DTHIRTY_TWO)
    message(STATUS "Using 32-bit types.")
else()
    message(STATUS "Using 64-bit types.")
endif()

# Always define ALPHABET as 5
add_definitions(-DALPHABET=5)

# ---- Restrict Windows to non-RLC mode ----
if(WIN32)
    if(RUN_LENGTH_COMPRESSION)
        message(FATAL_ERROR "Run-length compression is not supported on Windows. Please disable RUN_LENGTH_COMPRESSION.")
    endif()
endif()

# ---- Source Files ----

set(SOURCE_FILES
    src/indexhelpers.cpp
    src/bitparallelmatrix.cpp
    src/searchstrategy.cpp
    src/nucleotide.cpp
    src/search.cpp
    src/logger.cpp
    src/parameters/alignparameters.cpp
    src/indexinterface.cpp
)

# Add compression-specific or FM-index-specific source files
if(RUN_LENGTH_COMPRESSION)
    list(APPEND SOURCE_FILES src/bmove/bmove.cpp src/bmove/moverepr.cpp)
else()
    list(APPEND SOURCE_FILES src/fmindex/fmindex.cpp)
endif()

# Create the auxiliary library
add_library(AuxiliaryClasses ${SOURCE_FILES})

# ---- Update Git Submodules ----

option(COLUMBA_AUTO_UPDATE_SUBMODULES "Run git submodule update during CMake configure" OFF)

find_package(Git QUIET)
if(COLUMBA_AUTO_UPDATE_SUBMODULES)
    message(STATUS "Checking and updating submodules...")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_SUBMOD_RESULT
    )

    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "Failed to initialize git submodules. Please run 'git submodule update --init --recursive' manually.")
    endif()
else()
    message(STATUS "Using currently checked-out submodules (automatic updates disabled)")
endif()

# Update the Big-BWT submodule
if(RUN_LENGTH_COMPRESSION)
    find_package(Git QUIET)
    set(BIG_BWT_USABLE TRUE)

    # Find Python 3.8 or greater
    find_package(Python3 3.8 COMPONENTS Interpreter)

    if (Python3_FOUND)

        message(STATUS "Python executable detected by CMake: ${Python3_EXECUTABLE}")

        # Check Python version
        execute_process(
            COMMAND ${Python3_EXECUTABLE} -c "import sys; print(f'{sys.version_info[0]} {sys.version_info[1]}')"
            OUTPUT_VARIABLE PYTHON_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        # Split the major and minor version
        string(REPLACE " " ";" PYTHON_VERSION_LIST "${PYTHON_VERSION}")
        list(GET PYTHON_VERSION_LIST 0 PYTHON_MAJOR_VERSION)
        list(GET PYTHON_VERSION_LIST 1 PYTHON_MINOR_VERSION)

        # Check if the Python version is >= 3.8
        if (PYTHON_MAJOR_VERSION GREATER 3 OR (PYTHON_MAJOR_VERSION EQUAL 3 AND PYTHON_MINOR_VERSION GREATER_EQUAL 8))
            message(STATUS "Found Python: ${Python3_EXECUTABLE} (Version ${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION})")

            # Check if psutil is installed
            execute_process(
                COMMAND ${Python3_EXECUTABLE} -c "import psutil"
                RESULT_VARIABLE PSUTIL_INSTALLED
                ERROR_QUIET
                OUTPUT_QUIET
            )

            if (NOT PSUTIL_INSTALLED EQUAL "0")
                message(WARNING "Python package 'psutil' is required for prefix-free parsing, but not installed. Please install it using pip (pip install psutil) if you want to use prefix-free parsing.")
                set(BIG_BWT_USABLE FALSE)
            else()
                message(STATUS "Python 3.8 or greater and psutil are installed.")
                set(BIG_BWT_USABLE TRUE)
            endif()

        else()
            message(WARNING "Found Python version ${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION}, but Python 3.8 or greater is required for prefix-free parsing.")
            set(BIG_BWT_USABLE FALSE)
        endif()
    else()
        message(WARNING "Python not found, but is required for prefix-free parsing")
        set(BIG_BWT_USABLE FALSE)
    endif()

    if(BIG_BWT_USABLE)
        if(COLUMBA_AUTO_UPDATE_SUBMODULES AND GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
            message(STATUS "Updating Big-BWT submodule")
            execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init -q external/Big-BWT
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                            RESULT_VARIABLE GIT_SUBMOD_RESULT)
            if(NOT GIT_SUBMOD_RESULT EQUAL "0")
                message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
            endif()
        else()
            message(STATUS "Using local Big-BWT checkout")
        endif()
    else()
        message(WARNING "Prefix-free parsing is not usable due to missing dependencies (Python 3.8 and psutil), if you do not want to use prefix-free parsing you can ignore this message.")
    endif()
endif()

# ---- Executables and Linking ----

add_executable(columba src/parallel.cpp src/fastq.cpp src/seqfile.cpp src/util.cpp)
target_link_libraries(columba AuxiliaryClasses)

if(RUN_LENGTH_COMPRESSION)
    add_executable(columba_build src/logger.cpp src/buildindex.cpp src/bmove/moverepr.cpp src/parameters/buildparameters.cpp src/seqfile.cpp src/util.cpp)
    if(BIG_BWT_USABLE)
        add_definitions(-DBIG_BWT_USABLE)
        # Build the Big-BWT submodule if compression is enabled
        add_custom_command(
            OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT/bwtcheck
            COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT make
            COMMENT "Building Big-BWT executable"
        )
        add_custom_target(bigbwt_target ALL DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/external/Big-BWT/bwtcheck)
        add_dependencies(columba_build bigbwt_target)
    endif()
    # Copy the script to the build directory after columba_build is built
    add_custom_command(
        TARGET columba_build POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/src/bmove/columba_build_pfp.sh ${CMAKE_BINARY_DIR}/columba_build_pfp.sh
        COMMENT "Copying columba_build_pfp.sh to the build directory"
    )

else()
    add_executable(columba_build src/buildindex.cpp src/logger.cpp src/parameters/buildparameters.cpp src/seqfile.cpp )
endif()

# ---- External Dependencies ----

# overwrite the message function to suppress messages if MESSAGE_QUIET is set
function(message)
    if (NOT MESSAGE_QUIET)
        _message(${ARGN})
    endif()
endfunction()
set(MESSAGE_QUIET OFF)

set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "Suppress developer warnings")
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)

# FetchContent module for including external projects
include(FetchContent)
include(ExternalProject)

# libsais library (added to columba_build)

if(THIRTY_TWO)
    # Set specific build options for libsais
    message(STATUS "Fetching libsais library")
    set(LIBSAIS_USE_OPENMP OFF CACHE BOOL "Disable OpenMP support in libsais")
    set(LIBSAIS_BUILD_SHARED_LIB OFF CACHE BOOL "Build libsais as a static library")

    FetchContent_Declare(
        libsais
        GIT_REPOSITORY https://github.com/IlyaGrebnov/libsais.git
        GIT_TAG v2.8.5
    )
    FetchContent_MakeAvailable(libsais)
    target_link_libraries(columba_build libsais)
else()
    find_package(Git QUIET)

    if(COLUMBA_AUTO_UPDATE_SUBMODULES AND GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
        message(STATUS "Updating libdivsufsort submodule")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init external/libdivsufsort
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    else()
        message(STATUS "Using local libdivsufsort checkout")
    endif()

    set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build divsufsort64")
    set(BUILD_EXAMPLES OFF CACHE BOOL "Build examples")

    set(MESSAGE_QUIET ON)
    add_subdirectory(external/libdivsufsort)
    set(MESSAGE_QUIET OFF)

    # If needed, add dependencies to ensure libdivsufsort is built first
    add_dependencies(columba_build divsufsort64)

    # Link against the libdivsufsort library
    target_link_libraries(columba_build divsufsort64)

    # Specify include directories for the headers that are generated during build
    target_include_directories(columba_build PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/external/libdivsufsort/include
    )
endif()

# parallel-hashmap library
message(STATUS "Fetching parallel-hashmap library")
FetchContent_Declare(
    parallel-hashmap
    GIT_REPOSITORY https://github.com/greg7mdp/parallel-hashmap.git
    GIT_TAG e5b892baed478513adcb6425773cae1eda033057
)
FetchContent_MakeAvailable(parallel-hashmap)
include_directories(${parallel-hashmap_SOURCE_DIR})

# fmt library
message(STATUS "Fetching fmt library")
FetchContent_Declare(
    fmtlib
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG 10.2.1
)
FetchContent_MakeAvailable(fmtlib)
target_link_libraries(AuxiliaryClasses fmt::fmt)
target_link_libraries(columba fmt::fmt)
if(RUN_LENGTH_COMPRESSION)
    target_link_libraries(columba_build fmt::fmt)
endif()

# Enable IWYU (include-what-you-use) if desired
option(ENABLE_IWYU "Enable IWYU" OFF)
if(ENABLE_IWYU)
    find_program(IWYU_EXECUTABLE NAMES include-what-you-use iwyu)
    if(IWYU_EXECUTABLE)
        message(STATUS "Found IWYU: ${IWYU_EXECUTABLE}")
        set_property(TARGET AuxiliaryClasses PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${IWYU_EXECUTABLE})
        set_property(TARGET columba PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${IWYU_EXECUTABLE})
        set_property(TARGET columba_build PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${IWYU_EXECUTABLE})
    else()
        message(FATAL_ERROR "IWYU not found. Please install IWYU or disable ENABLE_IWYU.")
    endif()
endif()

# ---- Build Configuration ----

# Set default build type to Release if not specified
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# ---- Installed Dependencies and Linking ----

# Check for zlib
find_package(ZLIB)
if(ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIR})
    target_link_libraries(columba ${ZLIB_LIBRARY})
    target_link_libraries(columba_build ${ZLIB_LIBRARY})
endif()

# Link Threads package
find_package(Threads REQUIRED)
target_link_libraries(columba Threads::Threads)

# SDLS library (if RUN_LENGTH_COMPRESSION is enabled)
set(MY_SDSL_LIB "")
set(MY_SDSL_INC "")

if(RUN_LENGTH_COMPRESSION)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

    if(DEFINED SDSL_LIBRARY AND DEFINED SDSL_INCLUDE_DIR
       AND EXISTS "${SDSL_LIBRARY}" AND EXISTS "${SDSL_INCLUDE_DIR}")
        message(STATUS "Using caller-provided SDSL library: ${SDSL_LIBRARY}")
        message(STATUS "Using caller-provided SDSL include directory: ${SDSL_INCLUDE_DIR}")

        set(MY_SDSL_LIB "${SDSL_LIBRARY}")
        set(MY_SDSL_INC "${SDSL_INCLUDE_DIR}")

        include_directories(${MY_SDSL_INC})
        target_link_libraries(columba ${MY_SDSL_LIB})
        target_link_libraries(columba_build ${MY_SDSL_LIB})

    else()
        find_package(SDSL QUIET)

        if(SDSL_FOUND)
            message(STATUS "Found SDSL library at ${SDSL_LIBRARY}")
            message(STATUS "SDSL include directory: ${SDSL_INCLUDE_DIR}")

            include_directories(${SDSL_INCLUDE_DIR})
            target_link_libraries(columba ${SDSL_LIBRARY})
            target_link_libraries(columba_build ${SDSL_LIBRARY})

            set(MY_SDSL_LIB "${SDSL_LIBRARY}")
            set(MY_SDSL_INC "${SDSL_INCLUDE_DIR}")

        else()
            message(STATUS "SDSL-lite not found, fetching and installing it locally...")
            # existing fallback block unchanged
        endif()
    endif()
endif()
# ------------------------------------------------------------------------------
# HTSLIB - local pinned submodule at external/htslib (version 1.14)
# Built into the main build tree, then reused by PFP and LevioSAM
# ------------------------------------------------------------------------------

set(HTSLIB_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/htslib")
set(HTSLIB_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")

# Use source headers for configure-time consumers such as PFP
set(HTS_BUILD_INC "${HTSLIB_SOURCE_DIR}")

# Use installed headers for external projects configured later, such as LevioSAM
set(HTS_INSTALL_INC "${HTSLIB_INSTALL_PREFIX}/include")

set(HTS_LIB "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libhts.so")

if(NOT EXISTS "${HTSLIB_SOURCE_DIR}")
    message(FATAL_ERROR "HTSlib submodule not found at ${HTSLIB_SOURCE_DIR}")
endif()

find_library(HTS_LIB_FOUND hts PATHS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} NO_DEFAULT_PATH)

if(NOT HTS_LIB_FOUND)
    message(STATUS "htslib library not found. Building from submodule.")

    ExternalProject_Add(
        htslib_proj
        SOURCE_DIR ${HTSLIB_SOURCE_DIR}
        BINARY_DIR ${HTSLIB_SOURCE_DIR}
        UPDATE_COMMAND autoreconf -i
        CONFIGURE_COMMAND ./configure --disable-bz2 --prefix=${HTSLIB_INSTALL_PREFIX}
        BUILD_COMMAND $(MAKE)
        INSTALL_COMMAND $(MAKE) install
    )
else()
    message(STATUS "htslib library found at ${HTS_LIB_FOUND}.")
endif()

if (NOT TARGET htslib)
    message(STATUS "Creating htslib global library accessible to all CMakeLists.txt")
    add_library(htslib SHARED IMPORTED GLOBAL)
    set_property(TARGET htslib PROPERTY IMPORTED_LOCATION ${HTS_LIB})
    target_include_directories(htslib INTERFACE ${HTS_BUILD_INC})
endif()

target_include_directories(AuxiliaryClasses PUBLIC "${HTS_BUILD_INC}")
target_include_directories(AuxiliaryClasses PUBLIC "${HTS_BUILD_INC}/htslib")
target_link_libraries(AuxiliaryClasses htslib)
if(TARGET htslib_proj)
    add_dependencies(AuxiliaryClasses htslib_proj)
endif()
if(RUN_LENGTH_COMPRESSION AND MY_SDSL_INC)
    target_include_directories(AuxiliaryClasses PUBLIC "${MY_SDSL_INC}")
endif()

# ------------------------------------------------------------------------------
# PFP - local submodule, linked against the same HTSlib as above
# ------------------------------------------------------------------------------

set(LEVIOSAM_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/leviosam")
set(PFP_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/pfp")
if(EXISTS "${PFP_SOURCE_DIR}/CMakeLists.txt")
    message(STATUS "Adding PFP submodule")

    set(HTSlib_INCLUDE_DIRS "${HTS_BUILD_INC}" CACHE PATH "HTSlib include directory for PFP" FORCE)
    set(HTSlib_LIBRARIES "${HTS_LIB}" CACHE FILEPATH "HTSlib library for PFP" FORCE)
    set(HTSlib_INCLUDE_DIR "${HTS_BUILD_INC}" CACHE PATH "Legacy HTSlib include directory for PFP" FORCE)
    set(HTSlib_LIBRARY "${HTS_LIB}" CACHE FILEPATH "Legacy HTSlib library for PFP" FORCE)

    if(RUN_LENGTH_COMPRESSION AND EXISTS "${LEVIOSAM_SOURCE_DIR}/src/leviosam.hpp")
        message(STATUS "Providing local LevioSAM headers to PFP")
        set(PFP_USE_INTERNAL_LEVIOSAM OFF CACHE BOOL "Build/fetch LevioSAM inside PFP" FORCE)

        if(NOT TARGET lvsam)
            add_library(lvsam INTERFACE)
        endif()

        target_include_directories(lvsam INTERFACE "${LEVIOSAM_SOURCE_DIR}/src")
        target_link_libraries(lvsam INTERFACE htslib)

        if(MY_SDSL_INC)
            target_include_directories(lvsam INTERFACE "${MY_SDSL_INC}")
        endif()

        if(MY_SDSL_LIB)
            target_link_libraries(lvsam INTERFACE "${MY_SDSL_LIB}")
        endif()
    endif()
    
    set(_PFP_BUILD_SHARED_LIBS "${BUILD_SHARED_LIBS}")
    set(BUILD_SHARED_LIBS OFF)
    add_subdirectory("${PFP_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/external/pfp_build")
    set(BUILD_SHARED_LIBS "${_PFP_BUILD_SHARED_LIBS}")
    unset(_PFP_BUILD_SHARED_LIBS)

    foreach(_pfp_test IN ITEMS unit_tests unit_tests64)
        if(TARGET ${_pfp_test})
            set_target_properties(${_pfp_test} PROPERTIES EXCLUDE_FROM_ALL TRUE)
        endif()
    endforeach()

    set(PFP_TARGETS
        pfp
        pfp64
        pfp++
        pfp++64
        AuPair
        AuPair64
        mpfp++
        mpfp++64
        vcf_to_fa
        pfp++_unparse
        pfp++64_unparse
        check
        check64
    )

    foreach(_t IN LISTS PFP_TARGETS)
        if(TARGET ${_t})
            target_include_directories(${_t} PUBLIC "${HTS_BUILD_INC}/htslib")
            target_include_directories(${_t} PUBLIC "${HTS_BUILD_INC}")

            if(TARGET htslib_proj)
                add_dependencies(${_t} htslib_proj)
            endif()

            target_link_libraries(${_t} htslib)

            if(RUN_LENGTH_COMPRESSION AND MY_SDSL_LIB)
                target_include_directories(${_t} PUBLIC "${MY_SDSL_INC}")
                target_link_libraries(${_t} "${MY_SDSL_LIB}")
            endif()

            if(OpenMP_CXX_FOUND)
                target_link_libraries(${_t} OpenMP::OpenMP_CXX)
            endif()
        endif()
    endforeach()

    set(PFP_EXECUTABLES
        pfp++
        pfp++64
        AuPair
        AuPair64
        mpfp++
        mpfp++64
        vcf_to_fa
        pfp++_unparse
        pfp++64_unparse
        check
        check64
    )

    foreach(PFP_EXE IN LISTS PFP_EXECUTABLES)
        if(TARGET ${PFP_EXE})
            add_custom_target(copy_${PFP_EXE} ALL
                COMMAND ${CMAKE_COMMAND} -E copy
                    $<TARGET_FILE:${PFP_EXE}>
                    "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PFP_EXE}"
                DEPENDS ${PFP_EXE}
                COMMENT "Copying ${PFP_EXE} to main bin directory"
            )
        endif()
    endforeach()
else()
    message(WARNING "PFP source not found at ${PFP_SOURCE_DIR}; skipping.")
endif()

# ------------------------------------------------------------------------------
# LevioSAM - external build, forced to use the same HTSlib and SDSL
# ------------------------------------------------------------------------------

if(EXISTS "${LEVIOSAM_SOURCE_DIR}/CMakeLists.txt")
    set(LEVIOSAM_DEPS "")
    if(TARGET sdsl-lite_ext)
        list(APPEND LEVIOSAM_DEPS sdsl-lite_ext)
    endif()
    if(TARGET htslib_proj)
        list(APPEND LEVIOSAM_DEPS htslib_proj)
    endif()

    message(STATUS "Configuring LevioSAM against local HTSlib build")
    message(STATUS "LevioSAM will use HTSlib library: ${HTS_LIB}")
    message(STATUS "LevioSAM will use HTSlib includes: ${HTS_INSTALL_INC}")

    get_filename_component(MY_SDSL_LIB_DIR "${MY_SDSL_LIB}" DIRECTORY)

    _message(STATUS "DEBUG MY_SDSL_LIB      = [${MY_SDSL_LIB}]")
    _message(STATUS "DEBUG MY_SDSL_INC      = [${MY_SDSL_INC}]")
    _message(STATUS "DEBUG HTS_LIB          = [${HTS_LIB}]")
    _message(STATUS "DEBUG HTS_INSTALL_INC  = [${HTS_INSTALL_INC}]")
    
    ExternalProject_Add(leviosam_proj
        SOURCE_DIR "${LEVIOSAM_SOURCE_DIR}"
        PREFIX "${CMAKE_BINARY_DIR}/external/leviosam_wrapper"
        DEPENDS ${LEVIOSAM_DEPS}
        LIST_SEPARATOR "|"

        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/external/leviosam_install
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}

            "-DCMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}|${MY_SDSL_LIB_DIR}"
            "-DCMAKE_INCLUDE_PATH=${HTS_INSTALL_INC}|${MY_SDSL_INC}"

            -DHTS_LIB=${HTS_LIB}
            -DSDSL_LIB=${MY_SDSL_LIB}

        INSTALL_COMMAND make install
    )

    ExternalProject_Add_Step(leviosam_proj copy_binary
        COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_BINARY_DIR}/external/leviosam_install/bin/leviosam"
            "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/leviosam"
        DEPENDEES install
        COMMENT "Copying leviosam executable to project bin directory"
    )
else()
    message(WARNING "LevioSAM source not found at ${LEVIOSAM_SOURCE_DIR}; skipping.")
endif()

# ---- Installation ----

# Install targets to the bin directory
install(TARGETS columba columba_build DESTINATION bin)

if(EXISTS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/leviosam")
    install(PROGRAMS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/leviosam" DESTINATION bin)
endif()

foreach(PFP_EXE IN ITEMS
    pfp++ pfp++64 AuPair AuPair64 mpfp++ mpfp++64
    vcf_to_fa pfp++_unparse pfp++64_unparse check check64
)
    if(EXISTS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PFP_EXE}")
        install(PROGRAMS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PFP_EXE}" DESTINATION bin)
    endif()
endforeach()
