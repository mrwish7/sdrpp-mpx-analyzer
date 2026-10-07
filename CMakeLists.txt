cmake_minimum_required(VERSION 3.13)
project(mpx_analyzer)

file(GLOB SRC "src/*.cpp")

if (TARGET sdrpp_core)
    # Built from within the SDR++ build (add_subdirectory from SDR++'s CMakeLists.txt)
    include(${SDRPP_MODULE_CMAKE})
else ()
    # Standalone build against an SDR++ source tree.
    # The plugin must be built against the same SDR++ version/commit that will load it.
    set(SDRPP_SOURCE_DIR "" CACHE PATH "Path to the SDR++ source tree")
    if (NOT EXISTS "${SDRPP_SOURCE_DIR}/core/CMakeLists.txt")
        message(FATAL_ERROR "SDRPP_SOURCE_DIR must point to the SDR++ source tree, e.g. -DSDRPP_SOURCE_DIR=$HOME/SDRPlusPlus")
    endif ()

    # Options normally defined by the SDR++ root CMakeLists.txt and used by core/
    option(OPT_BACKEND_GLFW "Use the GLFW backend" ON)
    option(OPT_BACKEND_ANDROID "Use the Android backend" OFF)
    option(OPT_OVERRIDE_STD_FILESYSTEM "Use a local version of std::filesystem on systems that don't have it yet" OFF)
    option(USE_INTERNAL_LIBCORRECT "Use an internal version of libcorrect" ON)
    option(USE_BUNDLE_DEFAULTS "Set the default resource and module directories to the right ones for a MacOS .app" OFF)

    if (NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    endif ()

    # Same compiler flags as the SDR++ root CMakeLists.txt
    if (CMAKE_BUILD_TYPE MATCHES "Debug")
        if (MSVC)
            set(SDRPP_COMPILER_FLAGS /std:c++17 /EHsc)
        elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            set(SDRPP_COMPILER_FLAGS -g -Og -std=c++17 -Wno-unused-command-line-argument -undefined dynamic_lookup)
        else ()
            set(SDRPP_COMPILER_FLAGS -g -Og -std=c++17)
        endif ()
    else ()
        if (MSVC)
            set(SDRPP_COMPILER_FLAGS /O2 /Ob2 /std:c++17 /EHsc)
        elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            set(SDRPP_COMPILER_FLAGS -O3 -std=c++17 -Wno-unused-command-line-argument -undefined dynamic_lookup)
        else ()
            set(SDRPP_COMPILER_FLAGS -O3 -std=c++17)
        endif ()
    endif ()
    set(SDRPP_MODULE_COMPILER_FLAGS ${SDRPP_COMPILER_FLAGS})

    # Same default install prefix as SDR++ so the plugin lands in its plugin directory
    if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
            set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "..." FORCE)
        else ()
            set(CMAKE_INSTALL_PREFIX "/usr" CACHE PATH "..." FORCE)
        endif ()
    endif ()

    # Build sdrpp_core only to compile and link against; it is not installed
    add_subdirectory("${SDRPP_SOURCE_DIR}/core" "${CMAKE_CURRENT_BINARY_DIR}/sdrpp_core" EXCLUDE_FROM_ALL)

    # Equivalent of SDR++'s sdrpp_module.cmake
    add_library(${PROJECT_NAME} SHARED ${SRC})
    target_link_libraries(${PROJECT_NAME} PRIVATE sdrpp_core)
    target_compile_options(${PROJECT_NAME} PRIVATE ${SDRPP_MODULE_COMPILER_FLAGS})
    set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
    install(TARGETS ${PROJECT_NAME} DESTINATION lib/sdrpp/plugins)
endif ()
