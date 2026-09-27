pybind11_add_module(${PROJECT_NAME}
    ${CMAKE_CURRENT_LIST_DIR}/PostProcessorWrapper_py.cpp
)

target_link_libraries(
    ${PROJECT_NAME}
    PRIVATE
        Eigen3::Eigen
)

target_compile_definitions(${PROJECT_NAME}
    PRIVATE BUILD_PYTHON_BINDINGS
)

# -----------------------------------------------------------------------------
# Link Boost libraries
# Modern Boost (1.70+) is mostly header-only, but filesystem, timer, etc. still
# need linking
if(TARGET Boost::headers)
  # Modern Boost CONFIG mode (1.70+)
  # Need to manually find and link non-header-only libraries
  target_link_libraries(${PROJECT_NAME} PRIVATE Boost::headers)

  # Find Boost library directory - search in system paths and Boost_LIBRARY_DIRS
  find_library(BOOST_FILESYSTEM_LIB
    NAMES boost_filesystem libboost_filesystem
    PATHS ${Boost_LIBRARY_DIRS} /usr/lib64 /usr/lib /usr/local/lib
    PATH_SUFFIXES boost169
  )
  find_library(BOOST_SYSTEM_LIB
    NAMES boost_system libboost_system
    PATHS ${Boost_LIBRARY_DIRS} /usr/lib64 /usr/lib /usr/local/lib
    PATH_SUFFIXES boost169
  )

  if(BOOST_FILESYSTEM_LIB)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${BOOST_FILESYSTEM_LIB})
    message(STATUS "Linked Boost filesystem: ${BOOST_FILESYSTEM_LIB}")
  endif()
  if(BOOST_SYSTEM_LIB)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${BOOST_SYSTEM_LIB})
    message(STATUS "Linked Boost system: ${BOOST_SYSTEM_LIB}")
  endif()
elseif(TARGET Boost::boost)
  # Alternative modern Boost target
  target_link_libraries(${PROJECT_NAME} PRIVATE Boost::boost)

  # Still need to link non-header-only libraries
  find_library(BOOST_FILESYSTEM_LIB
    NAMES boost_filesystem libboost_filesystem
    PATHS ${Boost_LIBRARY_DIRS} /usr/lib64 /usr/lib /usr/local/lib
    PATH_SUFFIXES boost169
  )
  find_library(BOOST_SYSTEM_LIB
    NAMES boost_system libboost_system
    PATHS ${Boost_LIBRARY_DIRS} /usr/lib64 /usr/lib /usr/local/lib
    PATH_SUFFIXES boost169
  )

  if(BOOST_FILESYSTEM_LIB)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${BOOST_FILESYSTEM_LIB})
    message(STATUS "Linked Boost filesystem: ${BOOST_FILESYSTEM_LIB}")
  endif()
  if(BOOST_SYSTEM_LIB)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${BOOST_SYSTEM_LIB})
    message(STATUS "Linked Boost system: ${BOOST_SYSTEM_LIB}")
  endif()
else()
  # Legacy Boost with components (pre-1.70 or FindBoost module)
  # This is used on Linux manylinux with boost169
  target_link_libraries(
        ${PROJECT_NAME}
        PRIVATE
            Boost::date_time
            Boost::filesystem
            Boost::program_options
            Boost::system
            Boost::regex
            Boost::timer
            Boost::chrono
    )
  message(STATUS "Using legacy Boost component targets")
endif()

target_include_directories(
    ${PROJECT_NAME}
    PRIVATE
        ${OPENSMOKE_INCLUDE}
        ${CMAKE_CURRENT_LIST_DIR}/source
)

# Compiler-specific warning suppression
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(${PROJECT_NAME}
        PRIVATE
            -Wno-unused-value
            -Wno-mismatched-new-delete
            -Wno-switch
    )
elseif(MSVC)
  target_compile_options(${PROJECT_NAME}
        PRIVATE
            /wd4018  # signed/unsigned mismatch
            /wd4267  # conversion warnings
    )
endif()

target_compile_definitions(${PROJECT_NAME}
    PRIVATE
        VERSION_INFO="${PROJECT_VERSION}"
)

# --------------------------------------------------------------------
# Installation
# --------------------------------------------------------------------
install(
    TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION .
    RUNTIME DESTINATION .
)

# -----------------------------------------------------------------------------
# Standalone testfiles/*.cpp drivers (PYSMOKE_BUILD_TESTFILES, off by default)
# -----------------------------------------------------------------------------
# Bypass pybind entirely (no BUILD_PYTHON_BINDINGS define, so
# PostProcessorWrapper_py.h's pybind11 includes are skipped) - the fast
# iteration loop for the C++ backend used while refactoring it. Fixture data
# path is baked in via PYSMOKE_EXAMPLES_DATA_DIR instead of the hardcoded
# per-machine absolute paths the drivers used to carry, so they work from any
# checkout regardless of the current working directory.
if(PYSMOKE_BUILD_TESTFILES)
  set(PYSMOKE_TESTFILES
      ropa_homo
      ropa_het
      sensitivity_homo
      sensitivity_het
  )

  foreach(_test_name ${PYSMOKE_TESTFILES})
    add_executable(${_test_name} ${CMAKE_CURRENT_LIST_DIR}/testfiles/${_test_name}.cpp)

    target_include_directories(${_test_name}
        PRIVATE
            ${OPENSMOKE_INCLUDE}
            ${CMAKE_CURRENT_LIST_DIR}/source
    )

    target_compile_definitions(${_test_name}
        PRIVATE PYSMOKE_EXAMPLES_DATA_DIR="${CMAKE_SOURCE_DIR}/examples/data"
    )

    target_link_libraries(${_test_name} PRIVATE Eigen3::Eigen)

    if(TARGET Boost::headers)
        target_link_libraries(${_test_name} PRIVATE Boost::headers)
        if(BOOST_FILESYSTEM_LIB)
            target_link_libraries(${_test_name} PRIVATE ${BOOST_FILESYSTEM_LIB})
        endif()
        if(BOOST_SYSTEM_LIB)
            target_link_libraries(${_test_name} PRIVATE ${BOOST_SYSTEM_LIB})
        endif()
    elseif(TARGET Boost::boost)
        target_link_libraries(${_test_name} PRIVATE Boost::boost)
        if(BOOST_FILESYSTEM_LIB)
            target_link_libraries(${_test_name} PRIVATE ${BOOST_FILESYSTEM_LIB})
        endif()
        if(BOOST_SYSTEM_LIB)
            target_link_libraries(${_test_name} PRIVATE ${BOOST_SYSTEM_LIB})
        endif()
    else()
        target_link_libraries(${_test_name}
            PRIVATE
                Boost::date_time
                Boost::filesystem
                Boost::program_options
                Boost::system
                Boost::regex
                Boost::timer
                Boost::chrono
        )
    endif()
  endforeach()
endif()