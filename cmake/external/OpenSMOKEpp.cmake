# -----------------------------------------------------------------------------
# Configuration for OpenSMOKEpp
# -----------------------------------------------------------------------------
message(STATUS "Configuring OpenSMOKEpp...")

# First, check if OpenSMOKEpp is vendored in the repository
set(VENDORED_OPENSMOKE_PATH "${CMAKE_SOURCE_DIR}/external/opensmoke/source")

if(EXISTS "$ENV{OpenSMOKEpp_ROOT}")
    message(STATUS "Using system-installed OpenSMOKEpp from: $ENV{OpenSMOKEpp_ROOT}/source")
    set(OPENSMOKE_INCLUDE "$ENV{OpenSMOKEpp_ROOT}/source" CACHE PATH "Path to OpenSMOKE library")
elseif(EXISTS ${VENDORED_OPENSMOKE_PATH})
    message(STATUS "Using vendored OpenSMOKEpp from: ${VENDORED_OPENSMOKE_PATH}")
    set(OPENSMOKE_INCLUDE ${VENDORED_OPENSMOKE_PATH} CACHE PATH "Path to OpenSMOKE library")
else()
    # Fallback: Try to fetch from GitHub (requires authentication)
    message(STATUS "Vendored OpenSMOKEpp not found, attempting to fetch from GitHub...")
    message(WARNING "Note: OpenSMOKEpp repository may require authentication!")

    include(FetchContent)

    FetchContent_Declare(
      opensmoke
      GIT_REPOSITORY https://github.com/acuoci/OpenSMOKEpp.git
      GIT_TAG v0.24.0
      GIT_SHALLOW TRUE
      GIT_PROGRESS TRUE
      SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/_deps/opensmoke-src
    )

    FetchContent_MakeAvailable(opensmoke)

    set(OPENSMOKE_INCLUDE
      ${opensmoke_SOURCE_DIR}/source
      CACHE PATH "Path to OpenSMOKE library"
    )
endif()

message(STATUS "OpenSMOKEpp include path: ${OPENSMOKE_INCLUDE}")
