# Install script for directory: /home/josh/Desktop/563pba/FluidSim/CIS563_SmokeBaseCode/cmake-build-debug/partio-src/src/lib

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/josh/Desktop/563pba/FluidSim/CIS563_SmokeBaseCode/Linux-4.10.0-x86_64")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/josh/Desktop/563pba/FluidSim/CIS563_SmokeBaseCode/cmake-build-debug/partio-build/lib/libpartio.a")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/josh/Desktop/563pba/FluidSim/CIS563_SmokeBaseCode/cmake-build-debug/partio-src/src/lib/Partio.h"
    "/home/josh/Desktop/563pba/FluidSim/CIS563_SmokeBaseCode/cmake-build-debug/partio-src/src/lib/PartioAttribute.h"
    "/home/josh/Desktop/563pba/FluidSim/CIS563_SmokeBaseCode/cmake-build-debug/partio-src/src/lib/PartioIterator.h"
    )
endif()

