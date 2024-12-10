cmake_minimum_required(VERSION 3.10.0)

project(geometric-flow)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Print the build type
# set(CMAKE_BUILD_TYPE "Debug")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

### Configure the compiler
# This is a basic, decent setup that should do something sane on most compilers

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++11 -Wall -Wextra")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3 -fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()
  
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  
  # https://gitlab.kitware.com/cmake/cmake/-/issues/22484
  include(CheckCXXCompilerFlag)
  unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
  CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")
  else()
    SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -mcpu=apple-m1 -DNDEBUG") # Apple M1
  endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

# Lets add the openMP







# Edit files to build here.
set(SRCS 
  src/main.cpp
  src/Mem-3dg.cpp
  src/mem-3dg_implicit.cpp
  ../../utils/src/colormap.cpp
)

set(SRCSparallel 
  src/main_parallel.cpp

  src/Mem-3dg_parallel.cpp
  ../../utils/src/colormap.cpp
)

set(SRCScluster 
  src/main_cluster.cpp
  src/Mem-3dg.cpp
  ../../utils/src/colormap.cpp
)

set(SRCSclusterparallel 
  src/main_cluster_parallel.cpp

  src/Mem-3dg_parallel.cpp
  ../../utils/src/colormap.cpp
)



set(WindowSRCS 
  src/Window.cpp
  src/mean-curvature-flow.cpp
  src/Mem-3dg.cpp
  src/normal-flow.cpp
  src/mem-3dg_implicit.cpp
  ../../utils/src/colormap.cpp
)
set(TEST_Grad_SRCS 
  src/E_tests.cpp
  # src/mean-curvature-flow.cpp
  src/Mem-3dg.cpp
  ../../utils/src/colormap.cpp
)

# set(simple_SRCS
#   src/Simple_script.cpp
#   # src/mean-curvature-flow.cpp
#   src/Mem-3dg.cpp
#   ../../utils/src/colormap.cpp
# )


# set(TEST_SRCS
#   src/mean-curvature-flow.cpp
#   src/modified-mean-curvature-flow.cpp
#   src/normal-flow.cpp
#   src/gauss-curvature-flow.cpp
#   src/willmore-flow.cpp
#   src/willmore-flow-2.cpp
#   src/willmore-flow-scho.cpp
#   src/membrane-flow.cpp
#   src/test-flow.cpp
# )

add_subdirectory(../../deps/geometry-central deps/geometry-central)
add_subdirectory(../../deps/polyscope deps/polyscope)
add_subdirectory(../../deps/googletest deps/googletest)


find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()




# Change the executable name by changing "main" in the lines below.
add_executable(main "${SRCS}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main)
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# add the args.hxx project which we use for command line args
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_link_libraries(main PUBLIC geometry-central polyscope)
# find_package(OpenMP)
# if(OPENMP_CXX_FOUND)
#   target_link_libraries(main PUBLIC )
# endif()


add_executable(main_cluster "${SRCScluster}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main_cluster)
target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# add the args.hxx project which we use for command line args
target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_link_libraries(main_cluster PUBLIC geometry-central polyscope)


add_executable(main_parallel "${SRCSparallel}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main_parallel)
target_include_directories(main_parallel PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(main_parallel PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# add the args.hxx project which we use for command line args
target_include_directories(main_parallel PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_link_libraries(main_parallel PUBLIC geometry-central polyscope OpenMP::OpenMP_CXX)



add_executable(main_cluster_parallel "${SRCSclusterparallel}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main_cluster_parallel)
target_include_directories(main_cluster_parallel PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(main_cluster_parallel PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# add the args.hxx project which we use for command line args
target_include_directories(main_cluster_parallel PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_link_libraries(main_cluster_parallel PUBLIC geometry-central polyscope OpenMP::OpenMP_CXX)



# add_executable(mainbunny "${SRCSfun}")
# set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main)
# target_include_directories(mainbunny PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(mainbunny PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# # add the args.hxx project which we use for command line args
# target_include_directories(mainbunny PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
# target_link_libraries(mainbunny geometry-central polyscope)

add_executable(Window "${WindowSRCS}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Window)
target_include_directories(Window PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(Window PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# add the args.hxx project which we use for command line args
target_include_directories(Window PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_link_libraries(Window geometry-central polyscope)


add_executable(Grad_tests "${TEST_Grad_SRCS}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Grad_tests)
target_include_directories(Grad_tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(Grad_tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# add the args.hxx project which we use for command line args
target_include_directories(Grad_tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_link_libraries(Grad_tests geometry-central polyscope)




# add_executable(test-flow "${TEST_SRCS}")
# target_include_directories(test-flow PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_link_libraries(test-flow geometry-central polyscope gtest)







# add_executable(Simple_script "${simple_SRCS}")
# set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Simple_script)
# target_include_directories(Simple_script PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(Simple_script PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# # add the args.hxx project which we use for command line args
# target_include_directories(Simple_script PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
# target_link_libraries(Simple_script PUBLIC geometry-central polyscope OpenMP::OpenMP_CXX)
















  # src/mean-curvature-flow.cpp
  # src/modified-mean-curvature-flow.cpp
  # src/normal-flow.cpp
  # src/gauss-curvature-flow.cpp
  # src/willmore-flow.cpp
  # src/willmore-flow-2.cpp
  # src/willmore-flow-scho.cpp
  # src/membrane-flow.cpp