cmake_minimum_required(VERSION 3.10.0)

project(geometric-flow)



### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


FIND_PACKAGE(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
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
  src/Beads.cpp
  src/Mem-3dg.cpp
  src/Mem-3dg_implicit.cpp
  ../../utils/src/colormap.cpp
)



set(SRCScluster 
  src/main_cluster.cpp  
  src/Mem-3dg.cpp
  src/Beads.cpp
  src/Energy_Handler.cpp
  src/Interaction.cpp
  ../../utils/src/colormap.cpp
  src/ShapeNLP.cpp
)

set(SRCSvisualize 
  src/main_visualize.cpp  
  src/Mem-3dg.cpp
  src/Beads.cpp
  src/Energy_Handler.cpp
  src/Interaction.cpp
  ../../utils/src/colormap.cpp
)

set(SRCSmainpreserv 
  src/main_preserving_mean_curv.cpp  
  src/Mem-3dg.cpp
  src/Beads.cpp
  src/Interaction.cpp
  ../../utils/src/colormap.cpp
)

# set(SRCSclusterbeads 
#  src/main_cluster_bead.cpp
#  src/Mem-3dg.cpp
#  src/Beads.cpp
#  ../../utils/src/colormap.cpp
# )

# set(SRCSclusterpulling 
#   src/main_cluster_bead_pulling.cpp  
#   src/Mem-3dg.cpp
#   src/Beads.cpp
#   ../../utils/src/colormap.cpp
# )


# set(SRCSclusterpulling2 
#   src/main_pulling_two_beads.cpp  
#   src/Mem-3dg.cpp
#   src/Beads.cpp
#   ../../utils/src/colormap.cpp
# )

set(SRCSprocessing 
  src/Post_processing.cpp  
  src/Mem-3dg.cpp
  src/Beads.cpp
  src/Interaction.cpp
  ../../utils/src/colormap.cpp
)

set(SRCSprocessingtube 
  src/Processing_tube.cpp  
  src/Mem-3dg.cpp
  src/Beads.cpp
  src/Interaction.cpp
  ../../utils/src/colormap.cpp
)

# set(SRCSclusterforcefield 
#  src/main_force_field.cpp
#  src/Mem-3dg.cpp
#  src/Beads.cpp
#  ../../utils/src/colormap.cpp
# )
# set(WindowSRCS 
#   src/Window.cpp

#   src/Mem-3dg.cpp
#   ../../utils/src/colormap.cpp
# )
# set(TEST_Grad_SRCS 
#   src/E_tests.cpp
#   # src/mean-curvature-flow.cpp
#   src/Mem-3dg.cpp
#   src/Beads.cpp
#   ../../utils/src/colormap.cpp
# )

set(TEST_hessian_SRCS
  src/main_test_hess.cpp
  src/Mem-3dg.cpp
  src/Beads.cpp
  src/Energy_Handler.cpp
  src/Interaction.cpp
  ../../utils/src/colormap.cpp
)

set(TEST_Grad_th_SRCS 
  src/E_test_2.cpp
  src/Mem-3dg.cpp
  src/Beads.cpp
  src/Interaction.cpp
  ../../utils/src/colormap.cpp
)

set(TEST_fan_SRCS 
  src/Test_fan_basic.cpp
  src/Mem-3dg.cpp
  src/Beads.cpp
  src/Interaction.cpp
  ../../utils/src/colormap.cpp
)


# set(TEST_SRCS
#   src/mean-curvature-flow.cpp
#   src/modified-mean-curvcature-flow.cpp
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
add_subdirectory(../../deps/libarcsim deps/libarcsim)




find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -I../../deps/include -I/opt/local/include/")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()





# # Change the executable name by changing "main" in the lines below.
add_executable(main "${SRCS}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main)
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps")
# add the args.hxx project which we use for command line args
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")

# find_package(EigenRand REQUIRED)

target_link_libraries(main PUBLIC geometry-central polyscope arcsim)


add_executable(main_cluster "${SRCScluster}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main_cluster)
target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")

target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps")

# target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/EigenRand/")
# add the args.hxx project which we use for command line args
target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")

target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/include/")
target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/")
target_link_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/lib/")

# add_library(ipopt_interface ipopt_interface.cpp)
# target_link_libraries(main_cluster PUBLIC ipopt_interface)
# target_link_libraries(main_cluster PUBLIC ipopt)


add_library(EigenRand INTERFACE)
target_include_directories(EigenRand INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/EigenRand-0.5.1/")

link_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../deps/lib")
link_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../deps/include")

target_link_libraries(main_cluster PUBLIC geometry-central polyscope arcsim EigenRand ipopt)


add_executable(main_visualize "${SRCSvisualize}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main_visualize)
target_include_directories(main_visualize PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(main_visualize PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")

target_include_directories(main_visualize PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
target_include_directories(main_visualize PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps")

# add the args.hxx project which we use for command line args
target_include_directories(main_visualize PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")

target_include_directories(main_visualize PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")

target_link_libraries(main_visualize PUBLIC geometry-central polyscope arcsim)


# add_executable(main_cluster "${SRCScluster}")
# set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main_cluster)
# target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")
# # add the args.hxx project which we use for command line args
# target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
# target_include_directories(main_cluster PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")
# target_link_libraries(main_cluster PUBLIC geometry-central polyscope arcsim)

add_executable(main_preserving_mean_curv "${SRCSmainpreserv}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main_preserving_mean_curv)
target_include_directories(main_preserving_mean_curv PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(main_preserving_mean_curv PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# add the args.hxx project which we use for command line args
target_include_directories(main_preserving_mean_curv PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_link_libraries(main_preserving_mean_curv PUBLIC geometry-central polyscope)


# add_executable(main_cluster_beads "${SRCSclusterbeads}")
# set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main_cluster_beads)
# target_include_directories(main_cluster_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(main_cluster_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# target_include_directories(main_cluster_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")
# target_include_directories(main_cluster_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps")
# # add the args.hxx project which we use for command line args
# target_include_directories(main_cluster_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
# target_include_directories(main_cluster_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")
# target_link_libraries(main_cluster_beads PUBLIC geometry-central polyscope arcsim)

# add_executable(main_cluster_pulling "${SRCSclusterpulling}")
# set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main_cluster_pulling)
# target_include_directories(main_cluster_pulling PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(main_cluster_pulling PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")

# target_include_directories(main_cluster_pulling PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# target_include_directories(main_cluster_pulling PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps")
# target_include_directories(main_cluster_pulling PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")

# target_include_directories(main_cluster_pulling PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")

# target_link_libraries(main_cluster_pulling PUBLIC geometry-central polyscope arcsim)


# add_executable(main_cluster_pulling_beads "${SRCSclusterpulling2}")
# set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main_cluster_pulling_beads)
# target_include_directories(main_cluster_pulling_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(main_cluster_pulling_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")

# target_include_directories(main_cluster_pulling_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# target_include_directories(main_cluster_pulling_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps")

# # add the args.hxx project which we use for command line args
# target_include_directories(main_cluster_pulling_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")

# target_include_directories(main_cluster_pulling_beads PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")

# target_link_libraries(main_cluster_pulling_beads PUBLIC geometry-central polyscope arcsim)



# add_executable(main_cluster_field "${SRCSclusterforcefield}")
# set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main_cluster_field)
# target_include_directories(main_cluster_field PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(main_cluster_field PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# target_include_directories(main_cluster_field PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")
# target_include_directories(main_cluster_field PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps")
# # add the args.hxx project which we use for command line args
# target_include_directories(main_cluster_field PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
# target_include_directories(main_cluster_field PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")
# target_link_libraries(main_cluster_field PUBLIC geometry-central polyscope arcsim)



add_executable(Post_processing "${SRCSprocessing}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Post_processing)
target_include_directories(Post_processing PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(Post_processing PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")

target_include_directories(Post_processing PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
target_include_directories(Post_processing PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps")

# add the args.hxx project which we use for command line args
target_include_directories(Post_processing PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")

# target_include_directories(Post_processing PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")

target_link_libraries(Post_processing PUBLIC geometry-central polyscope)



add_executable(Processing_tube "${SRCSprocessingtube}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Processing_tube)
target_include_directories(Processing_tube PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(Post_processing PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")

target_include_directories(Processing_tube PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
target_include_directories(Processing_tube PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps")

# add the args.hxx project which we use for command line args
target_include_directories(Processing_tube PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")

# target_include_directories(Post_processing PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")

target_link_libraries(Processing_tube PUBLIC geometry-central polyscope)

















# add_executable(mainbunny "${SRCSfun}")
# set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT main)
# target_include_directories(mainbunny PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(mainbunny PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# # add the args.hxx project which we use for command line args
# target_include_directories(mainbunny PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
# target_link_libraries(mainbunny geometry-central polyscope)

# add_executable(Window "${WindowSRCS}")
# set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Window)
# target_include_directories(Window PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(Window PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# # add the args.hxx project which we use for command line args
# target_include_directories(Window PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
# target_link_libraries(Window geometry-central polyscope)


# add_executable(Grad_tests "${TEST_Grad_SRCS}")
# set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Grad_tests)
# target_include_directories(Grad_tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(Grad_tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# # add the args.hxx project which we use for command line args
# target_include_directories(Grad_tests PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
# target_link_libraries(Grad_tests geometry-central polyscope)


add_executable(Grad_tests_th "${TEST_Grad_th_SRCS}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Grad_tests_th)
target_include_directories(Grad_tests_th PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(Grad_tests_th PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
target_include_directories(Grad_tests_th PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")

target_include_directories(Grad_tests_th PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps")
# add the args.hxx project which we use for command line args
target_include_directories(Grad_tests_th PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_include_directories(Grad_tests_th PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")

target_link_libraries(Grad_tests_th PUBLIC geometry-central polyscope arcsim)


add_executable(Test_hessian "${TEST_hessian_SRCS}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Test_hessian)
target_include_directories(Test_hessian PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(Test_hessian PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
target_include_directories(Test_hessian PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src_arcsim/")

target_include_directories(Test_hessian PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps")
# add the args.hxx project which we use for command line args
target_include_directories(Test_hessian PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
target_include_directories(Test_hessian PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/libarcsim/include")

target_link_libraries(Test_hessian PUBLIC geometry-central polyscope arcsim)


add_executable(Test_fan "${TEST_fan_SRCS}")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Test_fan)
target_include_directories(Test_fan PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")

target_include_directories(Test_fan PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../utils/include")
# add the args.hxx project which we use for command line args
target_include_directories(Test_fan PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")


target_link_libraries(Test_fan geometry-central polyscope)


# TEST_fan_SRCS

# add_executable(test-flow "${TEST_SRCS}")
# target_include_directories(test-flow PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_link_libraries(test-flow geometry-central polyscope gtest)























  # src/mean-curvature-flow.cpp
  # src/modified-mean-curvature-flow.cpp
  # src/normal-flow.cpp
  # src/gauss-curvature-flow.cpp
  # src/willmore-flow.cpp
  # src/willmore-flow-2.cpp
  # src/willmore-flow-scho.cpp
  # src/membrane-flow.cpp