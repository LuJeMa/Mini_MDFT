#
# Some help to write CMakeLists.cmake for fortran project:
#
# - http://www.cmake.org/Wiki/CMakeForFortranExample
# - http://www.ngssc.se/courses/specialized-courses/advanced-programming/make_cmake.pdf
# - https://github.com/SethMMorton/cmake_fortran_template
#

# require cmake 2.8.0
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0 FATAL_ERROR)

##################
project("mini-mdft")
##################

# Fortran is disabled by default, so we need to wake-up
# note: string given to enable_language is case-sensitive
enable_language(Fortran)

#option to build FFTW-3.3.4 with MDFT
option(BUILD_FFTW "build FFTW-3.3.4 with MDFT" OFF)

## guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. Follow instructions in README.md.")
endif()

# Build type :
# - RELEASE (with optimization flags, without openmp)
# - DEBUG   (with debug flags, without openmp)
# make sure that the default is RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# allow cmake to use custom modules located in cmake subdir
# (e.g. FindFFTW)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
#find_package(MPI REQUIRED)

#
# check compiler (gfortran, ifort, ...)
# set FFLAGS depend on the compiler
# set where fortran module files will be located
#
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
get_filename_component (CC_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
#SET (CMAKE_Fortran_COMPILER  gfortran)
message(STATUS "Using compiler ${Fortran_COMPILER_NAME}")
message(STATUS "Using compiler ${CC_COMPILER_NAME}")

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/mod")

if (Fortran_COMPILER_NAME MATCHES "gfortran*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-march=native -ffree-line-length-none -O3 -ffast-math -pedantic -std=f2008 -fopenmp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-march=native -ffree-line-length-none -O0             -pedantic -std=f2008 -ffpe-trap=zero,underflow,overflow -Wall -fcheck=bounds -pg -fopenmp")
  #-ftree-vectorizer-verbose=10 -fopt-info
  #set (CMAKE_Fortran_FLAGS_DEBUG   "-g -fbacktrace -pedantic -fwhole-file -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -fbounds-check -pg -frecursive -fcheck=all -Wall -ffpe-trap=zero,underflow,overflow")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -xhost" )
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "pgf90")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif()

if (DEFINED Custom_Fortran_FLAGS_RELEASE)
  set (CMAKE_Fortran_FLAGS_RELEASE "${Custom_Fortran_FLAGS_RELEASE}")
endif()

message(STATUS "CMAKE_Fortran_FLAGS set to ${CMAKE_Fortran_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE set to ${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG set to ${CMAKE_Fortran_FLAGS_DEBUG}")


#add_subdirectory(vendor/bpp bpp EXCLUDE_FROM_ALL)
#
#bpp_preprocess(SRC_Bpp src/module_energy_cproj_mrso.F90.bpp src/module_rotation.F90.bpp)



#############################
# Now build mdft executable #
#############################

# set source files list
set(mdft_SRCS
  src/module_blas_etc.f90
  src/module_lecture.f90
  src/module_read_c_luc.f90
  src/module_rotation.f90
  src/module_energy_cproj_mrso.f90
  src/module_thermo.f90
  src/module_proc_ptrs.f90
  src/module_energy_and_gradient.f90
  src/energy_minimization.f90
  src/main.f90
  src/module_constants.f90
  src/module_density.f90
  src/module_fft.f90
  src/module_grid.f90
  src/module_init_simu.f90
  src/module_input.f90
  src/module_mathematica.f90
  src/module_precision_kinds.f90
  src/module_solute.f90
  src/module_solvent.f90
  src/module_system.f90
  src/module_time.f90
  src/module_wigner_d.f90
  src/module_orientation_projection_transform.f90
#  src/module_mpi.F90
)
#
# Create a symlink to the input and data folders
#
EXECUTE_PROCESS(COMMAND "${CMAKE_COMMAND}" -E create_symlink "${CMAKE_SOURCE_DIR}/data"  "${CMAKE_BINARY_DIR}/data")

#
# Build FFTW
#
include(ExternalProject)


if(NOT BUILD_FFTW)
	find_package(FFTW QUIET)
    if(NOT FFTW_FOUND)
		message(STATUS "WARNING: FFTW not found! FFTW-3.3.4 will be compiled from sources before MDFT.")
		set(BUILD_FFTW ON)
    endif()
endif()


if(BUILD_FFTW)
	ExternalProject_Add(
	  fftw
	  URL ${PROJECT_SOURCE_DIR}/src/fftw-3.3.4.tar.gz
	  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw3-cmake
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/libs/fftw3 CC=${CC_COMPILER_NAME}  --enable-threads  --enable-sse2 --enable-avx
	  BUILD_COMMAND make -j4
	  INSTALL_COMMAND make install -j4
	)
	ExternalProject_Add(
	  fftwf
	  URL ${PROJECT_SOURCE_DIR}/src/fftw-3.3.4.tar.gz
	  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw3-cmake
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${CMAKE_CURRENT_BINARY_DIR}/libs/fftw3 CC=${CC_COMPILER_NAME} --enable-float --enable-threads --enable-sse  --enable-sse2 --enable-avx
	  BUILD_COMMAND make -j4
	  INSTALL_COMMAND make install -j4
	)
	SET(FFTW_INCLUDES ${CMAKE_CURRENT_BINARY_DIR}/libs/fftw3/include)
	SET(FFTW_LIBRARIES ${CMAKE_CURRENT_BINARY_DIR}/libs/fftw3/lib/libfftw3.a;${CMAKE_CURRENT_BINARY_DIR}/libs/fftw3/lib/libfftw3f.a)
	include_directories(${FFTW_INCLUDES})
endif()

#get git commit
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
configure_file("${PROJECT_SOURCE_DIR}/src/module_git.f90.in" "git_module/module_git.f90")
list(APPEND mdft_SRCS "git_module/module_git.f90")

if (UNIX)
  message(STATUS "Building on a UNIX-like platform")


  include_directories(${FFTW_INCLUDES})
  add_executable(mini-mdft ${mdft_SRCS})

  message(STATUS "Compiling with fftw ${FFTW_LIBRARIES}")
  target_link_libraries(mini-mdft ${FFTW_LIBRARIES})

  if(BUILD_FFTW) 
    add_dependencies(mini-mdft fftw fftwf)
  endif()

#  if("${MPI_Fortran_HAS_MODULE}")
#        set_property(TARGET  mini-mdft
#                APPEND PROPERTY
#                COMPILE_DEFINITIONS "MPI_FORTRAN_HAS_MODULE")
#  endif()

#  include_directories(${MPI_Fortran_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH})
#  set_property(TARGET  mini-mdft APPEND_STRING PROPERTY COMPILE_FLAGS
#        " ${MPI_Fortran_COMPILE_FLAGS}  ${MPI_C_COMPILE_FLAGS} ")
#  set_property(TARGET  mini-mdft APPEND_STRING PROPERTY LINK_FLAGS
#        " ${MPI_Fortran_LINK_FLAGS} ${MPI_C_LINK_FLAGS} ")
#  target_link_libraries(mini-mdft ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES})

elseif(WIN32)
  message(STATUS "WIN32 platform not yet supported...")
endif()


#####################################
# Tell how to install this executable
#####################################

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${FOOEXE} RUNTIME DESTINATION bin)
