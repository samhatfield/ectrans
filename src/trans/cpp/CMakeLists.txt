# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define C++ common library, which contains all the C++ sources and requires no precision
# preprocessing
ecbuild_list_add_pattern( LIST ectrans_cpp_common_src GLOB external_cpp/*.cpp internal/*.cpp QUIET )
ecbuild_add_library(
  TARGET          ectrans_cpp_common
  SOURCES         ${ectrans_cpp_common_src}
  LINKER_LANGUAGE CXX
  PUBLIC_LIBS     fiat
)
target_include_directories(
  ectrans_cpp_common
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/internal
)

function( generate_backend_sources )
  set( options )
  set( oneValueArgs BACKEND DESTINATION OUTPUT )
  set( multiValueArgs )

  cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  set( backend ${_PAR_BACKEND} )
  set( destination ${_PAR_DESTINATION} )
  file( MAKE_DIRECTORY ${destination}/external )

  ecbuild_list_add_pattern( LIST files GLOB external/*.F90 QUIET )

  set( outfiles )
  foreach( file_i ${files} )
    get_filename_component( outfile_name    ${file_i} NAME )
    get_filename_component( outfile_name_we ${file_i} NAME_WE )
    get_filename_component( outfile_ext     ${file_i} EXT )
    get_filename_component( outfile_dir     ${file_i} DIRECTORY )
    set( outfile "${destination}/${file_i}" )
    ecbuild_debug( "Generate ${outfile}" )
    generate_file(
      BACKEND ${backend}
      INPUT   ${CMAKE_CURRENT_SOURCE_DIR}/${file_i}
      OUTPUT  ${outfile}
    )
    list( APPEND outfiles ${outfile} )
  endforeach( file_i )
  set( ${_PAR_OUTPUT} ${outfiles} PARENT_SCOPE )
endfunction( generate_backend_sources )

set( BUILD_INTERFACE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include/ectrans )

foreach( prec dp sp )
  if( HAVE_${prec} )
    set( GENERATED_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated/ectrans_cpp_${prec} )

    generate_backend_includes(
      BACKEND     cpp_${prec}
      TARGET      ectrans_cpp_${prec}_includes
      DESTINATION ${BUILD_INTERFACE_INCLUDE_DIR}
      INCLUDE_DIRECTORY ${PROJECT_SOURCE_DIR}/src/trans/include
    )
    generate_backend_sources(
      BACKEND     cpp_${prec}
      OUTPUT      ectrans_cpp_${prec}_src
      DESTINATION ${GENERATED_SOURCE_DIR}
    )

    # Define a library containing the Fortran interface files which call the C++ common library
    ecbuild_add_library(
      TARGET          ectrans_cpp_${prec}
      SOURCES         ${ectrans_cpp_${prec}_src}
      LINKER_LANGUAGE Fortran
      PUBLIC_INCLUDES $<INSTALL_INTERFACE:include/ectrans>
      PUBLIC_LIBS     ectrans_cpp_common ectrans_cpp_${prec}_includes
    )

    if( prec STREQUAL sp )
      target_compile_definitions( ectrans_cpp_${prec} PRIVATE TRANS_SINGLE )
    endif()

    # This interface library is for backward compatibility, and provides the older includes
    ecbuild_add_library( TARGET trans_cpp_${prec} TYPE INTERFACE )
    target_include_directories(
      trans_cpp_${prec}
      INTERFACE $<BUILD_INTERFACE:${BUILD_INTERFACE_INCLUDE_DIR}/trans_cpp_${prec}>
    )
    target_include_directories(
      trans_cpp_${prec}
      INTERFACE $<INSTALL_INTERFACE:include/ectrans/trans_cpp_${prec}>
    )
    target_link_libraries( trans_cpp_${prec} INTERFACE fiat ectrans_cpp_${prec} parkind_${prec})
  endif()
endforeach()

## Install trans interface
install( DIRECTORY ${BUILD_INTERFACE_INCLUDE_DIR}/ DESTINATION include/ectrans )
