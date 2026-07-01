# (using CMAKE_CURRENT_SOURCE_DIR is necessary because sources are in a different directory than the target library (trans_${prec}))
ecbuild_list_add_pattern(
    LIST ectrans4py_src
    GLOB ${CMAKE_CURRENT_SOURCE_DIR}/*.F90
    QUIET
)

# Build one ectrans4py library per enabled precision: ectrans4py_sp (linking
# trans_sp) and/or ectrans4py_dp (linking trans_dp). The Python layer loads
# whichever is installed and matches its field-data dtype accordingly (see the
# _REAL / ECTRANS4PY_PRECISION handling in __init__.py). The serial LAM wrappers
# use etrans, so etrans_${prec} is linked when ETRANS is enabled.
foreach( prec sp dp )
    if( prec STREQUAL "sp" AND NOT HAVE_SINGLE_PRECISION )
        continue()
    endif()
    if( prec STREQUAL "dp" AND NOT HAVE_DOUBLE_PRECISION )
        continue()
    endif()

    set( ectrans4py_libs fiat trans_${prec} parkind_${prec} )
    if( HAVE_ETRANS )
        list( APPEND ectrans4py_libs etrans_${prec} )
    endif()

    set( ectrans4py_src_${prec} ${ectrans4py_src} )
    if( prec STREQUAL "sp" )
        # sp2gp_fft1d4py uses the double-precision FFTW module (tpm_fftw_dp) and
        # REAL64 data; it is a niche (self-described "not tested") LAM helper not
        # needed for the single-precision interface, so it is left out of _sp.
        list( FILTER ectrans4py_src_${prec} EXCLUDE REGEX "sp2gp_fft1d4py" )
    endif()

    ecbuild_add_library(
          TARGET           ectrans4py_${prec}
          LINKER_LANGUAGE  Fortran
          SOURCES          ${ectrans4py_src_${prec}}
          PUBLIC_LIBS      ${ectrans4py_libs}
    )
    # Per-precision Fortran module directory so the two targets do not clash on
    # setup_core4py_mod.mod (they compile the same sources at different precision).
    set_target_properties( ectrans4py_${prec} PROPERTIES
          Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module_${prec} )
endforeach()
