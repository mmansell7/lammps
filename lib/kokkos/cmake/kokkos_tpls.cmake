KOKKOS_CFG_DEPENDS(TPLS OPTIONS)
KOKKOS_CFG_DEPENDS(TPLS DEVICES)

FUNCTION(KOKKOS_TPL_OPTION PKG DEFAULT)
  KOKKOS_ENABLE_OPTION(${PKG} ${DEFAULT} "Whether to enable the ${PKG} library")
  KOKKOS_OPTION(${PKG}_DIR "" PATH "Location of ${PKG} library")
  SET(KOKKOS_ENABLE_${PKG} ${KOKKOS_ENABLE_${PKG}} PARENT_SCOPE)
  SET(KOKKOS_${PKG}_DIR  ${KOKKOS_${PKG}_DIR} PARENT_SCOPE)
ENDFUNCTION()

KOKKOS_TPL_OPTION(HWLOC   Off)
KOKKOS_TPL_OPTION(LIBNUMA Off)
KOKKOS_TPL_OPTION(MEMKIND Off)
KOKKOS_TPL_OPTION(CUDA    Off)
KOKKOS_TPL_OPTION(LIBRT   Off)
KOKKOS_TPL_OPTION(LIBDL   On)

IF(Trilinos_ENABLE_Kokkos AND TPL_ENABLE_HPX)
SET(HPX_DEFAULT ON)
ELSE()
SET(HPX_DEFAULT OFF)
ENDIF()
KOKKOS_TPL_OPTION(HPX ${HPX_DEFAULT})

IF(Trilinos_ENABLE_Kokkos AND TPL_ENABLE_PTHREAD)
SET(PTHREAD_DEFAULT ON)
ELSE()
SET(PTHREAD_DEFAULT OFF)
ENDIF()
KOKKOS_TPL_OPTION(PTHREAD ${PTHREAD_DEFAULT})


#Make sure we use our local FindKokkosCuda.cmake
KOKKOS_IMPORT_TPL(HPX INTERFACE)
KOKKOS_IMPORT_TPL(CUDA INTERFACE)
KOKKOS_IMPORT_TPL(HWLOC)
KOKKOS_IMPORT_TPL(LIBNUMA)
KOKKOS_IMPORT_TPL(LIBRT)
KOKKOS_IMPORT_TPL(LIBDL)
KOKKOS_IMPORT_TPL(MEMKIND)
KOKKOS_IMPORT_TPL(PTHREAD INTERFACE)

#Convert list to newlines (which CMake doesn't always like in cache variables)
STRING(REPLACE ";" "\n" KOKKOS_TPL_EXPORT_TEMP "${KOKKOS_TPL_EXPORTS}")
#Convert to a regular variable
UNSET(KOKKOS_TPL_EXPORTS CACHE)
SET(KOKKOS_TPL_EXPORTS ${KOKKOS_TPL_EXPORT_TEMP})
