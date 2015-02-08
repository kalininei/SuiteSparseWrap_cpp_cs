cmake_minimum_required(VERSION 2.8)

IF(WIN32)
	#Force only release vs build since precompiled suitesparse libs have release compilation
	#and will not allow to compile in a debug mode
	SET(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "" FORCE)
ENDIF(WIN32)

PROJECT(CHOLMOD_WRAPPER)

# ------------- Detect SuiteSparse libraries:
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

IF(WIN32)
	#use archived libraries for visual studio compillation
	OPTION(SUITESPARSE_LOCAL_DISTR "Use custom path to suitesparse and metis" ON)
	#check if we need 64 bit default files
	if (CMAKE_CL_64)
		SET(WINLIB_DIR lib64)
	else (CMAKE_CL_64)
		SET(WINLIB_DIR lib)
	endif(CMAKE_CL_64)
	SET(SUITESPARSE_LOCAL_INCDIR ${CMAKE_SOURCE_DIR}/winlib/include CACHE
		PATH "Search for suitesparse and metis here if SUITESPARSE_LOCAL_DISTR is ON")
	SET(SUITESPARSE_LOCAL_LIBDIR ${CMAKE_SOURCE_DIR}/winlib/${WINLIB_DIR} CACHE
		PATH "Search for suitesparse headers here if SUITESPARSE_LOCAL_DISTR is ON")
ELSE(WIN32)
	#unix optionally uses libraries which have been already installed in the system
	OPTION(SUITESPARSE_LOCAL_DISTR "Use custom path to suitesparse and metis" OFF)
	SET(SUITESPARSE_LOCAL_INCDIR ~/bin/include CACHE
		PATH "Search for suitesparse headers here if SUITESPARSE_LOCAL_DISTR is ON")
	SET(SUITESPARSE_LOCAL_LIBDIR ~/bin/lib CACHE
		PATH "Search for suitesparse and metis here if SUITESPARSE_LOCAL_DISTR is ON")
ENDIF(WIN32)


find_package(SuiteSparse REQUIRED)

# ------------- define output paths
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
SET(ARCHIVE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# ------------- cholmod_wrapper
set(SOURCES cholmod_wrapper.cpp)
if (WIN32)
	set(SOURCES ${SOURCES} dll.def)
endif()
set (HEADERS cholmod_wrapper.h)

ADD_LIBRARY(cholmod_wrapper SHARED ${SOURCES} ${HEADERS})
INCLUDE_DIRECTORIES(cholmod_wrapper ${SUITESPARSE_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(cholmod_wrapper ${SUITESPARSE_LIBRARIES})

# -------------- cholmod installation
SET_TARGET_PROPERTIES(cholmod_wrapper PROPERTIES PUBLIC_HEADER "${HEADERS}")
INSTALL(TARGETS cholmod_wrapper
		RUNTIME DESTINATION .
		LIBRARY DESTINATION .
		ARCHIVE DESTINATION .
		PUBLIC_HEADER DESTINATION include
)

# --------------- install blas, lapack, quadmath, fortran-3, gcc_s_dw2-1
if (WIN32)
	OPTION(COPY_PRECOMPILED_BLAS_LAPACK_LIBS
		"Copy default blas/lapack binaries to installation directory" ON)
	if (CMAKE_CL_64)
		SET(LIBGCC_S_LIBNAME libgcc_s_sjlj-1.dll)
	else(CMAKE_CL_64)
		SET(LIBGCC_S_LIBNAME libgcc_s_dw2-1.dll)
	endif(CMAKE_CL_64)
	SET(LAPACK_FRAMEWORK
		${SUITESPARSE_LOCAL_LIBDIR}/libblas.dll;
		${SUITESPARSE_LOCAL_LIBDIR}/${LIBGCC_S_LIBNAME};
		${SUITESPARSE_LOCAL_LIBDIR}/libgfortran-3.dll;
		${SUITESPARSE_LOCAL_LIBDIR}/liblapack.dll;
		${SUITESPARSE_LOCAL_LIBDIR}/libquadmath-0.dll;
	)

	if(COPY_PRECOMPILED_BLAS_LAPACK_LIBS)
		INSTALL(FILES ${LAPACK_FRAMEWORK}
			DESTINATION .)
	endif(COPY_PRECOMPILED_BLAS_LAPACK_LIBS)

endif(WIN32)


# -------------- test program
add_subdirectory(cpptest)
