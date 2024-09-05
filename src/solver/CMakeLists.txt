#
file(GLOB srcs ./*.cpp)

if(ENABLED_COUNT EQUAL 0)
    set(MY_LIB_DIR "${CMAKE_SOURCE_DIR}/3rd_lib/SuiteSparse/lib")
    set(MY_LIBS
        "${MY_LIB_DIR}/libklu.a"
        "${MY_LIB_DIR}/libcholmod.a"
        "${MY_LIB_DIR}/libccolamd.a"
        "${MY_LIB_DIR}/libcolamd.a"
        "${MY_LIB_DIR}/libcamd.a"
        "${MY_LIB_DIR}/libbtf.a"
        "${MY_LIB_DIR}/libamd.a"
        "${MY_LIB_DIR}/libsuitesparseconfig.a"
    )
    # add_library(solver STATIC ${srcs} ${MY_LIBS} m)
    add_library(solver STATIC ${srcs})
    target_link_libraries(solver PRIVATE ${MY_LIBS} m)

elseif(USE_EIGEN)
    include_directories(${CMAKE_SOURCE_DIR}/3rd_lib/eigen-3.4.0)
    add_library(solver STATIC ${srcs})

elseif(USE_KLU)
    set(MY_LIB_DIR "${CMAKE_SOURCE_DIR}/3rd_lib/SuiteSparse/lib")
    set(MY_LIBS
        "${MY_LIB_DIR}/libklu.a"
        "${MY_LIB_DIR}/libcholmod.a"
        "${MY_LIB_DIR}/libccolamd.a"
        "${MY_LIB_DIR}/libcolamd.a"
        "${MY_LIB_DIR}/libcamd.a"
        "${MY_LIB_DIR}/libbtf.a"
        "${MY_LIB_DIR}/libamd.a"
        "${MY_LIB_DIR}/libsuitesparseconfig.a"
    )
    # add_library(solver STATIC ${srcs} ${MY_LIBS} m)
    add_library(solver STATIC ${srcs})
    target_link_libraries(solver PRIVATE ${MY_LIBS})

elseif(USE_SUPERLU)
    # include_directories(${CMAKE_SOURCE_DIR}/3rd_lib/superlu_dist-master/include)
    include_directories("/root/local/include")
    set(SUPERLU_LIB "${CMAKE_SOURCE_DIR}/3rd_lib/superlu_dist/lib")
    set(PARMETIS_LIB "/root/local/lib")
    set(MY_LIBS1
        "${SUPERLU_LIB}/libsuperlu_dist.a"
        "${PARMETIS_LIB}/libparmetis.a"
        "${PARMETIS_LIB}/libmetis.a"
        "${PARMETIS_LIB}/libGKlib.a"
    )

    if(DEFINED ENV{MKLROOT})
        set(MKL_INCLUDE_DIR $ENV{MKLROOT}/include)
        # set(MKL_LIBRARY_DIR $ENV{MKLROOT}/lib/intel64)
        set(MKL_LIBRARY_DIR $ENV{MKLROOT}/lib)

        include_directories(${MKL_INCLUDE_DIR})

        set(MY_LIBS2
            "${MKL_LIBRARY_DIR}/libmkl_intel_lp64.a"
            "${MKL_LIBRARY_DIR}/libmkl_intel_thread.a"
            "${MKL_LIBRARY_DIR}/libmkl_core.a"
            "/opt/intel/oneapi/compiler/2024.2/lib/libiomp5.a"
        )
    else()
        message(FATAL_ERROR "MKLROOT environment variable is not set.")
    endif()

    # add_library(solver STATIC ${srcs} ${MY_LIBS1} ${MY_LIBS2})
    add_library(solver STATIC ${srcs})
    target_link_libraries(solver PRIVATE ${MY_LIBS1} ${MY_LIBS2})

elseif(USE_MUMPS)
    
elseif(USE_MINISPICE)
    
endif()
