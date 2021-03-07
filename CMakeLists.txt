cmake_minimum_required(VERSION 3.17)
project(MD2 C)

set(CMAKE_C_STANDARD 11)

include_directories(src)
find_package(OpenMP)

add_executable(MD2
        src/main.c)



find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})


target_link_libraries(MD2 ${MPI_LIBRARIES})

target_link_libraries(MD2 m)

set(CMAKE_C_FLAGS "-O3 -Wall -Wextra -fopenmp")
set(CMAKE_C_FLAGS_DEBUG "-g -Wall -Wextra -fopenmp --coverage")

if("${CMAKE_C_COMPILER_ID}" MATCHES "(Apple)?[Cc]lang"
        OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "(Apple)?[Cc]lang")
    message("Building debug with llvm Code Coverage")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fprofile-instr-generate -fcoverage-mapping")
elseif(CMAKE_COMPILER_IS_GNUC)
    message("Building debug with lcov Code Coverage Tools")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} --coverage")
endif()

if(MPI_COMPILE_FLAGS)
    set_target_properties(MD2 PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(MD2 PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")

endif()



