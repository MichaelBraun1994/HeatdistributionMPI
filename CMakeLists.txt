cmake_minimum_required(VERSION 3.17)
project (heatdistribution)

find_package(MPI)

include_directories(include)
file(GLOB SOURCES "src/*.c")

add_executable(main ${SOURCES})
target_link_libraries(main PRIVATE ${MPI_C_LIBRARIES} m)
