cmake_minimum_required(VERSION 3.21)
project(CZ)

set(CMAKE_CXX_STANDARD 14)

add_executable(CZ main.cpp function.h VHSS.h GenData.h PKE.h tool.h)

target_link_libraries(CZ libgmp.so libflint.so libpbc.so m.so librelic_s.a /usr/local/lib/libntl.so)