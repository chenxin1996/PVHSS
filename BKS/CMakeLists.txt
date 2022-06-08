cmake_minimum_required(VERSION 3.21)
project(BKS)

set(CMAKE_CXX_STANDARD 14)

add_executable(BKS main.cpp function.h HSS.h GenData.h tool.h PKE.h)

target_link_libraries(BKS libgmp.so libflint.so libpbc.so m.so librelic_s.a /usr/local/lib/libntl.so)