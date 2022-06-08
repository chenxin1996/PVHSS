cmake_minimum_required(VERSION 3.22)
project(DecisionTree_PVHSS)

set(CMAKE_CXX_STANDARD 14)

add_executable(DecisionTree_PVHSS main.cpp Eval_DT.h PKE.h PVHSS.h tool.h)

target_link_libraries(DecisionTree_PVHSS libgmp.so libflint.so libpbc.so m.so librelic.so /usr/local/lib/libntl.so)