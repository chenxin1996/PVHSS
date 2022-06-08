cmake_minimum_required(VERSION 3.22)
project(DecisionTree_BGV_CKKS)

set(CMAKE_CXX_STANDARD 17)

add_executable(DecisionTree_BGV_CKKS main.cpp Eval_DT.h)

target_link_libraries(DecisionTree_BGV_CKKS libgmp.so libhelib.so /usr/local/lib/libntl.so libpthread.a libflint.so)





