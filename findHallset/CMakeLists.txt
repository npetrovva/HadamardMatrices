cmake_minimum_required(VERSION 3.10)
project(findHallset)

add_executable(setfindmin setfindmin.cpp)
add_executable(Hallset Hallset.cpp)


add_custom_target(run_all
    COMMAND ./setfindmin "ls Had*" "ls minHad*"
    COMMAND ./Hallset "ls minHad*" "ls swminHad*"
    COMMAND ./setfindmin "ls swminHad*" "ls resswminHad*"
)
