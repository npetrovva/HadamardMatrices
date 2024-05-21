cmake_minimum_required(VERSION 3.10)
project(findQuadruples)

add_executable(setfindmin setfindmin.cpp)
add_executable(quadruples quadruples.cpp)
add_executable(quadruples_transpose quadruples_transpose.cpp)

add_custom_target(run_all
    COMMAND ./setfindmin ''ls *.txt" "ls min*.txt"
    COMMAND ./quadruples "ls min*.txt" "ls swmin*.txt"
    COMMAND ./setfindmin "ls swmin*.txt" "ls res.txt"
    COMMAND ./setfindmin "ls *.txt" "ls min*.txt"
    COMMAND ./quadruples_transpose "ls min*.txt" "ls swtransp_min*.txt"
    COMMAND ./setfindmin "ls swtransp_min*.txt" "ls transpres*.txt"
)
