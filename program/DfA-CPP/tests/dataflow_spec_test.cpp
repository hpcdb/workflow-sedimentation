/* 
 * File:   dataflow_specification.cpp
 * Author: vitor
 *
 * Created on March 31, 2018, 1:15 PM
 */

#include "dataflow.h"

#include <stdlib.h>
#include <iostream>

/*
 * Simple C++ Test Suite
 */

void test_dataflow(){
    Dataflow dataflow = Dataflow("clothing");
    dataflow.save();
}

void test1() {
    std::cout << "newsimpletest test 1" << std::endl;
}

void test2() {
    std::cout << "newsimpletest test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (newsimpletest) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% newsimpletest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;
    
    std::cout << "%TEST_STARTED% test_dataflow" << std::endl;
    test_dataflow();
    std::cout << "%TEST_FINISHED% time=0 test_dataflow" << std::endl;

//    std::cout << "%TEST_STARTED% test1 (newsimpletest)" << std::endl;
//    test1();
//    std::cout << "%TEST_FINISHED% time=0 test1 (newsimpletest)" << std::endl;

//    std::cout << "%TEST_STARTED% test2 (newsimpletest)\n" << std::endl;
//    test2();
//    std::cout << "%TEST_FINISHED% time=0 test2 (newsimpletest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}
