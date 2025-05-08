#include <catch2/catch.hpp>
#include "OptiNLC_Data.h"

TEST_CASE("TEST Datastructures")
{
    int SIZE = 10;
    VECTOR<double, 10> input = { 1.0, 0.0};
    
    SECTION("Size")
    {
        REQUIRE(input.data()[0] == 1.0);
        REQUIRE(input.data()[1] == 0.0);
        REQUIRE(input.size() == SIZE);
    }
    SECTION("Setter")
    {
        REQUIRE(input.data()[0] == 1.0);
        REQUIRE(input.data()[1] == 0.0);
        input.setConstant(30.0);
        REQUIRE(input.data()[0] == 30.0);
        REQUIRE(input.data()[1] == 30.0);
        REQUIRE(input.data()[SIZE-1] == 30.0);
        input.set(10.0, 1);
        REQUIRE(input.data()[0] == 30.0);
        REQUIRE(input.data()[1] == 10.0);
        REQUIRE(input.data()[SIZE-1] == 30.0);
    }
    SECTION("Vector Setter")
    {
        VECTOR<double, 10> input2 = { 2.0 };
        input.set(input2);
        REQUIRE(input.data()[0] == 2.0);
        REQUIRE(input.data()[1] == 0.0);
        VECTOR<double, 10> input3 = { 2.0 , 2.0, 2.0 };
        VECTOR<double, 10> input4 = { 4.0 };
        input3.set(input4);
        REQUIRE(input3.data()[0] == 4.0);
        REQUIRE(input3.data()[1] > -0.1);
        REQUIRE(input3.data()[1] < 0.1);
    }
    SECTION("Eigen Matrix")
    {
        auto eigenMatrix = input.getDataAsEigen();
        REQUIRE(eigenMatrix.size() == 10);
    }
    SECTION("Printer")
    {
        input.print();
        input.print_T();
        //input.print(0,0);
        REQUIRE(true);
    }
    SECTION("Operators")
    {
        //how to test? +-*+=-=
        REQUIRE(true);
    }
    SECTION("Methods")
    {
        VECTOR<double, 10> input2 = { 1.0, 0.0, 10.6};
        REQUIRE(input2.sum() == 11.6);
        REQUIRE(input2.norm_inf() == 10.6);
    }
    SECTION("minMAX")
    {
        VECTOR<double, 10> input2 = { 1.0, 0.0, 10.6};
        input2.cwiseMin(5.0);
        REQUIRE(input2[0] == 1.0);
        REQUIRE(input2[2] == 5.0);
        input2.cwiseMax(3.0);
        REQUIRE(input2[0] == 3.0);
        REQUIRE(input2[2] == 5.0);
    }
    SECTION("minMAXVector")
    {
        VECTOR<double, 10> input2 = { 0.0, 10.6};
        input.cwiseMin(input2);
        REQUIRE(input[0] == 0.0);
        REQUIRE(input[1] == 0.0);
        input2.cwiseMax(input);
        REQUIRE(input2[0] == 0.0);
        REQUIRE(input2[1] == 10.6);
    }
    SECTION("norm"){
        //positive
        VECTOR<double, 2> input2 = { 3.0, 4.0};
        REQUIRE(input2.norm(1) == 7.0);
        REQUIRE(input2.norm() == 5.0);
        REQUIRE(input2.norm(2) == 5.0);
        VECTOR<double, 3> input3 = { 3.0, 4.0, 5.0};
        REQUIRE(input3.norm(3) > 5.9);
        REQUIRE(input3.norm(3) < 6.1);

        //zero
        VECTOR<double, 2> input4 = { 0, 0};
        REQUIRE(input4.norm(1) == 0);
        REQUIRE(input4.norm() == 0);
        REQUIRE(input4.norm(2) == 0);
        VECTOR<double, 10> input5 = { 0, 0, 0};
        REQUIRE(input5.norm(3) == 0);

        //negative
        VECTOR<double, 2> input6 = { -3.0, -4.0};
        REQUIRE(input6.norm(1) == 7.0);
        REQUIRE(input6.norm() == 5.0);
        REQUIRE(input6.norm(2) == 5.0);
        VECTOR<double, 3> input7 = { -3, -4, 5};
        REQUIRE(input7.norm(3) > 5.9);
        REQUIRE(input7.norm(3) < 6.1);
    }
}
