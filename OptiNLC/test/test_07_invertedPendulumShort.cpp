#include <fstream>
#include <iostream>
#include <vector>

#include "OptiNLC_Solver.h"

#include "fileCompare.h"
#include <catch2/catch.hpp>

constexpr int numberOfRuns = 100;

// Define the state and input sizes
constexpr int StateSize       = 4;
constexpr int InputSize       = 1;
constexpr int ConstraintsSize = 0;
constexpr int ControlPoints   = 40;

OptiNLC_Solver<double, 1, 4, 0, ControlPoints> solveExample(){

  OptiNLC_Options options;
  options.setDefaults();
  options.intermediateIntegration = 6;
  options.OSQP_max_iter           = 10000;
  options.maxNumberOfIteration    = 10000;
  options.OSQP_time_limit         = 1.;
  options.OptiNLC_ACC             = 1.e-10;
  options.OSQP_eps_abs            = 1.e-12;
  options.OSQP_eps_rel            = 1.e-12;
  options.timeStep                = 0.09;
  options.debugPrint = false;

  // Create an instance of the OCP class with the appropriate horizon and time step

  OptiNLC_OCP<double, InputSize, StateSize, ConstraintsSize, ControlPoints> ocp( &options );


  VECTOR<double, StateSize> initialState = { 0.0, 0.0, -1.0, 0.0 };
  VECTOR<double, InputSize> initialInput = { 0.0 };

  // Define a simple dynamic model as a lambda function
  ocp.setDynamicModel( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input,
                            VECTOR<double, StateSize>& derivative, double currentTime, void* userData ) {
    const double m2a     = 1.0;
    double       epsilon = m2a / ( m2a + 10. );
    derivative[0]        = state[1];
    derivative[1]        = -epsilon * state[2] + input[0];
    derivative[2]        = state[3];
    derivative[3]        = state[2] - input[0];
  } );

  // Define an objective function as a lambda function, states and input at tf
  ocp.setObjectiveFunction( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input, double current_time ) {
    // Define your objective function here
    // For example, you can minimize the squared sum of state variables
    // Compute squared norm manually
    double stateSquaredNorm = 0.0;
    for( int i = 0; i < StateSize; ++i )
    {
      stateSquaredNorm += state[i] * state[i];
    }

    double inputSquaredNorm = 0.0;
    for( int i = 0; i < InputSize; ++i )
    {
      inputSquaredNorm += input[i] * input[i];
    }

    return stateSquaredNorm + 0.001 * inputSquaredNorm;
  } );


  // Define a simple input update method
  ocp.setInputUpdate(
    [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input, double currentTime, void* userData ) {
      VECTOR<double, InputSize> updatedInput = { 0.0 };
      return updatedInput;
    } );
  // Define a states constraints method
  ocp.setUpdateStateLowerBounds( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input ) {
    VECTOR<double, StateSize> stateConstraints;
    stateConstraints.setConstant( -1e3 );
    return stateConstraints;
  } );
  // Define a states constraints method
  ocp.setUpdateStateUpperBounds( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input ) {
    VECTOR<double, StateSize> stateConstraints;
    stateConstraints.setConstant( 1e3 );
    return stateConstraints;
  } );
  // Define a inputs constraints method
  ocp.setUpdateInputLowerBounds( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input ) {
    VECTOR<double, InputSize> inputConstraints;
    inputConstraints.setConstant( -1e3 );
    return inputConstraints;
  } );
  // Define a inputs constraints method
  ocp.setUpdateInputUpperBounds( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input ) {
    VECTOR<double, InputSize> inputConstraints;
    inputConstraints.setConstant( 1e3 );
    return inputConstraints;
  } );

  // Define a functions constraints method
  //
  ocp.setUpdateFunctionConstraints( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input ) {
    VECTOR<double, ConstraintsSize> functionsConstraints;
    functionsConstraints.setConstant( 0.0 );
    return functionsConstraints;
  } );
  ocp.setUpdateFunctionConstraintsLowerBounds( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input ) {
    VECTOR<double, ConstraintsSize> functionsConstraints;
    functionsConstraints.setConstant( -std::numeric_limits<double>::infinity() );
    return functionsConstraints;
  } );
  ocp.setUpdateFunctionConstraintsUpperBounds( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input ) {
    VECTOR<double, ConstraintsSize> functionsConstraints;
    functionsConstraints.setConstant( std::numeric_limits<double>::infinity() );
    return functionsConstraints;
  } );
  // user defined function for all constraints


  OptiNLC_Solver<double, InputSize, StateSize, ConstraintsSize,ControlPoints> _solver(ocp);
  _solver.solve(0.0, initialState, initialInput);
  return _solver;
}

void writeData(std::string filename, std::array<std::array<float, StateSize>, ControlPoints + 1> data) { 
  std::ofstream dataFile(filename);
  if (dataFile.is_open()) {
    for (int i = 0; i < ControlPoints + 1; ++i) {
        dataFile << i << " " << data[i][0] << " " << data[i][1]
                << " " << data[i][2] << " " << data[i][3] << std::endl;
    }
    dataFile.close();
    std::cout << "Data saved to " << filename << std::endl;
  } else {
      std::cerr << "Unable to open file for writing." << std::endl;
      REQUIRE( false );
  }
}

void writeRawData(std::string filename, OptiNLC_Solver<double, 1, 4, 0, ControlPoints> solver){
std::ofstream dataFile(filename);
auto opt_x = solver.get_optimal_states();
auto time = solver.getTime();
  // Check if the file is open
  if (dataFile.is_open()) {
      for (int i = 0; i < time.size(); ++i) {
          dataFile << time[i]<< " " << opt_x[StateSize*i + 0] << " " << opt_x[StateSize*i + 1]
                  << " " << opt_x[StateSize*i + 2]
                  << " " << opt_x[StateSize*i + 3] << std::endl;
      }
      dataFile.close();
      std::cout << "Data saved to eigen_data_01.txt" << std::endl;
  } else {
      std::cerr << "Unable to open file for writing." << std::endl;
      REQUIRE( false );
  }
}


TEST_CASE( "TEST Inverted Pendulum short:" )
{
  std::array<std::array<float, StateSize>, ControlPoints + 1> min;
  std::array<std::array<float, StateSize>, ControlPoints + 1> max;
  std::array<std::array<float, StateSize>, ControlPoints + 1> maxDiff;
  for (auto& row : maxDiff)
    row.fill(FLT_MAX);
  for (auto& row : min)
    row.fill(FLT_MAX);
  for (auto& row : max)
    row.fill(-FLT_MAX);

  std::vector<OptiNLC_Solver<double, 1, 4, 0, ControlPoints>> solvers;
  for (int i = 0; i < numberOfRuns; i++){
    solvers.push_back(solveExample());
    {
      //writeRawData("eigen_data_0_" + std::to_string(i) + ".txt", solvers[i]);
    }

  }
  for (int i = 0; i < numberOfRuns; i++){
    auto _solver = solvers[i];
  
    auto opt_x =  _solver.get_optimal_states();
    auto time = _solver.getTime();
    REQUIRE(time.size() == ControlPoints + 1);
    for (int i = 0; i < time.size(); ++i) {
      auto a = opt_x[StateSize*i + 0];
      auto b = opt_x[StateSize*i + 1];
      auto c = opt_x[StateSize*i + 2];
      auto d = opt_x[StateSize*i + 3];
      if (a < min[i][0])
        min[i][0] = a;
      if (b < min[i][1])
        min[i][1] = b;
      if (c < min[i][2])
        min[i][2] = c;
      if (d < min[i][3])
        min[i][3] = d;
      if (a > max[i][0])
        max[i][0] = a;
      if (b > max[i][1])
        max[i][1] = b;
      if (c > max[i][2])
        max[i][2] = c;
      if (d > max[i][3])
        max[i][3] = d;
    }
  }

  for (int i = 0; i < ControlPoints + 1; i++){
    for (int j = 0; j < StateSize; j++) {
      maxDiff[i][j] = max[i][j] - min[i][j];
    }
  }

  writeData("eigen_data_max.txt", max);
  writeData("eigen_data_min.txt", min);
  writeData("eigen_data_07_diff.txt", maxDiff);

  //debug
  //find max different solutions
  int ii = 0;
  int jj = 0;
  float maxDiffOneValue = 0;
  for (int i = 0; i < ControlPoints + 1; i++){
    for (int j = 0; j < StateSize; j++) {
      if (maxDiff[i][j] > maxDiffOneValue){
        ii = i;
        jj = j;
      }
    }
  }
  bool ma = false;
  bool mi = false;
  for (int i = 0; i < numberOfRuns; i++){
    auto opt_x =  solvers[i].get_optimal_states();
    if (opt_x[StateSize*ii + jj] == max[ii][jj]){
      writeRawData("eigen_data_MAX.txt", solvers[i]);
      ma = true;
    }
    if (opt_x[StateSize*ii + jj] == min[ii][jj]){
      writeRawData("eigen_data_MIN.txt", solvers[i]);
      mi = true;
    }
    if (ma && mi){
      break;
    }
  }
}
