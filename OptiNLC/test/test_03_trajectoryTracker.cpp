#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "OptiNLC_Data.h"
#include "OptiNLC_Options.h"
#include "OptiNLC_Solver.h"
#include <catch2/catch.hpp>

class ReadData
{
public:

  ReadData() {}

  struct piecewisePolynomialVector
  {
    std::vector<double> breaks;
    std::vector<double> coef1;
    std::vector<double> coef2;
    std::vector<double> coef3;
    std::vector<double> coef4;
  } pp_x, pp_y, pp_psi;

  piecewisePolynomialVector
  get_pp_x()
  {
    return pp_x;
  }

  piecewisePolynomialVector
  get_pp_y()
  {
    return pp_y;
  }

  piecewisePolynomialVector
  get_pp_psi()
  {
    return pp_psi;
  }

  bool
  readTextData()
  {
    std::string filename_x   = "/tmp/OptiNLC/course_x.txt";
    std::string filename_y   = "/tmp/OptiNLC/course_y.txt";
    std::string filename_psi = "/tmp/OptiNLC/course_phi.txt";
    // Open the file
    std::ifstream file_x( filename_x );
    std::ifstream file_y( filename_y );
    std::ifstream file_psi( filename_psi );
    // Check if the file is open
    if( !file_x.is_open() || !file_y.is_open() | !file_psi.is_open() )
    {
      std::cerr << "Error opening file: " << std::endl;
      return false; // Return an error code
    }
    // Read data from the file and push it into the vector

    double value;
    while( file_x >> value )
    {
      pp_x.breaks.push_back( value );
      // Read the remaining values for the other vectors
      file_x >> value;
      pp_x.coef1.push_back( value );
      file_x >> value;
      pp_x.coef2.push_back( value );
      file_x >> value;
      pp_x.coef3.push_back( value );
      file_x >> value;
      pp_x.coef4.push_back( value );
    }
    file_x.close();
    value = 0;
    while( file_y >> value )
    {
      pp_y.breaks.push_back( value );
      // Read the remaining values for the other vectors
      file_y >> value;
      pp_y.coef1.push_back( value );
      file_y >> value;
      pp_y.coef2.push_back( value );
      file_y >> value;
      pp_y.coef3.push_back( value );
      file_y >> value;
      pp_y.coef4.push_back( value );
    }
    file_y.close();
    value = 0;
    while( file_psi >> value )
    {
      pp_psi.breaks.push_back( value );
      // Read the remaining values for the other vectors
      file_psi >> value;
      pp_psi.coef1.push_back( value );
      file_psi >> value;
      pp_psi.coef2.push_back( value );
      file_psi >> value;
      pp_psi.coef3.push_back( value );
      file_psi >> value;
      pp_psi.coef4.push_back( value );
    }
    file_psi.close();
    return true;
  }

  void
  test()
  {
    // std::cout<<pp_x.breaks.size()<<"\n";
  }
};

TEST_CASE( "OCP Test NEW FORMAT Tracker:" )
{
  std::cout << "OCP Test NEW FORMAT Tracker:" << std::endl;


  OptiNLC_Options options;
  options.setDefaults();

  options.intermediateIntegration = 2;
  options.OptiNLC_ACC             = 1.e-7;
  options.maxNumberOfIteration    = 500;
  options.OSQP_verbose            = false;
  options.OSQP_max_iter           = 400;
  options.OptiNLC_time_limit      = 0.1;

  // options.OSQP_time_limit            = 0.1;
  enum STATES
  {
    X,
    Y,
    PSI,
    V,
    S,
    DELTA,
    dDELTA,
    L
  };

  enum CONTROLS
  {
    ddDELTA,
    F
  };

  // Define the state and input sizes
  constexpr int StateSize       = 8;
  constexpr int InputSize       = 2;
  constexpr int ConstraintsSize = 0;
  constexpr int ControlPoints   = 10;
  double        simTime         = 0.50;
  options.timeStep              = ( simTime / ControlPoints );

  // Create an instance of the OCP class with the appropriate horizon and time step
  ReadData test;
  if( !test.readTextData() )
  {
    return;
  }
  double                    x_init       = test.pp_x.coef4[14];
  double                    y_init       = test.pp_y.coef4[14];
  double                    psi_init     = test.pp_psi.coef4[14];
  double                    s_init       = test.pp_x.breaks[14];
  VECTOR<double, StateSize> initialState = { x_init, y_init, psi_init, 0.1, s_init, 0.0, 0.0, 0.0 };
  VECTOR<double, InputSize> initialInput = { 0.00, 0.0 };
  OptiNLC_OCP<double, InputSize, StateSize, ConstraintsSize, ControlPoints> ocp( &options );
  // Define a simple dynamic model as a lambda function


  ocp.setDynamicModel( [&test]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input,
                                VECTOR<double, StateSize>& derivative, double currentTime, void* userData ) {
    const double l         = 4.2;
    const double tau       = 1.0;
    const double desired_v = 3.60;
    double       k         = desired_v / tau;
    derivative[X]          = state[V] * cos( state[PSI] );
    derivative[Y]          = state[V] * sin( state[PSI] );
    derivative[PSI]        = state[V] * tan( state[DELTA] ) / l;
    derivative[V]          = ( -( 1.0 / tau ) * state[V] ) + ( k * input[F] );
    derivative[S]          = state[V];
    derivative[DELTA]      = state[dDELTA];
    derivative[dDELTA]     = input[ddDELTA];
    int    index           = std::floor( state[S] );
    double ds              = state[S] - test.pp_x.breaks[index];
    double reference_x     = test.pp_x.coef1[index] * ds * ds * ds + test.pp_x.coef2[index] * ds * ds + test.pp_x.coef3[index] * ds
                       + test.pp_x.coef4[index];
    double reference_y = test.pp_y.coef1[index] * ds * ds * ds + test.pp_y.coef2[index] * ds * ds + test.pp_y.coef3[index] * ds
                       + test.pp_y.coef4[index];
    derivative[L] = ( state[V] - 3.250 ) * ( state[V] - 3.250 ) * 0.1 + ( state[Y] - reference_y ) * ( state[Y] - reference_y )
                  + ( state[X] - reference_x ) * ( state[X] - reference_x );
  } );

  // Define an objective function as a lambda function, states and input at tf
  ocp.setObjectiveFunction( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input, double current_time ) {
    // Define your objective function here
    // For example, you can minimize the squared sum of state variables
    return state[L];
  } );


  // Define a simple input update method
  ocp.setInputUpdate(
    [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input, double currentTime, void* userData ) {
      VECTOR<double, InputSize> updatedInput = { 0.0, 0.90 };
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
    inputConstraints.setConstant( -1.0 );
    inputConstraints[F] = -2.0;
    return inputConstraints;
  } );
  // Define a inputs constraints method
  ocp.setUpdateInputUpperBounds( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input ) {
    VECTOR<double, InputSize> inputConstraints;
    inputConstraints.setConstant( 1.0 );
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
  // user defined function for all constreaintzs


  /*
 OptiNLC_Solver<double, InputSize, StateSize, ConstraintsSize,ControlPoints> _solver(ocp);
  _solver.solve(0.0, initialState, initialInput);
 auto opt_x =  _solver.get_optimal_states();
 auto time = _solver.getTime();
 auto opt_u = _solver.get_optimal_inputs();
  // Open a file for writing

      std::ofstream dataFile("eigen_data.txt");

      // Check if the file is open
      if (dataFile.is_open()) {
          for (int i = 0; i < time.size(); ++i) {
              dataFile << time[i]<< " " << opt_x[StateSize*i + X] << " " << opt_x[StateSize*i + Y]
                      << " " << opt_x[StateSize*i + PSI]
                      << " " << opt_x[StateSize*i + V]
                      << " " << opt_x[StateSize*i + DELTA]<< " " << opt_u[InputSize*i] <<std::endl;
          }
          dataFile.close();
          std::cout << "Data saved to eigen_data.txt" << std::endl;
      } else {
          std::cerr << "Unable to open file for writing." << std::endl;
      }

*/
  REQUIRE( false );
}
