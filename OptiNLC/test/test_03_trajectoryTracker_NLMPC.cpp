#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "OptiNLC_Data.h"
#include "OptiNLC_Options.h"
#include "OptiNLC_Solver.h"
#include <catch2/catch.hpp>

class ReadOfflineData
{
public:

  ReadOfflineData()
  {
    // You can handle exceptions here or remove the throw statement
    // depending on your desired behavior.
    // throw std::runtime_error("Born");
    std::cout << "\nCALLED\n";
  }

  struct PiecewisePolynomialVector
  {
    std::vector<double> coef1;
    std::vector<double> coef2;
  } pp_X, pp_Y, pp_PSI, pp_V, pp_D, pp_A;

  std::vector<double> Time;

  bool
  readTextDataOffline()
  {
    std::string filename_x     = "/tmp/OptiNLC/pp_x.txt";
    std::string filename_y     = "/tmp/OptiNLC/pp_y.txt";
    std::string filename_psi   = "/tmp/OptiNLC/pp_psi.txt";
    std::string filename_v     = "/tmp/OptiNLC/pp_v.txt";
    std::string filename_delta = "/tmp/OptiNLC/pp_delta.txt";
    std::string filename_a     = "/tmp/OptiNLC/pp_a.txt";
    std::string filename_t     = "/tmp/OptiNLC/time.txt";

    // Open the files
    std::ifstream file_x( filename_x );
    std::ifstream file_y( filename_y );
    std::ifstream file_psi( filename_psi );
    std::ifstream file_v( filename_v );
    std::ifstream file_delta( filename_delta );
    std::ifstream file_a( filename_a );
    std::ifstream file_t( filename_t );

    // Check if the files are open
    if( !file_x.is_open() || !file_y.is_open() || !file_psi.is_open() || !file_v.is_open() || !file_delta.is_open() || !file_a.is_open()
        || !file_t.is_open() )
    {
      std::cerr << "Error opening one or more files" << std::endl;
      return false;
    }

    double value;

    // Read data from the files and push it into the vectors
    while( file_x >> value )
    {
      pp_X.coef1.push_back( value );
      file_x >> value;
      pp_X.coef2.push_back( value );
    }
    file_x.close();
    value = 0.;
    while( file_y >> value )
    {
      pp_Y.coef1.push_back( value );
      file_y >> value;
      pp_Y.coef2.push_back( value );
    }
    file_y.close();
    value = 0.;
    while( file_psi >> value )
    {
      pp_PSI.coef1.push_back( value );
      file_psi >> value;
      pp_PSI.coef2.push_back( value );
    }
    file_psi.close();
    value = 0.;
    while( file_v >> value )
    {
      pp_V.coef1.push_back( value );
      file_v >> value;
      pp_V.coef2.push_back( value );
    }
    file_v.close();
    value = 0.;
    while( file_delta >> value )
    {
      pp_D.coef1.push_back( value );
      file_delta >> value;
      pp_D.coef2.push_back( value );
    }
    file_delta.close();
    value = 0.;
    while( file_a >> value )
    {
      pp_A.coef1.push_back( value );
      file_a >> value;
      pp_A.coef2.push_back( value );
    }
    file_a.close();

    value = 0.;
    while( file_t >> value )
    {
      Time.push_back( value );
    }
    file_t.close();
    // Similar code for file_y and file_psi

    return true;
  }

  int
  findIndex( double data, std::vector<double>& ref )
  {
    int index = -1;
    for( int i = 0; i < ref.size() - 1; i++ )
    {
      if( data >= ref[i] && data <= ref[i + 1] )
      {
        index = i;
        break;
      }
    }
    if( index == -1 )
      std::cout << "\nERROR\n";
    return index;
  }

  double
  interpolation( double data, int index, PiecewisePolynomialVector& pp )
  {
    return pp.coef1[index] * data + pp.coef2[index];
  }

  std::vector<double>
  generateData( double start, double end, double step )
  {


    // Calculate the number of elements in the sequence
    size_t numElements = static_cast<size_t>( ( end - start ) / step ) + 1;
    // Initialize the vector with the calculated number of elements
    std::vector<double> data( numElements );
    // Generate the sequence and store it in the vector
    for( size_t i = 0; i < numElements; ++i )
    {
      data[i] = start + i * step;
    }
    return data;
  }
};

TEST_CASE( "OCP Test NEW FORMAT Tracker NLMPC:" )
{
  std::cout << "OCP Test NEW FORMAT Tracker NLMPC:" << std::endl;


  ReadOfflineData test;
  bool            stausData = test.readTextDataOffline();
  if( !stausData )
  {
    std::cout << "\nOffline data are NOT loaded";
    return;
  }
  std::cout << "\n****** Offline data are loaded ******\n";
  // std::cout<<"\n"<<test.pp_X.coef1[5]<<"\t"<<test.pp_Y.coef1[5]<<"\t"<<test.pp_PSI.coef1[5]<<"\t"<<test.pp_V.coef1[5]<<"\t"<<test.pp_D.coef1[5]<<"\t"<<test.pp_A.coef1[5]<<"\t"<<test.Time[5];

  // auto hr_data = test.generateData(test.Time[0], test.Time.back(), 0.01);
  // for(int i=0; i<hr_data.size(); i++)
  // {
  //    int index = test.findIndex(hr_data[i],test.Time);
  //    std::cout<<"\n"<<test.interpolation(hr_data[i],index,test.pp_X);

  //}

  OptiNLC_Options options;
  options.setDefaults();

  options.intermediateIntegration = 2;
  options.OptiNLC_ACC             = 1.e-7;
  options.maxNumberOfIteration    = 500;
  options.OSQP_verbose            = false;
  options.OSQP_max_iter           = 400;
  options.OptiNLC_time_limit      = 0.05;
  options.OSQP_time_limit         = 0.1;

  enum STATES
  {
    X,
    Y,
    PSI,
    V,
    L
  };

  enum CONTROLS
  {
    DELTA,
    A
  };

  // Define the state and input sizes
  constexpr int StateSize       = 5;
  constexpr int InputSize       = 2;
  constexpr int ConstraintsSize = 0;
  constexpr int ControlPoints   = 10;
  double        simTime         = 0.5;
  options.timeStep              = ( simTime / ControlPoints );

  // Create an instance of the OCP class with the appropriate horizon and time step
  double current_time = test.Time[20];
  int    index        = test.findIndex( current_time, test.Time );
  // to be replaced with ego position
  double                    x_init       = test.interpolation( current_time, index, test.pp_X );
  double                    y_init       = test.interpolation( current_time, index, test.pp_Y );
  double                    psi_init     = test.interpolation( current_time, index, test.pp_PSI );
  double                    v_init       = test.interpolation( current_time, index, test.pp_V );
  double                    a_init       = test.interpolation( current_time, index, test.pp_A );
  double                    delta_init   = test.interpolation( current_time, index, test.pp_D );
  VECTOR<double, StateSize> initialState = { x_init, y_init + 0.1, psi_init, v_init, 0.0 };
  VECTOR<double, InputSize> initialInput = { delta_init, 0 };
  OptiNLC_OCP<double, InputSize, StateSize, ConstraintsSize, ControlPoints> ocp( &options );
  // Define a simple dynamic model as a lambda function


  ocp.setDynamicModel( [&test]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input,
                                VECTOR<double, StateSize>& derivative, double currentTime, void* userData ) {
    derivative[X]        = state[V] * cos( state[PSI] );
    derivative[Y]        = state[V] * sin( state[PSI] );
    derivative[PSI]      = state[V] * tan( input[DELTA] ) / 4.0;
    derivative[V]        = input[A];
    int    index         = test.findIndex( currentTime, test.Time );
    double reference_x   = test.interpolation( currentTime, index, test.pp_X );
    double reference_y   = test.interpolation( currentTime, index, test.pp_Y );
    double reference_psi = test.interpolation( currentTime, index, test.pp_PSI );
    double reference_v   = test.interpolation( currentTime, index, test.pp_V );
    double c_psi         = cos( reference_psi );
    double s_psi         = sin( reference_psi );
    double d_lat         = -s_psi * ( reference_x - state[X] ) + c_psi * ( reference_y - state[Y] );
    double d_long        = c_psi * ( reference_x - state[X] ) + s_psi * ( reference_y - state[Y] );
    derivative[L]        = ( d_lat * d_lat + d_long * d_long ) * 10.0 + 0.0 * ( input[A] * input[A] + input[DELTA] * input[DELTA] * 0.5 )
                  + 0.1 * ( state[V] - reference_v ) * ( state[V] - reference_v );
    // derivative[L] = (state[V]-reference_v)*(state[V]-reference_v)*0.0 + (state[Y]-reference_y)*(state[Y]-reference_y)
    // + (state[X]-reference_x)*(state[X]-reference_x) + (state[PSI]-reference_psi)*(state[PSI]-reference_psi);
  } );

  // Define an objective function as a lambda function, states and input at tf
  ocp.setObjectiveFunction( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input, double t_end ) {
    // Define your objective function here
    // For example, you can minimize the squared sum of state variables
    int    index         = test.findIndex( t_end, test.Time );
    double reference_x   = test.interpolation( t_end, index, test.pp_X );
    double reference_y   = test.interpolation( t_end, index, test.pp_Y );
    double reference_psi = test.interpolation( t_end, index, test.pp_PSI );
    return state[L] + ( ( state[X] - reference_x ) * ( state[X] - reference_x ) + ( state[Y] - reference_y ) * ( state[Y] - reference_y ) );
  } );


  // Define a simple input update method
  ocp.setInputUpdate(
    [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input, double currentTime, void* userData ) {
      int                       index        = test.findIndex( currentTime, test.Time );
      VECTOR<double, InputSize> updatedInput = { test.interpolation( currentTime, index, test.pp_D ),
                                                 test.interpolation( currentTime, index, test.pp_A ) };
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
    inputConstraints.setConstant( -3.0 );
    inputConstraints[DELTA] = -.70;
    return inputConstraints;
  } );
  // Define a inputs constraints method
  ocp.setUpdateInputUpperBounds( [&]( const VECTOR<double, StateSize>& state, const VECTOR<double, InputSize>& input ) {
    VECTOR<double, InputSize> inputConstraints;
    inputConstraints.setConstant( 3.0 );
    inputConstraints[DELTA] = 0.70;
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
  _solver.solve(current_time, initialState, initialInput);
 auto opt_x =  _solver.get_optimal_states();
 auto time = _solver.getTime();
 auto opt_u = _solver.get_optimal_inputs();
 double c_psi, s_psi;
 std::vector<double> x_ref, y_ref, v_ref, psi_ref, d_lat, d_long, d_v;
  for(int i=0; i<time.size(); i++)
  {
     int index = test.findIndex(time[i],test.Time);
     x_ref.push_back(test.interpolation(time[i],index,test.pp_X));
     y_ref.push_back(test.interpolation(time[i],index,test.pp_Y));
     v_ref.push_back(test.interpolation(time[i],index,test.pp_V));
     psi_ref.push_back(test.interpolation(time[i],index,test.pp_PSI));
     c_psi = cos(psi_ref.back());
     s_psi = sin(psi_ref.back());
     d_lat.push_back(-s_psi * (x_ref.back() - opt_x[StateSize*i + X] ) + c_psi * (y_ref.back() - opt_x[StateSize*i + Y]));
     d_long.push_back( c_psi * (x_ref.back() - opt_x[StateSize*i + X]) +  s_psi * (y_ref.back() - opt_x[StateSize*i + Y]));
     d_v.push_back(v_ref.back() - opt_x[StateSize*i + V]);
 }
  // Open a file for writing

      std::ofstream dataFile("eigen_data.txt");

      // Check if the file is open
      if (dataFile.is_open()) {
          for (int i = 0; i < time.size(); ++i) {
              dataFile << time[i]<< " "
                      << opt_x[StateSize*i + X]
                      << " " << opt_x[StateSize*i + Y]
                      << " " << opt_x[StateSize*i + PSI]
                      << " " << opt_x[StateSize*i + V]
                      << " " << opt_u[InputSize*i + DELTA]
                      << " " << opt_u[InputSize*i + A]
                      << " " << x_ref[i]
                      << " " << y_ref[i]
                      << " " << v_ref[i]
                      << " " << d_lat[i]
                      << " " << d_long[i]
                      << " " << d_v[i]
                      <<std::endl;
          }
          dataFile.close();
          std::cout << "Data saved to eigen_data.txt" << std::endl;
      } else {
          std::cerr << "Unable to open file for writing." << std::endl;
      }

*/
  REQUIRE( false );
}
