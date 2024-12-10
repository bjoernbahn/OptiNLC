#pragma once
#include <functional>
#include <iomanip>
#include <memory>

#include "OptiNLC_Data.h"
#include "OptiNLC_Integrator.h"
#include "OptiNLC_Options.h"

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
class OptiNLC_OCP
{
public:

  using InputVector           = VECTOR<Scalar, InputSize>;
  using StateVector           = VECTOR<Scalar, StateSize>;
  using IntegratedStateVector = VECTOR<Scalar, StateSize*( NumberOfControlPoints + 1 )>;
  using InputVectorHorizon    = VECTOR<Scalar, InputSize*( NumberOfControlPoints )>;
  ;
  using ConstraintVector  = VECTOR<Scalar, ConstraintSize>;
  using ConstraintHorizon = VECTOR<Scalar, ConstraintSize*( NumberOfControlPoints )>;
  using TimeVector        = VECTOR<Scalar, NumberOfControlPoints + 1>;
  using DynamicModelFunction
    = std::function<void( const StateVector&, const InputVector&, StateVector&, double current_time, void* userData )>;
  using ObjectiveFunction = std::function<Scalar( const StateVector&, const InputVector&, double current_time )>;
  using UpdateInputMethod = std::function<InputVector( const StateVector&, const InputVector&, double current_time, void* userData )>;
  using StateLowerBoundsFunction       = std::function<StateVector( const StateVector&, const InputVector& )>;
  using StateUpperBoundsFunctions      = std::function<StateVector( const StateVector&, const InputVector& )>;
  using InputLowerBoundsFunction       = std::function<InputVector( const StateVector&, const InputVector& )>;
  using InputUpperBoundsFunction       = std::function<InputVector( const StateVector&, const InputVector& )>;
  using FunctionConstraints            = std::function<ConstraintVector( const StateVector&, const InputVector& )>;
  using FunctionConstraintsLowerBounds = std::function<ConstraintVector( const StateVector&, const InputVector& )>;
  using FunctionConstraintsUpperBounds = std::function<ConstraintVector( const StateVector&, const InputVector& )>;
  OptiNLC_Options* options;

  OptiNLC_OCP( OptiNLC_Options* options ) :
    num_steps_( NumberOfControlPoints )
  {

    this->options           = options;
    time_step_              = options->timeStep;
    intermediateIntegration = options->intermediateIntegration;
    numerical_integrator_   = std::make_shared<OptiNLC_Integrator<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>>(
      time_step_ );
  }

  ~OptiNLC_OCP()
  {
    // Clean up the NumericalIntegration object using delete
  }

  int
  getIntermediateIntegration()
  {
    return intermediateIntegration;
  }

  // Pass the dynamic model function or functor
  void
  setDynamicModel( const DynamicModelFunction& dynamic_model )
  {
    dynamic_model_ = dynamic_model;
    numerical_integrator_->setDynamicModel( dynamic_model_ );
  }

  // Set the user-defined objective function
  void
  setObjectiveFunction( const ObjectiveFunction& objective_function )
  {
    objective_function_ = objective_function;
  }

  // Set the user-defined input update method
  void
  setInputUpdate( const UpdateInputMethod& update_input_method )
  {
    update_input_method_ = update_input_method;
  }

  void
  setUpdateStateLowerBounds( const StateLowerBoundsFunction& update_state_lower_bounds )
  {
    update_state_lower_bounds_ = update_state_lower_bounds;
  }

  StateVector
  getStateLowerBounds( const StateVector& states, const InputVector& inputs )
  {
    return update_state_lower_bounds_( states, inputs );
  }

  void
  setUpdateStateUpperBounds( const StateUpperBoundsFunctions& update_state_upper_bounds )
  {
    update_state_upper_bounds_ = update_state_upper_bounds;
  }

  StateVector
  getStateUpperBounds( const StateVector& states, const InputVector& inputs )
  {
    return update_state_upper_bounds_( states, inputs );
  }

  void
  setUpdateInputLowerBounds( const InputLowerBoundsFunction& update_input_lower_bounds )
  {
    update_input_lower_bounds_ = update_input_lower_bounds;
  }

  InputVector
  getInputLowerBounds( const StateVector& states, const InputVector& inputs )
  {
    return update_input_lower_bounds_( states, inputs );
  }

  void
  setUpdateInputUpperBounds( const InputUpperBoundsFunction& update_input_upper_bounds )
  {
    update_input_upper_bounds_ = update_input_upper_bounds;
  }

  InputVector
  getInputUpperBounds( const StateVector& states, const InputVector& inputs )
  {
    return update_input_upper_bounds_( states, inputs );
  }

  void
  setUpdateFunctionConstraints( const FunctionConstraints& update_function_constraints )
  {
    update_function_constraints_ = update_function_constraints;
  }

  ConstraintVector
  computeFunctionConstraints( const StateVector& states, const InputVector& inputs )
  {
    return update_function_constraints_( states, inputs );
  }

  void
  setUpdateFunctionConstraintsLowerBounds( const FunctionConstraintsLowerBounds& update_function_constraints_lower_bounds )
  {
    update_function_constraints_lower_bounds_ = update_function_constraints_lower_bounds;
  }

  ConstraintVector
  getFunctionConstraintsLowerBounds( const StateVector& states, const InputVector& inputs )
  {
    return update_function_constraints_lower_bounds_( states, inputs );
  }

  void
  setUpdateFunctionConstraintsUpperBounds( const FunctionConstraintsUpperBounds& update_function_constraints_upper_bounds )
  {
    update_function_constraints_upper_bounds_ = update_function_constraints_upper_bounds;
  }

  ConstraintVector
  getFunctionConstraintsUpperBounds( const StateVector& states, const InputVector& inputs )
  {
    return update_function_constraints_upper_bounds_( states, inputs );
  }

  InputVector
  computeInputUpdate( const StateVector& states, const InputVector& inputs, double current_time, void* userData = nullptr )
  {
    return update_input_method_( states, inputs, current_time, userData );
  }

  ConstraintVector
  evaluateConstraintFunction( const StateVector& states, const InputVector& inputs )
  {
    return update_function_constraints_( states, inputs );
  }

  ConstraintHorizon
  evaluateConstraintsHorizon( const IntegratedStateVector& integratedStates, const InputVectorHorizon& inputHorizon )
  {
    ConstraintHorizon constraintHorizon;
    InputVector       tmp_input;
    StateVector       tmp_state;
    for( int i = 0; i < NumberOfControlPoints; ++i )
    {
      tmp_state.set( &integratedStates[i * StateSize] );
      tmp_input.set( &inputHorizon[i * InputSize] );
      constraintHorizon.setSubVector( evaluateConstraintFunction( tmp_state, tmp_input ), i * ConstraintSize );
    }
    return constraintHorizon;
  }

  // Compute the objective value based on the user-defined function
  Scalar
  computeObjective( const StateVector& states, const InputVector& inputs, double currentTime )
  {
    if( objective_function_ )
    {
      return objective_function_( states, inputs, currentTime );
    }
    // Default to zero if no objective function is defined
    return Scalar( 0 );
  }

  void
  createTimeVector( Scalar current_time = 0.0 )
  {
    time_vector.clear();
    for( int i = 0; i <= num_steps_; ++i )
    {
      time_vector.set( current_time + static_cast<Scalar>( i ) * time_step_, i );
    }
  }

  TimeVector
  getTime()
  {
    return time_vector;
  }

  StateVector
  integrateDynamicModel( const StateVector& states, const InputVector& inputs, const double current_time )
  {
    return numerical_integrator_->integrateRK4( states, inputs, current_time );
  }

  void
  setIntegratorTimeStep( double ts )
  {
    numerical_integrator_->setTimeStep( ts );
  }

  IntegratedStateVector
  integrate( const StateVector& x, const InputVectorHorizon& inputs, const TimeVector& time )
  {

    IntegratedStateVector integratedStates;
    integratedStates.setSubVector( x, 0 );
    StateVector current_state;
    InputVector current_input;
    double      current_time;

    double dt_in = time_step_ / intermediateIntegration;
    // intermediate integration
    numerical_integrator_->setTimeStep( time_step_ / intermediateIntegration );
    current_state.set( x );
    for( int i = 1; i <= NumberOfControlPoints; ++i )
    {
      current_time = time[i - 1];
      current_input.set( &inputs[( i - 1 ) * InputSize] );
      for( int j = 0; j < intermediateIntegration; j++ )
      {
        current_state.set( numerical_integrator_->integrateRK4( current_state, current_input, current_time + j * dt_in ) );
      }
      integratedStates.setSubVector( current_state, i * StateSize );
    }

    return integratedStates;
  }

  double
  getTimeStep()
  {
    return time_step_;
  }

private:

  Scalar                         time_step_;
  int                            intermediateIntegration;
  int                            num_steps_;
  TimeVector                     time_vector;
  DynamicModelFunction           dynamic_model_;
  ObjectiveFunction              objective_function_;
  UpdateInputMethod              update_input_method_;
  StateLowerBoundsFunction       update_state_lower_bounds_;
  StateUpperBoundsFunctions      update_state_upper_bounds_;
  InputLowerBoundsFunction       update_input_lower_bounds_;
  InputUpperBoundsFunction       update_input_upper_bounds_;
  FunctionConstraints            update_function_constraints_;
  FunctionConstraintsLowerBounds update_function_constraints_lower_bounds_;
  FunctionConstraintsUpperBounds update_function_constraints_upper_bounds_;

  std::shared_ptr<OptiNLC_Integrator<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>> numerical_integrator_;
};
