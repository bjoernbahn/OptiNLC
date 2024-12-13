#pragma once
#include <cstring>
#include <memory>

#include "OptiNLC_Data.h"
#include "OptiNLC_OCP.h"
#include "OptiNLC_Options.h"
#include "OptiNLC_QP.h"
#include "OptiNLC_SQP.h"

/**
 * @brief Solver for Nonlinear Optimal Control Problems using Sequential Quadratic Programming (SQP).
 *
 * @tparam Scalar Data type for numerical operations.
 * @tparam InputSize Dimension of the input vector.
 * @tparam StateSize Dimension of the state vector.
 * @tparam ConstraintSize Dimension of the constraint vector.
 * @tparam NumberOfControlPoints Number of control points.
 */
template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
class OptiNLC_Solver
{
public:

  using OCPType               = OptiNLC_OCP<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>;
  using StateVector           = typename OCPType::StateVector;
  using IntegratedStateVector = typename OCPType::IntegratedStateVector;
  using InputVector           = typename OCPType::InputVector;
  using InputVectorHorizon    = typename OCPType::InputVectorHorizon;
  using ConstraintVector      = typename OCPType::ConstraintVector;
  using ConstraintHorizon     = typename OCPType::ConstraintHorizon;
  using TimeVector            = typename OCPType::TimeVector;
  using BoundaryVector
    = VECTOR<Scalar, InputSize * NumberOfControlPoints + StateSize*( NumberOfControlPoints + 1 ) + ConstraintSize * NumberOfControlPoints>;

  explicit OptiNLC_Solver( OCPType& ocp ) :
    ocp_( std::make_shared<OCPType>( ocp ) ),
    options_( ocp_->options ),
    final_obj_function_( 0.0 )
  {}

  void solve( double currentTime, const StateVector& state, const InputVector& input );

  TimeVector
  getTime() const
  {
    return ocp_->getTime();
  }

  double
  get_final_objective_function() const
  {
    return final_obj_function_;
  }

  InputVectorHorizon
  get_optimal_inputs() const
  {
    return optimal_inputs_;
  }

  IntegratedStateVector
  get_optimal_states() const
  {
    return optimal_states_;
  }

private:

  std::shared_ptr<OCPType> ocp_;
  OptiNLC_Options*         options_;
  double                   final_obj_function_;
  InputVectorHorizon       optimal_inputs_;
  IntegratedStateVector    optimal_states_;

  InputVectorHorizon                                     input_lb_, input_ub_;
  IntegratedStateVector                                  state_lb_, state_ub_;
  VECTOR<Scalar, ConstraintSize * NumberOfControlPoints> constraint_lb_, constraint_ub_;

  void           evaluateBoundaries();
  BoundaryVector getL() const;
  BoundaryVector getU() const;
  BoundaryVector constraintsFunction( const IntegratedStateVector& integratedStates, const InputVectorHorizon& inputHorizon ) const;
  void           initializeAndSimulate( const StateVector& state, const InputVector& input );
};

// Implementation

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
void
OptiNLC_Solver<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::solve( double currentTime, const StateVector& state,
                                                                                            const InputVector& input )
{
  ocp_->createTimeVector( currentTime );
  initializeAndSimulate( state, input );
  evaluateBoundaries();

  auto objectiveFunctionLambda = [&]( const StateVector& x, const InputVector& u, double currentTime ) {
    return ocp_->computeObjective( x, u, currentTime );
  };
  auto constraintsFunctionLambda = [&]( const IntegratedStateVector& x, const InputVectorHorizon& u ) {
    return constraintsFunction( x, u );
  };
  auto integratorFunctionLambda = [&]( const StateVector& x, const InputVectorHorizon& u, const TimeVector& time ) {
    return ocp_->integrate( x, u, time );
  };

  constexpr int num_constraints = InputSize * NumberOfControlPoints + ConstraintSize * NumberOfControlPoints
                                + StateSize * ( NumberOfControlPoints + 1 );

  OptiNLC_Engine<Scalar, InputSize * NumberOfControlPoints, StateSize, num_constraints, NumberOfControlPoints> problem(
    optimal_inputs_, getL(), getU(), objectiveFunctionLambda, constraintsFunctionLambda, integratorFunctionLambda );


  OptiNLC_SQP<Scalar, InputSize * NumberOfControlPoints, StateSize, num_constraints, NumberOfControlPoints> sqp_solver( problem, options_ );
  sqp_solver.RUN( state, ocp_->getTime() );

  optimal_inputs_     = sqp_solver.get_u();
  optimal_states_     = ocp_->integrate( state, optimal_inputs_, ocp_->getTime() );
  final_obj_function_ = sqp_solver.obj_f;
}

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
void
OptiNLC_Solver<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::evaluateBoundaries()
{
  input_lb_.clear();
  input_ub_.clear();
  state_lb_.clear();
  state_ub_.clear();
  constraint_lb_.clear();
  constraint_ub_.clear();

  StateVector tmp_state;
  InputVector tmp_input;

  for( int i = 0; i <= NumberOfControlPoints; ++i )
  {
    tmp_state.set( &optimal_states_[i * StateSize] );
    tmp_input.set( &optimal_inputs_[i * InputSize] );

    if( i < NumberOfControlPoints )
    {
      input_lb_.setSubVector( ocp_->getInputLowerBounds( tmp_state, tmp_input ), i * InputSize );
      input_ub_.setSubVector( ocp_->getInputUpperBounds( tmp_state, tmp_input ), i * InputSize );
    }

    state_lb_.setSubVector( ocp_->getStateLowerBounds( tmp_state, tmp_input ), i * StateSize );
    state_ub_.setSubVector( ocp_->getStateUpperBounds( tmp_state, tmp_input ), i * StateSize );

    if( ConstraintSize )
    {
      constraint_lb_.setSubVector( ocp_->getFunctionConstraintsLowerBounds( tmp_state, tmp_input ), i * ConstraintSize );
      constraint_ub_.setSubVector( ocp_->getFunctionConstraintsUpperBounds( tmp_state, tmp_input ), i * ConstraintSize );
    }
  }
}

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
typename OptiNLC_Solver<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::BoundaryVector
OptiNLC_Solver<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::getL() const
{
  BoundaryVector L;
  L.setSubVector( input_lb_, 0 );
  L.setSubVector( state_lb_, input_lb_.size() );
  L.setSubVector( constraint_lb_, input_lb_.size() + state_lb_.size() );
  return L;
}

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
typename OptiNLC_Solver<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::BoundaryVector
OptiNLC_Solver<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::getU() const
{
  BoundaryVector U;
  U.setSubVector( input_ub_, 0 );
  U.setSubVector( state_ub_, input_ub_.size() );
  U.setSubVector( constraint_ub_, input_ub_.size() + state_ub_.size() );
  return U;
}

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
typename OptiNLC_Solver<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::BoundaryVector
OptiNLC_Solver<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::constraintsFunction(
  const IntegratedStateVector& integratedStates, const InputVectorHorizon& inputHorizon ) const
{
  BoundaryVector    result;
  ConstraintHorizon constraintHorizon;
  if( ConstraintSize )
  {
    constraintHorizon = ocp_->evaluateConstraintsHorizon( integratedStates, inputHorizon );
  }

  std::memcpy( &result[0], &inputHorizon[0], sizeof( inputHorizon ) );
  std::memcpy( &result[inputHorizon.size()], &integratedStates[0], sizeof( integratedStates ) );
  if( ConstraintSize )
  {
    std::memcpy( &result[inputHorizon.size() + integratedStates.size()], &constraintHorizon[0], sizeof( constraintHorizon ) );
  }
  return result;
}

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
void
OptiNLC_Solver<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::initializeAndSimulate( const StateVector& state,
                                                                                                            const InputVector& input )
{
  optimal_states_ = ocp_->integrate( state, InputVectorHorizon{}, ocp_->getTime() );
}
