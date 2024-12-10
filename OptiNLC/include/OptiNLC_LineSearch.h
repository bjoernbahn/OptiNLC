/*
 * OptiNLC_LineSearch.h
 *
 * This file defines the LineSearch class, which implements various line search algorithms for optimization problems within the OptiNLC
 * toolbox.
 *
 * Author: Reza Dariani
 * Email: reza.dariani@dlr.de
 */
#pragma once
#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>

#include "OptiNLC_Engine.h"
#include "OptiNLC_Solver.h"

/**
 * @brief Implements line search algorithms for optimization problems in the OptiNLC toolbox.
 *
 * @tparam Scalar The data type for numerical operations.
 * @tparam InputSize Dimension of the input vector.
 * @tparam StateSize Dimension of the state vector.
 * @tparam ConstraintSize Dimension of the constraint vector.
 * @tparam NumberOfControlPoints Number of control points.
 */
template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
class OptiNLC_LineSearch
{
public:

  using OptProblem = OptiNLC_Engine<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>;
  using TimeVector = VECTOR<Scalar, NumberOfControlPoints + 1>;

  OptiNLC_LineSearch( std::shared_ptr<OptProblem> OP, double perturbation = 1e-12 ) :
    OP_( std::move( OP ) ),
    perturbation( perturbation ),
    resetHessian( false )
  {}

  double max( const VECTOR<Scalar, InputSize>& u, const VECTOR<Scalar, InputSize>& direction, const VECTOR<Scalar, InputSize>& lb,
              const VECTOR<Scalar, InputSize>& ub ) const;

  double raydanLineSearch( const VECTOR<Scalar, StateSize>& x, const VECTOR<Scalar, InputSize>& u,
                           const VECTOR<Scalar, InputSize>& direction, double initialStepSize = 1.0, double alpha = 0.5, double beta = 0.5,
                           int maxIterations = 50 );

  double directionalDerivative( const VECTOR<Scalar, StateSize>& x, const VECTOR<Scalar, InputSize>& u,
                                const VECTOR<Scalar, InputSize>& direction, const Eigen::MatrixXd& hessian,
                                const VECTOR<Scalar, InputSize>& gradient, double initialObjective, const TimeVector& time,
                                int maxIter = 100, double rho = 0.5, double eta = 0.25, double tau = 0.5 );

  bool checkFullStep( const VECTOR<Scalar, StateSize>& x, const VECTOR<Scalar, InputSize>& u, const VECTOR<Scalar, InputSize>& direction,
                      const VECTOR<Scalar, ConstraintSize>& lambda, const TimeVector& time, double cn, double tol );

  bool
  reset() const
  {
    return resetHessian;
  }

  void
  deactivateReset()
  {
    resetHessian = false;
  }


private:

  std::shared_ptr<OptProblem> OP_;
  double                      perturbation;
  bool                        resetHessian;
};

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
bool
OptiNLC_LineSearch<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::checkFullStep(
  const VECTOR<Scalar, StateSize>& x, const VECTOR<Scalar, InputSize>& u, const VECTOR<Scalar, InputSize>& direction,
  const VECTOR<Scalar, ConstraintSize>& lambda, const TimeVector& time, double cn, double tol )
{
  VECTOR<Scalar, InputSize> newU = u + direction;
  OP_->evaluate( x, newU, time, perturbation, 3 );

  auto   constr = OP_->getConstraints();
  double cnNew  = std::max( ( OP_->getL() - constr ).maxCoeff(), ( constr - OP_->getU() ).maxCoeff() );

  auto                      gradObj      = OP_->getGradient();
  auto                      jacobian     = OP_->getJacobian();
  VECTOR<Scalar, InputSize> gradLagrange = gradObj;

  VECTOR<Scalar, ConstraintSize - InputSize>    lambdaHead;
  VECTOR<Scalar, InputSize>                     lambdaTail;
  MATRIX<ConstraintSize - InputSize, InputSize> jacobianPart;

  std::memcpy( &lambdaHead[0], &lambda[0], sizeof( lambdaHead ) );
  std::memcpy( &lambdaTail[0], &lambda[ConstraintSize - InputSize], sizeof( lambdaTail ) );

  for( int i = 0; i < ConstraintSize - InputSize; ++i )
  {
    for( int j = 0; j < InputSize; ++j )
    {
      jacobianPart[i][j] = jacobian[i][j];
    }
  }

  gradLagrange -= ( lambdaHead * jacobianPart );
  gradLagrange -= lambdaTail;

  double lgn        = gradLagrange.norm_inf();
  double lambdaNorm = lambda.norm_inf();
  double tolNew     = lgn / ( 1.0 + lambdaNorm );

  return ( cnNew + tolNew < 0.999999 * ( cn + tol ) );
}

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
double
OptiNLC_LineSearch<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::directionalDerivative(
  const VECTOR<Scalar, StateSize>& x, const VECTOR<Scalar, InputSize>& u, const VECTOR<Scalar, InputSize>& direction,
  const Eigen::MatrixXd& hessian, const VECTOR<Scalar, InputSize>& gradient, double initialObjective, const TimeVector& time, int maxIter,
  double rho, double eta, double tau )
{
  double stepSize               = 1.0;
  double constraintViolation_L1 = OP_->L1ConstrainsViolation( x, u, time );

  Eigen::VectorXd objGrad = gradient.getDataAsEigen();
  Eigen::VectorXd dir     = direction.getDataAsEigen();

  double mu               = ( objGrad.dot( dir ) + 0.5 * dir.dot( hessian * dir ) ) / ( ( 1 - rho ) * constraintViolation_L1 );
  double directionalDeriv = objGrad.dot( dir ) - mu * constraintViolation_L1;
  double merit            = initialObjective + mu * constraintViolation_L1;

  double alpha     = 1.0;
  double tolerance = 1e-6;

  for( int i = 1; i <= maxIter; ++i )
  {
    VECTOR<Scalar, InputSize> newU = u + alpha * direction;
    OP_->evaluate( x, newU, time, perturbation, 0 );

    double newObjective = OP_->getObjectiveFunctionValue();
    double newMerit     = newObjective + mu * constraintViolation_L1;

    if( newMerit <= merit + alpha * eta * directionalDeriv || std::abs( alpha * directionalDeriv ) < tolerance )
    {
      return alpha;
    }

    alpha *= tau;
  }

  resetHessian = true;
  return alpha;
}

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
double
OptiNLC_LineSearch<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::max( const VECTOR<Scalar, InputSize>& u,
                                                                                              const VECTOR<Scalar, InputSize>& direction,
                                                                                              const VECTOR<Scalar, InputSize>& lb,
                                                                                              const VECTOR<Scalar, InputSize>& ub ) const
{
  double step = std::numeric_limits<Scalar>::infinity(); // Start with an infinite step size.

  for( int i = 0; i < InputSize; i++ )
  {
    if( direction[i] > 0.0 ) // Moving in the positive direction.
    {
      if( u[i] < ub[i] ) // Check if we're inside the upper bound
        step = std::min( step, ( ub[i] - u[i] ) / direction[i] );
    }
    else if( direction[i] < 0.0 ) // Moving in the negative direction.
    {
      if( u[i] > lb[i] ) // Check if we're inside the lower bound
        step = std::min( step, ( lb[i] - u[i] ) / direction[i] );
    }
    // If direction[i] == 0.0, skip to the next element as it doesn't affect the step size.
  }

  return step;
}
