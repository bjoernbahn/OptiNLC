#pragma once

#include <chrono>
#include <iostream>
#include <memory>

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include "OptiNLC_Data.h"
#include "OptiNLC_LineSearch.h"
#include "OptiNLC_Options.h"
#include "OptiNLC_QP.h"

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
class OptiNLC_SQP
{
public:

  using InputVectorHorizon = VECTOR<Scalar, InputSize>;
  using StateVector        = VECTOR<Scalar, StateSize>;
  using OptProblem         = OptiNLC_Engine<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>;
  using TimeVector         = VECTOR<Scalar, NumberOfControlPoints + 1>;
  double obj_f;

  OptiNLC_SQP( OptProblem& op, OptiNLC_Options* options ) :
    Problem( op )
  {
    std::cerr << "\n********************************************************";
    std::cerr << "\nSequential Quadratic Problem has created an optimization problem";
    std::cerr << "\nNumber of Inputs:         " << InputSize;
    std::cerr << "\nNumber of States:         " << StateSize;
    std::cerr << "\nNumber of Constraints:    " << ConstraintSize;
    std::cerr << "\nNumber of Control Points: " << NumberOfControlPoints;
    std::cerr << "\n********************************************************\n";
    this->options                  = options;
    iteration                      = 1;
    previousObjectiveFunctionValue = std::numeric_limits<double>::infinity();
    tol_prev                       = std::numeric_limits<double>::infinity();
    cns_prev                       = std::numeric_limits<double>::infinity();
    perturbation                   = options->perturbation;
    objValueChange                 = 0.0;
    alpha                          = 1.0;
    convergenceThreshold           = options->convergenceThreshold;
    SQP_ACC                        = options->OptiNLC_ACC;
    maxNumberOfIteration           = options->maxNumberOfIteration;
    hessian.setZero( InputSize, InputSize );
    hessian.setIdentity();
    _QP = std::make_shared<OptiNLC_QP<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>>( options );
    _LS = std::make_shared<OptiNLC_LineSearch<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>>(
      std::make_shared<OptProblem>( Problem ), perturbation );

    // std::cout<<"\nNumber of Control Points: "<<NumberOfControlPoints;
  }

  ~OptiNLC_SQP() {}; //, const Eigen::VectorXd& x0) {

  void
  RUN( StateVector x0, TimeVector time )
  {
    auto startMainProcess = std::chrono::high_resolution_clock::now();
    bool Status           = false;
    DETERIORATION         = false;
    timeLimitReached      = false;
    deteriorationCounter  = 0;
    debug_ls_time         = 0.0;
    debug_qp_time         = 0.0;
    debug_hessian_time    = 0.0;
    debug_eval_time       = 0.0;
    totalTime             = 0.0;
    this->x0.set( x0 );
    std::cout << "\nRUN SQP\n";

    u0.set( Problem.getInitialGuess() );
    u.set( u0 );
    Problem.evaluate( x0, u, time, perturbation, 3 ); // complete evaluation
    obj_f                          = Problem.getObjectiveFunctionValue();
    localGradient                  = Problem.getGradient();
    localJacobian                  = Problem.getJacobian();
    objValueChange                 = std::abs( obj_f - previousObjectiveFunctionValue );
    convergence                    = objValueChange < convergenceThreshold;
    previousObjectiveFunctionValue = obj_f;
    totalTime += std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - startMainProcess ).count();
    // std::cout<<"\ntime milliseconds: "<<std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() -
    // start).count()<<"\n";
    for( iteration; iteration <= maxNumberOfIteration; iteration++ )
    {
      auto loopTic = std::chrono::high_resolution_clock::now();
      // std::cout<<"\ninteration: "<<iteration<<"\n";
      //----------------------------------projection to the boundaries
      // std::cout<<"\nprojection to the boundaries\n";
      // auto start = std::chrono::high_resolution_clock::now();
      VECTOR<Scalar, InputSize> lb_u, ub_u; // lower and upper boundary of the inputs
      std::memcpy( &lb_u[0], &Problem.getL()[0], sizeof( lb_u ) );
      std::memcpy( &ub_u[0], &Problem.getU()[0], sizeof( ub_u ) );
      lb_u.cwiseMin( ub_u );
      u.cwiseMax( lb_u );
      // std::cout<<"\ntime milliseconds: "<<std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() -
      // start).count()<<"\n";
      //----------------------------------------------------Solve QP
      // std::cout<<"\nSolve QP\n";
      auto start = std::chrono::high_resolution_clock::now();
      // start = std::chrono::high_resolution_clock::now();
      auto tmp_jac   = Problem.getJacobian();
      solution       = _QP->solveQP( tmp_jac, hessian, Problem.getGradient(), Problem.getConstraints(), Problem.getL(), Problem.getU() );
      debug_qp_time += std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - start ).count();
      // std::cout<<"\ntime milliseconds: "<<std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() -
      // start).count()<<"\n";
      P_LAMBDA.set( solution.dual );
      P_LAMBDA.set( P_LAMBDA - LAMBDA );
      //-------------------------------------------------- excecute maximum feasble step
      double max_step = _LS->max( u, solution.primal, lb_u, ub_u );
      // std::cout<<"\nMax step: "<<max_step;
      //-----------------------------------------------------------LINE SEARCH
      // std::cout<<"\nLine Search\n";
      start = std::chrono::high_resolution_clock::now();
      // alpha = _LS->raydanLineSearch(x0,u,solution.primal, std::max(1e-12,max_step));
      // alpha = _LS->ExactLineSearch(x0,u,solution.primal);
      alpha = _LS->directionalDerivative( x0, u, solution.primal, hessian, localGradient, obj_f, time );
      if( false ) //_LS->reset())
      {
        if( _LS->checkFullStep( x0, u, solution.primal, LAMBDA, time, cns_prev, tol_prev ) )
        {
          alpha = 1.000;
        }
      }
      //
      alpha          = std::min( alpha, 1.0 );
      debug_ls_time += std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - start ).count();
      // alpha = std::max(alpha, 1.e-8);
      // std::cout<<"\ntime milliseconds: "<<std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() -
      // start).count()<<"\n";
      //-----------------------------------------------------------UPDATE
      // std::cout<<"\n"<<alpha;
      u      = u + alpha * solution.primal;
      LAMBDA = LAMBDA + alpha * P_LAMBDA;
      //-----------------------------------------------------------Lagrange Gradient old
      calcLagrangeGradient( LAMBDA, localGradient, localJacobian, L_gradient, 0 );
      L_gradient_prev = L_gradient;
      //-----------------------------------------------------------EVALUATE
      // std::cout<<"\nStartof evaluate;\n";
      start = std::chrono::high_resolution_clock::now();
      Problem.evaluate( x0, u, time, perturbation, 3 ); // complete evaluation
      debug_eval_time += std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - start ).count();
      obj_f            = Problem.getObjectiveFunctionValue();
      localGradient    = Problem.getGradient();
      localJacobian    = Problem.getJacobian();
      objValueChange   = std::abs( obj_f - previousObjectiveFunctionValue );
      convergence      = objValueChange < convergenceThreshold;
      previousObjectiveFunctionValue = obj_f;

      //-----------------------------------------------------------Lagrange Gradient current
      calcLagrangeGradient( LAMBDA, localGradient, localJacobian, L_gradient, 1 );
      //-----------------------------------------------------------HESSIAN APPROXIMATION
      start = std::chrono::high_resolution_clock::now();
      if( !_LS->reset() )
      {

        hessian = calcBFGS( L_gradient, solution.primal, hessian, alpha );
        Eigen::LLT<Eigen::MatrixXd> llt( hessian );
        bool                        isPosDef = ( llt.info() != Eigen::NumericalIssue );
        if( !isPosDef )
        {
          hessian.setIdentity() * 1.e-3;
        }
        else
        {
          hessian = hessian.triangularView<Eigen::Upper>();
        }
      }
      else
      {
        hessian.setIdentity() * 1.e-3;
      }
      debug_hessian_time += std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - start ).count();
      //-----------------------------------------------------------CHECK CONVERGENCE
      // std::cout<<"\nobj_f: "<<obj_f;

      bool CONVERGENCE = termination( L_gradient, LAMBDA, u, Problem.getConstraints(), Problem.getL(), Problem.getU(), solution.primal,
                                      solution.dual, alpha, convergence );
      if( CONVERGENCE )
      {
        Status = true;
        std::cerr << greenCode << "\nConverged at iteration  : " << iteration;
        std::cerr << "\nObj function            : " << obj_f;
        std::cerr << "\nQP total time           : " << debug_qp_time;
        std::cerr << "\nline search total time  : " << debug_ls_time;
        std::cerr << "\nhessian calc total time : " << debug_hessian_time;
        std::cerr << "\nproblem eval total time : " << debug_eval_time;
        std::cerr << "\ntime milliseconds       : "
                  << std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - startMainProcess ).count()
                  << resetCode << "\n";
        break;
      }
      totalTime += std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - loopTic ).count();
      if( totalTime >= options->OptiNLC_time_limit * 1000. )
      {
        timeLimitReached = true;
        break;
      }
    } // for loop
    if( !Status )
    {
      std::cerr << redCode << "NOT Converged. Iteration  : " << iteration;
      if( timeLimitReached )
      {
        std::cerr << orangeCode << "\nTime limit reached       ";
      }
      std::cerr << redCode << "\nObj function            : " << obj_f;
      std::cerr << "\nQP total time           : " << debug_qp_time;
      std::cerr << "\nline search total time  : " << debug_ls_time;
      std::cerr << "\nhessian calc total time : " << debug_hessian_time;
      std::cerr << "\nproblem eval total time : " << debug_eval_time;
      Problem.evaluate( x0, u_best, time, perturbation, 0 ); // complete evaluation
      obj_f = Problem.getObjectiveFunctionValue();
      std::cerr << "\nObj function best backup: " << obj_f;
      std::cerr << "\ntime milliseconds       : "
                << std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - startMainProcess ).count()
                << resetCode << "\n";
    }
  }

  InputVectorHorizon
  get_u()
  {
    return u;
  }

private:

  double debug_ls_time, debug_qp_time, debug_hessian_time, debug_eval_time;

  void
  calcLagrangeGradient( const VECTOR<Scalar, ConstraintSize>& lambda, const VECTOR<Scalar, InputSize>& gradObj,
                        const MATRIX<ConstraintSize, InputSize>& jacobian, VECTOR<Scalar, InputSize>& gradLagrange, int flag )
  {

    // Objective gradient
    if( flag == 0 )
    {
      gradLagrange.set( gradObj );
    }
    else if( flag == 1 )
    {
      gradLagrange.set( gradObj - gradLagrange );
      // std::cout<<"\nFLAG 1\n";gradLagrange.print_T();
    }
    else
    {
      gradLagrange.clear();
    }

    VECTOR<Scalar, ConstraintSize - InputSize>    lambda_head;
    VECTOR<Scalar, InputSize>                     lambda_tail;
    MATRIX<ConstraintSize - InputSize, InputSize> jacobian_part;
    std::memcpy( &lambda_head[0], &lambda[0], sizeof( lambda_head ) );
    std::memcpy( &lambda_tail[0], &lambda[ConstraintSize - InputSize], sizeof( lambda_tail ) );

    for( int i = 0; i < ConstraintSize - InputSize; ++i )
    {
      for( int j = 0; j < InputSize; ++j )
      {
        jacobian_part[i][j] = jacobian[i][j];
      }
    }
    auto lambda_jac = ( lambda_head * jacobian_part );
    gradLagrange    = ( gradLagrange - lambda_jac );
    gradLagrange.set( gradLagrange - lambda_tail );
  }

  Eigen::MatrixXd
  calcBFGS( const VECTOR<Scalar, InputSize>& LGradient, const VECTOR<Scalar, InputSize>& primal, const Eigen::MatrixXd& hessian,
            double alpha )
  {

    double h1          = 0.0;
    double h2          = 0.0;
    double thetaPowell = 0.0;
    int    damped;
    auto   LGradient_     = LGradient.getDataAsEigen();
    auto   LGradient_back = LGradient_;
    auto   B              = hessian;
    auto   primal_        = primal.getDataAsEigen();
    auto   Bdelta         = B * primal_;
    // std::cout<<"\nB\n"<<B;
    // std::cout<<"\nprimal\n"<<primal_;
    // std::cout<<"\nBD\n"<<Bdelta;
    h1 = primal_.dot( Bdelta );

    // h1 = delta^T * B * delta
    // h2 = delta^T * gamma
    h2 = primal_.dot( LGradient_ );
    // std::cout<<"\nh: "<<h1<<"\t"<<h2;
    //  std::cout<<"\nLGradient_back\n"<<LGradient_back;
    /* Powell's damping strategy to maintain pos. def. (Nocedal/Wright p.537; SNOPT paper)
     * Interpolates between current approximation and unmodified BFGS */
    double hessDampFac = 0.2;
    if( h2 < hessDampFac * h1 / alpha && fabs( h1 - h2 ) > 1.0e-12 )
    { // At the first iteration h1 and h2 are equal due to COL scaling

      thetaPowell = ( 1.0 - hessDampFac ) * h1 / ( h1 - h2 );
      // Redefine gamma and h2 = delta^T * gamma
      h2             = 0.0;
      LGradient_back = thetaPowell * LGradient_back + ( 1.0 - thetaPowell ) * Bdelta;
      h2             = primal_.dot( LGradient_back );
      damped         = 1;
    }
    // B_k+1 = B_k - Bdelta * (Bdelta)^T / h1 + gamma * gamma^T / h2
    double myEps = 1.0e2 * 1.e-16;
    if( fabs( h1 ) < myEps || fabs( h2 ) < myEps )
    {
      std::cout << "\nThe Hessian matrix exhibits poor condition, prompting a reset";
      B.setIdentity() * 1e-3;
    }
    else
    {
      // std::cout<<"\nELSE*******************";
      // std::cout<<"\nh: "<<h1<<"\t"<<h2;
      // std::cout<<"\nLGradient_back\n"<<LGradient_back;
      B = B - Bdelta * Bdelta.transpose() / h1 + LGradient_back * LGradient_back.transpose() / h2;
    }
    // std::cout<<"\nHessian final\n"<<B;
    return B;
  }

  bool
  termination( const VECTOR<Scalar, InputSize>& lagrangeGradient, const VECTOR<Scalar, ConstraintSize>& lambda,
               const VECTOR<Scalar, InputSize>& u, const VECTOR<Scalar, ConstraintSize>& constr, const VECTOR<Scalar, ConstraintSize>& lb,
               const VECTOR<Scalar, ConstraintSize>& ub, const VECTOR<Scalar, InputSize>& direction,
               const VECTOR<Scalar, ConstraintSize>& dual, const double alpha, const bool convergence )
  {
    // std::cout<<"\nLG\n"; lagrangeGradient.print_T();
    double lgn = lagrangeGradient.norm_inf();
    // std::cout<<"\nlgn: "<<lgn;
    double lambda_n = lambda.norm_inf();
    double tol      = lgn / ( 1.0 + lambda_n );
    Scalar cn       = 0.0;
    if( ConstraintSize )
    {
      cn = fmax( cn, ( lb - constr ).maxCoeff() );
      cn = fmax( cn, ( constr - ub ).maxCoeff() );
    }
    double un  = u.norm_inf();
    double cns = cn / ( 1.0 + un );


    if( ( tol > tol_prev ) && ( cns > cns_prev ) )
    {
      // std::cout<<"\ndeterioration";
      deteriorationCounter++;
      if( deteriorationCounter > 5 )
      {
        DETERIORATION        = true;
        deteriorationCounter = 0;
      }
    }
    else
    {
      // std::cout<<"\nNormal";
      u_best               = u;
      LAMBDA_best          = LAMBDA;
      DETERIORATION        = false;
      deteriorationCounter = 0;
    }
    // std::cout << std::fixed << std::defaultfloat <<"\n"<<tol - tol_prev<<"\t"<<cns - cns_prev;
    // std::cout << std::fixed << std::defaultfloat <<"\n"<<tol<<"\t"<<cns;
    tol_prev = tol;
    cns_prev = cns;
    std::cout << "\033[32m";
    if( tol <= SQP_ACC && cns <= SQP_ACC )
    {

      std::cout << "\nOptimization problem has converged with " << SQP_ACC << "accuracy";
      std::cout << "\n" << tol << "\t" << cns;
      std::cout << "\033[0m";
      return true;
    }

    // Other approach based on step norm
    double primal_step_norm     = alpha * direction.norm_inf();
    double dual_step_norm       = alpha * dual.norm_inf();
    double constraints_criteria = cn;
    // std::cout<<"\n";
    std::cout << std::fixed << std::defaultfloat;

    // std::cout << std::setw(15)
    // <<alpha<<"\t"<<primal_step_norm<<"\t"<<dual_step_norm<<"\t"<<constraints_criteria<<"\t"<<tol<<"\t"<<convergence<<"\n";
    //   std::cout<<"\n"<<alpha<<"\t"<<primal_step_norm<<"\t"<<dual_step_norm<<"\t"<<constraints_criteria<<"\t"<<convergence<<"\t"<<tol<<"\t"<<cns;
    if( primal_step_norm < SQP_ACC && dual_step_norm < SQP_ACC && constraints_criteria < SQP_ACC && convergence && tol < 0.0010 )
    {
      std::cout << "\nOptimization problem has converged with " << SQP_ACC << " accuracy. [Step Norm]";
      std::cout << "\n"
                << alpha << "\t" << primal_step_norm << "\t" << dual_step_norm << "\t" << constraints_criteria << "\t" << convergence
                << "\t" << tol << "\t" << cns;
      std::cout << "\033[0m";
      return true;
    }
    std::cout << "\033[0m";
    return false;
  }

  bool
  hasGoodConditionNumber( const Eigen::MatrixXd& hessian, double threshold = 1e-8 )
  {


    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es( hessian );
    if( es.info() != Eigen::Success )
    {
      std::cerr << "Eigenvalue computation failed. Hessian may be singular or ill-conditioned." << std::endl;
      return false;
    }

    double maxEigenvalue = es.eigenvalues().maxCoeff();
    double minEigenvalue = es.eigenvalues().minCoeff();

    if( minEigenvalue < threshold )
    {
      std::cerr << "Hessian matrix is ill-conditioned (minimum eigenvalue < threshold)." << std::endl;
      return false;
    }

    double conditionNumber = maxEigenvalue / minEigenvalue;

    if( conditionNumber > 5.0 )
    {
      // std::cout << "\nHESSIAN RESET\n";
      return false;
      // You might perform an action here based on the high condition number
    }

    return ( conditionNumber < threshold );
  }

  OptProblem&                                                                                              Problem;
  OptiNLC_Options*                                                                                         options;
  std::shared_ptr<OptiNLC_QP<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>>         _QP;
  std::shared_ptr<OptiNLC_LineSearch<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>> _LS;

  typename OptiNLC_QP<Scalar, InputSize, StateSize, ConstraintSize, NumberOfControlPoints>::SOLUTION solution;


  InputVectorHorizon                u0;
  InputVectorHorizon                u, u_best;
  StateVector                       x0;
  unsigned int                      iteration;
  double                            previousObjectiveFunctionValue;
  double                            perturbation;
  double                            objValueChange;
  bool                              convergence;
  double                            convergenceThreshold;
  int                               maxNumberOfIteration;
  double                            alpha;
  Eigen::MatrixXd                   hessian;
  VECTOR<Scalar, InputSize>         localGradient;
  VECTOR<Scalar, InputSize>         L_gradient, L_gradient_prev;
  MATRIX<ConstraintSize, InputSize> localJacobian;
  VECTOR<Scalar, ConstraintSize>    LAMBDA, P_LAMBDA, LAMBDA_best;
  double                            tol_prev, cns_prev;
  double                            SQP_ACC;
  bool                              DETERIORATION;
  bool                              timeLimitReached;
  double                            totalTime;
  int                               deteriorationCounter;
  std::string                       orangeCode = "\033[38;5;208m";
  std::string                       redCode    = "\033[38;5;196m";
  std::string                       greenCode  = "\033[38;5;46m";
  std::string                       resetCode  = "\033[0m";
};
