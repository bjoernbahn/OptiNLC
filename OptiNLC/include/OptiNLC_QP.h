#pragma once


/*
 * OptiNLC_QP.h
 *
 * This file solves the QP problem for the OptiNLC toolbox.
 * OptiNLC is specialized for solving nonlinear optimal control problems using the line search Sequential Quadratic Programming (SQP)
 * method.
 *
 * Author: Reza Dariani
 * Email: reza.dariani@dlr.de
 */


#include <float.h>

#include <cstddef>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "osqp.h"
#include <eigen3/Eigen/Sparse>
#ifndef DFLOAT
typedef double c_float; /* for numerical values  */
#else
typedef float c_float; /* for numerical values  */
#endif

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
class OptiNLC_QP
{
public:

  OptiNLC_QP( OptiNLC_Options* options )
  {
    // Constructor code goes here (if any)
    debug_qp_time = 0.0;
    this->options = options;
  }

  struct SOLUTION
  {
    VECTOR<Scalar, InputSize>      primal;
    VECTOR<Scalar, ConstraintSize> dual;
    VECTOR<Scalar, ConstraintSize> constJ_p; // product of constraint Jacobian with p (direction) i.e. primal
    double                         run_time;
    int                            status;
  };

  SOLUTION
  solveQP( MATRIX<ConstraintSize, InputSize>& A, const Eigen::MatrixXd& P, const VECTOR<Scalar, InputSize>& q,
           const VECTOR<Scalar, ConstraintSize>& constr, const VECTOR<Scalar, ConstraintSize>& l,
           const VECTOR<Scalar, ConstraintSize>& u ) const
  {
    auto     start_0  = std::chrono::high_resolution_clock::now();
    OSQPInt  exitflag = 0;
    SOLUTION solution;
    solution.status = 11; // status initialized to unsolved

    // Create lower and upper bounds (lb, ub)
    const VECTOR<Scalar, ConstraintSize> lb = l - constr;
    const VECTOR<Scalar, ConstraintSize> ub = u - constr;

    // Convert to standard C++ arrays (no need for raw pointers)
    std::vector<double> tmp_q( q.data(), q.data() + InputSize );
    std::vector<double> tmp_l( lb.data(), lb.data() + ConstraintSize );
    std::vector<double> tmp_u( ub.data(), ub.data() + ConstraintSize );

    // Set up OSQP solver
    OSQPWorkspace* work     = nullptr;
    OSQPSolver*    solver   = nullptr;
    OSQPSettings*  settings = new OSQPSettings(); // Use `new` instead of `malloc`
    CompactArrays  P_ca     = get_ca( P );

    A.getSparsity();
    auto A_dense = A.sparse;

    OSQPCscMatrix* _P = new OSQPCscMatrix(); // Use `new` for OSQPCscMatrix
    OSQPCscMatrix* _A = new OSQPCscMatrix(); // Use `new` for OSQPCscMatrix

    csc_set_data( _P, InputSize, InputSize, P_ca.nnz, P_ca.x.data(), P_ca.i.data(), P_ca.o.data() );
    csc_set_data( _A, ConstraintSize, InputSize, A_dense.getNonZeros(), A_dense.valuePtr(), A_dense.innerIndexPtr(),
                  A_dense.outerIndexPtr() );

    osqp_set_default_settings( settings );
    settings->verbose       = options->OSQP_verbose;
    settings->max_iter      = options->OSQP_max_iter;
    settings->warm_starting = options->OSQP_warm_starting;
    settings->time_limit    = options->OSQP_time_limit;
    settings->eps_abs       = options->OSQP_eps_abs;
    settings->eps_rel       = options->OSQP_eps_rel;

    auto start_1 = std::chrono::high_resolution_clock::now();
    exitflag     = osqp_setup( &solver, _P, tmp_q.data(), _A, tmp_l.data(), tmp_u.data(), ConstraintSize, InputSize, settings );

    if( !exitflag )
      exitflag = osqp_solve( solver );

    auto start_2 = std::chrono::high_resolution_clock::now();

    // Check solver status
    if( solver->info->status_val == 9 )
    {
      std::cout << "\nsolver->info->status_val == 9\n";
      Eigen::MatrixXd P_copy = P; // Make a non-const copy
      P_copy.resize( InputSize, InputSize );
      P_copy.setIdentity();

      P_ca = get_ca( P_copy );
      csc_set_data( _P, InputSize, InputSize, P_ca.nnz, P_ca.x.data(), P_ca.i.data(), P_ca.o.data() );
      exitflag = osqp_setup( &solver, _P, tmp_q.data(), _A, tmp_l.data(), tmp_u.data(), ConstraintSize, InputSize, settings );
      if( solver->info->status_val == 9 )
        std::cout << "\nQP STATUS: problem non convex   ";
      if( solver->info->status_val == 1 )
        std::cout << "\nQP STATUS: solved";
      if( solver->info->status_val == 2 )
        std::cout << "\nQP STATUS: solved inaccurate";
      std::cout << "\n******************\n" << solver->info->status_val;
    }

    // Fill the solution
    for( int i = 0; i < InputSize; ++i )
    {
      solution.primal[i] = solver->solution->x[i];
    }
    for( int i = 0; i < ConstraintSize; ++i )
    {
      solution.dual[i] = solver->solution->y[i];
    }
    solution.constJ_p = ( A * solution.primal );
    solution.run_time = solver->info->run_time;

    osqp_cleanup( solver );

    // Clean up dynamically allocated memory
    delete _A;
    delete _P;
    delete settings;

    return solution;
  }

  struct CompactArrays
  {
    std::vector<c_float>       x;   // values
    std::vector<long long int> i;   // inner indices
    std::vector<long long int> o;   // outer indices
    int                        nnz; // number of non-zeros elements
    int                        nc;
    int                        nr;
  };

  static CompactArrays
  get_ca( Eigen::MatrixXd P )
  {
    CompactArrays ca;
    ca.nr                                   = P.rows();
    ca.nc                                   = P.cols();
    Eigen::SparseMatrix<double> P_spr       = P.sparseView( 0.00 );
    int                         NumNonZeros = P_spr.nonZeros();
    int                         innerSize   = P_spr.innerSize();
    int                         outerSize   = P_spr.outerSize();
    ca.nnz                                  = NumNonZeros;

    for( int i = 0; i < NumNonZeros; i++ )
    {
      ca.x.push_back( P_spr.valuePtr()[i] );
      ca.i.push_back( P_spr.innerIndexPtr()[i] );
    }
    for( int i = 0; i < outerSize; i++ )
    {
      ca.o.push_back( P_spr.outerIndexPtr()[i] );
    }
    ca.o.push_back( P_spr.outerIndexPtr()[outerSize] );
    return ca;
  }

  double
  getTotalTime()
  {
    return debug_qp_time;
  }

private:

  OptiNLC_Options* options;
  double           debug_qp_time;

  Eigen::MatrixXd
  makeUpperTriangular( const Eigen::MatrixXd& inputMatrix )
  {
    Eigen::MatrixXd upperTriangular = inputMatrix;

    int numRows = upperTriangular.rows();
    int numCols = upperTriangular.cols();

    for( int col = 0; col < numCols - 1; ++col )
    {
      for( int row = col + 1; row < numRows; ++row )
      {
        double factor               = upperTriangular( row, col ) / upperTriangular( col, col );
        upperTriangular.row( row ) -= factor * upperTriangular.row( col );
      }
    }

    return upperTriangular;
  }
};
