#pragma once

#include <chrono>
#include <string>

/**
 * @brief Configuration options for the OptiNLC toolbox.
 *
 * Provides parameters to control optimization, SQP solver, and OSQP solver.
 */
struct OptiNLC_Options
{
  // General optimization parameters
  double timeStep                = 0.1;    ///< Time step for the optimization
  int    intermediateIntegration = 10;     ///< Intermediate integration steps
  double perturbation            = 1.e-12; ///< Perturbation value for numerical differentiation
  double convergenceThreshold    = 1.e-10; ///< Convergence threshold for optimization
  double accuracy                = 1e-6;   ///< Accuracy for SQP optimization
  double timeLimit               = 0.500;  ///< Time limit for the OptiNLC solver
  int    maxNumberOfIteration    = 100;    ///< Maximum number of iterations

  // OSQP solver options
  bool   OSQP_verbose       = false; ///< Verbosity option for OSQP solver
  bool   OSQP_warm_starting = true;  ///< Warm-starting option for OSQP solver
  int    OSQP_max_iter      = 100;   ///< Maximum number of iterations for OSQP solver
  double OSQP_time_limit    = 0.010; ///< Time limit for OSQP solver
  double OSQP_eps_abs       = 1e-5;  ///< Absolute tolerance for OSQP solver
  double OSQP_eps_rel       = 1e-5;  ///< Relative tolerance for OSQP solver

  double SQP_ACC = 1e-6; // Accuracy for the Sequential Quadratic Programming (SQP) method

  double OptiNLC_ACC        = 1e-6; // Accuracy for the Sequential Quadratic Programming (SQP) method
  double OptiNLC_time_limit = 0.5;  // Time limit for OptiNLC solver

  bool debugPrint = false; // Print for debugging

  /**
   * @brief Loads default options for the OptiNLC toolbox.
   * This is already handled by default member initializers but kept for explicit usage.
   */
  void
  setDefaults()
  {
    timeStep                = 0.1;
    intermediateIntegration = 10;
    perturbation            = 1.e-12;
    convergenceThreshold    = 1.e-10;
    accuracy                = 1e-6;
    timeLimit               = 0.500;
    maxNumberOfIteration    = 100;

    OSQP_verbose       = false;
    OSQP_warm_starting = true;
    OSQP_max_iter      = 100;
    OSQP_time_limit    = 0.010;
    OSQP_eps_abs       = 1e-5;
    OSQP_eps_rel       = 1e-5;

    SQP_ACC            = 1e-6;
    OptiNLC_ACC        = 1e-6;
    OptiNLC_time_limit = 0.5;

    debugPrint = false;
  }
};
