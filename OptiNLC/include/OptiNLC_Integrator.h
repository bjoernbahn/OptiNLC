#pragma once

#include <functional>
#include <iostream>

#include "OptiNLC_Data.h"

/**
 * @brief Integrator class for numerically solving dynamic models.
 *
 * This class provides numerical integration using the Runge-Kutta 4th order (RK4) method. It is designed to work with
 * user-defined dynamic models within the OptiNLC toolbox.
 *
 * @tparam Scalar The data type for numerical operations.
 * @tparam InputSize The dimension of the input vector.
 * @tparam StateSize The dimension of the state vector.
 * @tparam ConstraintSize The dimension of the constraints vector (unused in integration but kept for consistency).
 * @tparam NumberOfControlPoints The number of control points (unused in integration but kept for consistency).
 */
template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
class OptiNLC_Integrator
{
public:

  using StateVector          = VECTOR<Scalar, StateSize>;
  using InputVector          = VECTOR<Scalar, InputSize>;
  using DynamicModelFunction = std::function<void( const StateVector&, const InputVector&, StateVector&, double, void* )>;

  /**
   * @brief Constructor for the integrator.
   * @param step_size The time step size for numerical integration.
   */
  explicit OptiNLC_Integrator( Scalar step_size ) noexcept :
    step_size_( step_size )
  {}

  /**
   * @brief Sets the dynamic model function for integration.
   * @param dynamic_model A user-defined dynamic model function.
   */
  void
  setDynamicModel( const DynamicModelFunction& dynamic_model ) noexcept
  {
    dynamic_model_ = dynamic_model;
  }

  /**
   * @brief Sets the integration time step.
   * @param time_step The time step size.
   */
  void
  setTimeStep( Scalar time_step ) noexcept
  {
    step_size_ = time_step;
  }

  /**
   * @brief Performs numerical integration using the Runge-Kutta 4th order (RK4) method.
   *
   * @param state The current state vector.
   * @param input The current input vector.
   * @param current_time The current simulation time.
   * @param user_data Optional user-defined data passed to the dynamic model.
   * @return The updated state vector after one time step.
   */
  StateVector
  integrateRK4( const StateVector& state, const InputVector& input, double current_time, void* user_data = nullptr ) const
  {
    Scalar      h_half = step_size_ / 2.0;
    StateVector k1, k2, k3, k4, next_state, temp_state;

    // Compute k1
    dynamic_model_( state, input, k1, current_time, user_data );

    // Compute k2
    for( int i = 0; i < StateSize; ++i )
    {
      temp_state[i] = state[i] + h_half * k1[i];
    }
    dynamic_model_( temp_state, input, k2, current_time + h_half, user_data );

    // Compute k3
    for( int i = 0; i < StateSize; ++i )
    {
      temp_state[i] = state[i] + h_half * k2[i];
    }
    dynamic_model_( temp_state, input, k3, current_time + h_half, user_data );

    // Compute k4
    for( int i = 0; i < StateSize; ++i )
    {
      temp_state[i] = state[i] + step_size_ * k3[i];
    }
    dynamic_model_( temp_state, input, k4, current_time + step_size_, user_data );

    // Update state using RK4 formula
    for( int i = 0; i < StateSize; ++i )
    {
      next_state[i] = state[i] + ( step_size_ / 6.0 ) * ( k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i] );
    }

    return next_state;
  }

private:

  DynamicModelFunction dynamic_model_; ///< User-defined dynamic model function.
  Scalar               step_size_;     ///< Time step size for numerical integration.
};
