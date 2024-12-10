#ifndef OPTINLC_NENGINE_H
#define OPTINLC_NENGINE_H

#include <atomic>
#include <functional>
#include <thread>

#include "OptiNLC_Data.h"

template<typename Scalar, int InputSize, int StateSize, int ConstraintSize, int NumberOfControlPoints>
class OptiNLC_Engine
{
public:

  using InputVectorHorizon    = VECTOR<Scalar, InputSize>;
  using InputVector           = VECTOR<Scalar, InputSize / NumberOfControlPoints>;
  using StateVector           = VECTOR<Scalar, StateSize>;
  using ConstraintVector      = VECTOR<Scalar, ConstraintSize>;
  using IntegratedStateVector = VECTOR<Scalar, StateSize*( NumberOfControlPoints + 1 )>;
  using TimeVector            = VECTOR<Scalar, NumberOfControlPoints + 1>;

  OptiNLC_Engine(
    const InputVectorHorizon& initialGuess, const ConstraintVector& l, const ConstraintVector& u,
    const std::function<double( const StateVector&, const InputVector&, double currentTime )>&                      objectiveFunction_ocp,
    const std::function<ConstraintVector( const IntegratedStateVector&, const InputVectorHorizon& )>&               constraintFunctions_ocp,
    const std::function<IntegratedStateVector( const StateVector&, const InputVectorHorizon&, const TimeVector& )>& integrator_ocp ) :
    initialGuess( initialGuess ),
    objectiveFunction_ocp( objectiveFunction_ocp ),
    constraintFunctions_ocp( constraintFunctions_ocp ),
    integrator_ocp( integrator_ocp ),
    l( l ),
    u( u )
  {

    std::cout << "\n********************************************************";
    std::cout << "\nOptimization Engine has created an optimization problem";
    std::cout << "\nNumber of Inputs:         " << InputSize;
    std::cout << "\nNumber of States:         " << StateSize;
    std::cout << "\nNumber of Constraints:    " << ConstraintSize;
    std::cout << "\nNumber of Control Points: " << NumberOfControlPoints;
    std::cout << "\n********************************************************\n";
    // std::cout<<"\nInitial Guess:\n"; initialGuess.print();
  }

  InputVectorHorizon
  getInitialGuess()
  {
    return initialGuess;
  }

  void
  evaluate( const StateVector& x, const InputVectorHorizon& u, const TimeVector& time, double perturbation = 1.e-12, int flag = 0 )
  {
    // ConstraintSequence and vector can be used now
    //  u is a vector from t0 --> tf
    //  x is only at t0
    const int             N = InputSize;
    IntegratedStateVector integratedStates;
    int NumOfInputVariables = N / NumberOfControlPoints; // number of variables origianly in vector u at t0, or any other t
    // integrate : initial state, input sequence
    integratedStates = integrator_ocp( x, u, time );
    StateVector lastState;
    lastState.set( &integratedStates[StateSize * NumberOfControlPoints] );
    InputVector lastInput;
    lastInput.set( &u[InputSize * ( NumberOfControlPoints - 1 )] );

    // if Flag == 0, only evaluate objective function
    // if flag == 1, only evaluate constraints
    // if flag == 2, objective function + gradient
    //  if flag == 3, all
    if( flag == 0 )
    {
      objFunction = evaluateObjectiveFunction( lastState, lastInput, time[NumberOfControlPoints - 1] );
    }
    else if( flag == 1 )
    {
      // auto start = std::chrono::high_resolution_clock::now();
      // Evaluate only the constraints
      // entire u , entire integrated states
      constraints.set( evaluateConstraintFunctions( integratedStates, u ) );
    }
    else if( flag == 2 )
    {
      // auto start = std::chrono::high_resolution_clock::now();
      objFunction      = evaluateObjectiveFunction( lastState, lastInput, time[NumberOfControlPoints - 1] );
      double obj_f_bck = objFunction;
      //-----------------create parallel inputs and perturbation
      PerturbationMatrix_.setIdentity( perturbation );
      for( int i = 0; i < InputSize; i++ )
      {
        ParallelInputs_.setCol( i, u );
      }
      PerturbedInputs_ = ParallelInputs_ + PerturbationMatrix_;


      for( int i = 0; i < InputSize; ++i )
      {
        currentPerturbedInput_    = PerturbedInputs_.getColumn( i );
        integratedStatesPerturbed = integrator_ocp( x, currentPerturbedInput_, time );
        lastState.set( &integratedStatesPerturbed[StateSize * NumberOfControlPoints] );
        gradient[i] = ( ( evaluateObjectiveFunction( lastState, lastInput, time[NumberOfControlPoints - 1] ) - obj_f_bck ) / perturbation );
      }
    }
    else if( flag == 3 )
    {

      objFunction      = evaluateObjectiveFunction( lastState, lastInput, time[NumberOfControlPoints - 1] );
      double obj_f_bck = objFunction;
      constraints.set( evaluateConstraintFunctions( integratedStates, u ) );
      constraints_bck.set( constraints );
      //-----------------create parallel inputs and perturbation

      PerturbationMatrix_.setIdentity( perturbation );
      for( int i = 0; i < InputSize; i++ )
      {
        ParallelInputs_.setCol( i, u );
      }
      PerturbedInputs_ = ParallelInputs_ + PerturbationMatrix_;
      if( true )
      {
        //-----------------------------------
        const int numThreads = std::thread::hardware_concurrency();
        auto      eval       = [&]( int start, int end ) {
          for( int i = start; i < end; ++i )
          {
            // std::cout<<"\nthread number: "<<i;
            VECTOR<double, InputSize> currentPerturbedInput      = PerturbedInputs_.getColumn( i );
            IntegratedStateVector     integratedStatesPerturbed_ = integrator_ocp( x, currentPerturbedInput, time );
            StateVector               lastState_;
            lastState_.set( &integratedStatesPerturbed_[StateSize * NumberOfControlPoints] );
            gradient[i] = ( ( evaluateObjectiveFunction( lastState_, lastInput, time[NumberOfControlPoints - 1] ) - obj_f_bck )
                            / perturbation );
            jacobian.setCol( i, ( evaluateConstraintFunctions( integratedStatesPerturbed_, currentPerturbedInput ) - constraints_bck )
                                             / perturbation );
          }
        };
        if( numThreads <= 1 || InputSize < numThreads )
        {
          // Single-threaded or not enough dimensions for parallelization
          eval( 0, InputSize );
        }
        else
        {
          // Multi-threaded parallelization
          std::vector<std::thread>      threads;
          std::vector<std::atomic<int>> threadRanges( numThreads );

          for( int i = 0; i < numThreads; ++i )
          {
            int start       = ( i * InputSize ) / numThreads;
            int end         = ( ( i + 1 ) * InputSize ) / numThreads;
            threadRanges[i] = start;
            threads.emplace_back( std::thread( eval, start, end ) );
          }
          for( int i = 0; i < numThreads; ++i )
          {
            threads[i].join();
          }
        }
        //-------------------------------------
      }
      else
      {
        for( int i = 0; i < InputSize; ++i )
        {
          currentPerturbedInput_    = PerturbedInputs_.getColumn( i );
          integratedStatesPerturbed = integrator_ocp( x, currentPerturbedInput_, time );
          lastState.set( &integratedStatesPerturbed[StateSize * NumberOfControlPoints] );
          gradient[i] = ( ( evaluateObjectiveFunction( lastState, lastInput, time[NumberOfControlPoints - 1] ) - obj_f_bck )
                          / perturbation );
          jacobian.setCol( i, ( evaluateConstraintFunctions( integratedStatesPerturbed, currentPerturbedInput_ ) - constraints_bck )
                                / perturbation );
        }
      }
    }
    else
    {
      // Handle invalid flag values
      std::cerr << "Invalid flag value. Please use 0, 1, 3, or 4." << std::endl;
      return;
    }
  }

  double
  evaluateObjectiveFunction( const StateVector& x, const InputVector& u, double currentTime ) const
  {
    return objectiveFunction_ocp( x, u, currentTime );
  }

  ConstraintVector
  evaluateConstraintFunctions( const IntegratedStateVector& x, const InputVectorHorizon& u ) const
  {
    return constraintFunctions_ocp( x, u );
  }

  double
  maxConstrainsViolation( const StateVector& _x, const InputVectorHorizon& _u, const TimeVector& time )
  {

    double                cn               = 0.0;
    IntegratedStateVector integratedStates = integrator_ocp( _x, _u, time );
    ConstraintVector      constr           = constraintFunctions_ocp( integratedStates, _u );
    cn                                     = fmax( cn, ( l - constr ).maxCoeff() );
    cn                                     = fmax( cn, ( constr - u ).maxCoeff() );
    return cn;
  }

  double
  L1ConstrainsViolation( const StateVector& _x, const InputVectorHorizon& _u, const TimeVector& time )
  {

    double                constraint_norm_l1 = std::numeric_limits<Scalar>::epsilon();
    IntegratedStateVector integratedStates   = integrator_ocp( _x, _u, time );
    ConstraintVector      constr             = constraintFunctions_ocp( integratedStates, _u );
    ConstraintVector      lower_deviations   = ( l - constr );
    lower_deviations.cwiseMax( 0.0 );
    ConstraintVector upper_deviations = ( constr - u );
    upper_deviations.cwiseMax( 0.0 );
    constraint_norm_l1 += lower_deviations.sum() + upper_deviations.sum();
    return constraint_norm_l1;
  }

  VECTOR<Scalar, InputSize>
  getGradient()
  {
    return gradient;
  }

  MATRIX<ConstraintSize, InputSize>
  getJacobian()
  {
    return jacobian;
  }

  double
  getObjectiveFunctionValue()
  {
    return objFunction;
  }

  ConstraintVector
  getConstraints()
  {
    return constraints;
  }

  ConstraintVector
  getL()
  {
    return l;
  }

  ConstraintVector
  getU()
  {
    return u;
  }

private:

  int                                         N;
  int                                         NC;
  InputVectorHorizon                          initialGuess;
  ConstraintVector                            l, u;
  VECTOR<Scalar, InputSize>                   gradient;
  VECTOR<double, InputSize>                   currentPerturbedInput_;
  IntegratedStateVector                       integratedStatesPerturbed;
  MATRIX<ConstraintSize, InputSize>           jacobian;
  MATRIX<InputSize, InputSize>                ParallelInputs_, PerturbationMatrix_, PerturbedInputs_;
  ConstraintVector                            constraints, constraints_bck;
  double                                      objFunction;
  std::function<double( const InputVector& )> objectiveFunction;
  std::function<double( const StateVector&, const InputVector&, double currentTime )> objectiveFunction_ocp;
  // std::function<double(const StateVector&, const InputVector&)> objectiveFunction_ocp;
  std::function<ConstraintVector( const InputVector& )>                                                    constraintFunctions;
  std::function<ConstraintVector( const IntegratedStateVector&, const InputVectorHorizon& )>               constraintFunctions_ocp;
  std::function<IntegratedStateVector( const StateVector&, const InputVectorHorizon&, const TimeVector& )> integrator_ocp;
};
#endif // OPTINLC_NENGINE_H
