#ifndef differential_equations_h
#define differential_equations_h



// compile flags are neccessary:
//   -I flag tells g++ where to look for Eigen
//USE_COMPILE_FLAG -I${PROJECTDIR}/include
//USE_COMPILE_FLAG -Wno-shadow -Wno-conversion


#include <stdexcept>
#include <list>
#include <vector>
#include <string>
#include <tuple>

#include <Eigen/Dense>
#include "automatic_derivatives/Fad.h"
#include "Trace.h"
#include "constants.h"



typedef Fad<double, 1> FDouble;

//#define USE_SPARSE_MATRICES

#ifdef USE_SPARSE_MATRICES
  #include <Eigen/SparseCore>
  typedef Eigen::SparseMatrix<double> Matrix_t;
  typedef Eigen::Triplet<double> Triplet_t;
#else
  typedef Eigen::MatrixXd Matrix_t;
#endif

class DifferentialEquationException : public std::runtime_error
{
  public:
    DifferentialEquationException() :
      std::runtime_error("DifferentialEquationException") {}
    DifferentialEquationException(const std::string & str) :
      std::runtime_error(str) {}
};



// DiffEqFunctor is a pure virtual class. To derive a class from it:
//   define getDerivatives()
//   if there are any supplementary variables that could be targeted by traces
//     call addTraceTarget() for each one
//   then call ONE of the two initializeDiffEq() functions
//     (depending on if you want default state names or custom ones)
//   then call solve(initialStates)


class DiffEqFunctor
{
  public:
    DiffEqFunctor() : tInitial (0.0), tFinal (Inf), accuracyGoal (0.1),
                    minDeltaT (0.0), stepSizeGain (0.5), algorithmOrder (4.0),
                    maxSolveFailures (10), maxStepFailures (100),
                    verbosity (0) {}
    virtual ~DiffEqFunctor() {}
    
    // initialize differential equation by attaching any traces
    void initializeDiffEq(std::list<Trace> & traces);
    // initialize differential equation by setting the number of states,
    //   and attaching any traces
    void initializeDiffEq(std::list<Trace> & traces, size_t numStates);
    // initialize differential equation by setting the state names, (and
    //   optionally state units), and attaching any traces
    void initializeDiffEq(std::list<Trace> & traces,
                          const std::vector<std::string> & stateNames,
                          const std::vector<std::string> & stateUnits = {});
    // add a state to the DiffEqFunctor (before calling initialize), and
    // return the index to that state
    size_t addStateValue(std::string stateName="", std::string stateUnit="1");
    
    // solve the differential equation, starting with initial states,
    //   and integrating until tFinal
    void solve(const std::vector<double> & initialStates);
    // return the states after solve() completes
    std::vector<double> getFinalStates(void) const;
    // return the initial states
    std::vector<double> getInitialStates(void) const;
    // return the state names
    const std::vector<std::string> & getStateNames(void) const {
      return stateNames;
    }
    // return the number of function evaluations
    size_t getNumFunctionEvaluations(void) const
      { return numDerivatives; }
    // return the number of matrix exponentials
    size_t getNumMatrixExponentials(void) const
      { return numMatrixExponentials; }
    
    // starting time for solve()
    double tInitial;
    // stopping time for solve()
    double tFinal;
    
    // maximum error of any state over entire solve time interval
    double accuracyGoal;
    // minimum deltaT allowed for any step
    double minDeltaT;
    // how closely to adjust step size to max estimated size
    //   in open interval (0, 1)
    double stepSizeGain;
    // what is the effective order of this algorithm?
    double algorithmOrder;
    
    // maximum number of times the solver can fail to converge
    size_t maxSolveFailures;
    // maximum number of times a timestep can fail to converge
    size_t maxStepFailures;
    // set level of output (0 = silent, increase for more output)
    int verbosity;
    
    // traces for perturbing and/or recording differential equation states
    std::list<Trace> traces;
  
  protected:
    // compute d/dt of each state
    virtual void getDerivatives(const std::vector<FDouble> & states,
                                const double & time,
                                std::vector<FDouble> & dStates) const = 0;
    // let traces target a variable for recording/clamping/fitting
    void addTraceTarget(double & target, const std::string & targetName,
                        const std::string & targetUnits);
    
  private:
    // attach traces to differential equation object, so that recording, etc
    //   can be carried out
    void attachTraces(std::list<Trace> & traces);
    // save the initial states in case solver must restart
    void saveInitialStates(const std::vector<double> & initialStates);
    // initialize bookkeeping for solve
    void initializeSolve(void);
    // allocate space for vectors, matrices, etc
    void allocateObjects(void);
    // return true if simulation is not complete, false otherwise
    bool continueSteps(void) const;
    // Repeatedly attempt to solve for the first time step
    inline void solveFirstStep(void);
    // Repeatedly attempt to solve for the next time step
    inline void solveStep(void);
    // update states, traces, bookkeeping
    void updateAfterSuccessfulTimestep(void);
    // set the value of deltaT for the first attempted step
    //  -in event of minor error, set deltaT = 0 to restart ::solve()
    void setFirstDeltaT(void);
    // set the value of deltaT for the next call to timestep()
    //  -in event of minor error, set deltaT = 0 to restart ::solve()
    void setDeltaT(void);
    // multiply deltaT by scale factor chosen based upon past error
    void scaleDeltaT(void);
    // predict states and derivatives of initial step
    void predictFirstStep(void);
    // predict state derivatives of next time step
    void predictStep(void);
    // advance the solution in time
    void timeStep(void);
    // get the largest time step that won't step past a trace action
    double getMaxTraceDeltaT(void) const;
    // predict the value of states at the middle and end of the time step
    void predictFutureStates(void);
    // get derivativesMat and derivativesVec, specifying instantaneous
    // linear approximation to differential equation:
    //    d/dt states = derivativesMat * states + derivativesVec
    void getLinearizedDerivatives(const Eigen::VectorXd & stateVec,
                                  const double & time,
                                  Eigen::VectorXd & dStateVec,
                                  Matrix_t & derivativesMat,
                                  Eigen::VectorXd & derivativesVec);
    // use Magnus expansion to solve differential equation
    void solveDynamicsEquations(void);
    // compute error (difference between predicted and actual states)
    double getPredictionError(void);
    // update the state predictor
    void updatePredictor(void);
    // predict the value of a given state specified by coefficients, at t + dT
    double predictState(const double & dT,
                        const std::vector<double> & coefs) const;
    // compute augmented magnus matrix to solve time endpoint at 2nd-order
    void fillMagnus2ndOrder(const Matrix_t & dMat,
                            const Eigen::VectorXd & dVec, const double & dT);
    // compute augmented magnus matrix to solve time endpoint at 4th-order
    void fillMagnusHalf4thOrder(const Matrix_t & dMat1, const Eigen::VectorXd & dVec1,
                            const Matrix_t & dMat2, const Eigen::VectorXd & dVec2,
                            const Matrix_t & dMat3, const Eigen::VectorXd & dVec3,
                            const double & dT);
    // compute augmented magnus matrix to solve time endpoint at 4th-order
    void fillMagnus4thOrder(const Matrix_t & dMat1, const Eigen::VectorXd & dVec1,
                            const Matrix_t & dMat2, const Eigen::VectorXd & dVec2,
                            const Matrix_t & dMat3, const Eigen::VectorXd & dVec3,
                            const double & dT);
    // compute augmented magnus matrix to solve time endpoint at 6th-order
    void fillMagnusPhony6thOrder(const Matrix_t & dMat1, const Eigen::VectorXd & dVec1,
                                 const Matrix_t & dMat2, const Eigen::VectorXd & dVec2,
                                 const Matrix_t & dMat3, const Eigen::VectorXd & dVec3,
                                 const double & dT);
    
    // estimate error in magnus expansion
    double getMagnusError(void);
    #ifndef USE_SPARSE_MATRICES
    // solve inhomogenous equation with augmented matrix:
    //   [statesFinal, 1] = exp(magnus) * [statesInitial, 1]
    // via Padé approximant
    void expSolvePade(Eigen::VectorXd & statesFinal,
                      const Eigen::VectorXd & statesInitial);
    // helper functions for expSolvePade:
    int scaleMagnus(const double & norm, const double & maxNorm);
    void pade3Solve(Eigen::VectorXd & statesFinal,
                    const Eigen::VectorXd & statesInitial,
                    size_t numSquares = 0);
    void pade5Solve(Eigen::VectorXd & statesFinal,
                    const Eigen::VectorXd & statesInitial,
                    size_t numSquares = 0);
    void pade7Solve(Eigen::VectorXd & statesFinal,
                    const Eigen::VectorXd & statesInitial,
                    size_t numSquares = 0);
    void pade9Solve(Eigen::VectorXd & statesFinal,
                    const Eigen::VectorXd & statesInitial,
                    size_t numSquares = 0);
    void pade13Solve(Eigen::VectorXd & statesFinal,
                     const Eigen::VectorXd & statesInitial,
                     size_t numSquares=0);
    void padeSolveSolution(Eigen::VectorXd & statesFinal,
                           const Eigen::VectorXd & statesInitial);
    void padeSolveSolutionWithSquares(Eigen::VectorXd & statesFinal,
                                      const Eigen::VectorXd & statesInitial,
                                      size_t numSquares=0);
    #endif
    // solve inhomogenous equation with augmented matrix:
    //   [statesFinal, 1] = exp(magnus) * [statesInitial, 1]
    // via Truncated Taylor Series
    //  Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
    // Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
    //  970-989, 2009.
    void expSolveTaylor(Eigen::VectorXd & statesFinal,
                        const Eigen::VectorXd & statesInitial);
    void expSolveTaylorInner(Eigen::VectorXd & statesFinal,
                             const double & stopTol,
                             const size_t & numSteps,
                             const size_t & taylorDegree);
    // Compute optimal number of steps and degree of truncated Taylor series
    std::tuple<size_t, size_t>
    getTruncatedParameters(const double & normMagnus, const double & stopTol);
    // Estimate the p-th root of the 1-norm of magnus^p, for an integer p
    double estimatePNorm(const size_t & p);
    
    // initial values of state variables that evolve according to differential
    // equation
    Eigen::VectorXd startStates;
    
    // states and predicted states at beginning (1), midpoint (2), and end (3)
    // of interval
    Eigen::VectorXd states1;
    Eigen::VectorXd states2;
    Eigen::VectorXd states3;
    Eigen::VectorXd states3Order6;  // states3 computed differently
    Eigen::VectorXd predictedStates;
    
    // matrix and vector linearization of differential equation
    Matrix_t derivativesMatM1;
    Eigen::VectorXd derivativesVecM1;
    Matrix_t derivativesMat0;
    Eigen::VectorXd derivativesVec0;
    Matrix_t derivativesMat1;
    Eigen::VectorXd derivativesVec1;
    Matrix_t derivativesMat1Left;
    Eigen::VectorXd derivativesVec1Left;
    Matrix_t derivativesMat2;
    Eigen::VectorXd derivativesVec2;
    Matrix_t derivativesMat3;
    Eigen::VectorXd derivativesVec3;
    Matrix_t commMat1;
    Eigen::VectorXd commVec1;
    Matrix_t commMat2;
    Eigen::VectorXd commVec2;
    
    // helper matrices and vectors for computing magnus
    Matrix_t meanDerivMat;
    Eigen::VectorXd meanDerivVec;
    Matrix_t deltaDerivMat;
    Eigen::VectorXd deltaDerivVec;
    Matrix_t quadraticDerivMat;
    Eigen::VectorXd quadraticDerivVec;
    #ifdef USE_SPARSE_MATRICES
      Matrix_t magnusMat;  // sparse numStates x numStates upper left section
      Eigen::VectorXd magnusVec;  // dense numStates x 1 upper right section
      double testVectorB; // bottom entry of testVector
    #else
      // augmented magnus matrix (homogenous + inhomogenous parts)
      Matrix_t magnus;
      // helper matrices for expSolvePade
      Matrix_t magnus2;
      Matrix_t magnus4;
      Matrix_t magnus6;
      Matrix_t magnus8;
      Matrix_t evenTerms;
      Matrix_t oddTerms;
      Matrix_t numerator;
      Matrix_t denominator;
      Matrix_t numeratorTimesStates1;
    #endif
    // helper matrices and vectors for expSolveTaylor
    double magnusBR;  // bottom right entry
    Eigen::VectorXd taylorTerm;
    Eigen::VectorXd testVector;
    Matrix_t testMatrix;
    Eigen::VectorXi nonNegVector;
    Eigen::VectorXi oldNonNegVec;
    Eigen::VectorXd magnusPNonNeg;
    
    // d/dt of states
    Eigen::VectorXd dStates;
    // saved value of d/dt of states2
    Eigen::VectorXd dStates2;
    // states vector that supports forward differentiation
    std::vector<FDouble> fStates;
    // d/dt states vector that supports forward differentitation
    std::vector<FDouble> dFStates;
    
    // simulation time
    long double t;
    // maximum remaining time until solve terminates
    long double tRemaining;
    // duration of current time step
    double deltaT;
    // discrepency between Magnus expansions of different orders
    double magnusError;
    // discrepency between predicted and computed states
    double predictionError;
    // specify maximum error per time for a time step
    double errPerT;
    // specify maximum error for this time step (now that deltaT is specified)
    double maxErr;
    // exponent used in setting deltaT given magnus error
    double errorPow;
    // minimum deltaT used by any time step in solution
    double minSolutionDeltaT;
    // maximum deltaT for this step (e.g. to not step past a trace action)
    double maxStepDeltaT;
    // deltaT from previous successful step
    double deltaTOld;
    // duration of the shortest trace
    double shortestTraceTime;
    // pointer to Trace with shortest time
    Trace const* shortestTracePtr;
    // pointer to Trace that terminates solve()
    Trace const* terminalTracePtr;
    // are the derivatives discontinuous at t1?
    bool discontinuous1;
    // do we need to restart due to error?
    bool restartSteps;
    
    // names of the state variables
    std::vector<std::string> stateNames;
    // units of the state variables (default to "none")
    std::vector<std::string> stateUnits;
    // map relating names of Trace targets to their underlying data
    std::map<std::string, TraceTarget> traceTargetMap;
    // number of states
    size_t numStates;
    // number of successful time steps
    size_t numTimeSteps;
    // total number of time steps, successful or no
    size_t totalTimeSteps;
    // number of times the solver failed to converge
    size_t numSolveFailures;
    // number of times this timestep failed to converge
    size_t numStepFailures;
    // greatest numStepFailures accross all steps in solution
    size_t greatestNumStepFailures;
    // number of times the derivatives have been evaluatated
    size_t numDerivatives;
    // number of times a matrix has been exponentiated
    size_t numMatrixExponentials;
};

#endif
