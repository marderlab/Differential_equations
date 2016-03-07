#include "../include/differential_equations.h"
#include "../include/differential_equations.h"

#include <iostream>


using namespace std;
using namespace Eigen;



#define DEBUG_LEVEL 0
// estimate error by using a phony 6th-order magnus expansion
#define USE_PHONY_6TH
// Try to speed up simulation by using less accuracy than machine precision
#define LOW_ACCURACY_MATRIX_EXP
// Use Pade' approximants 
//#define USE_PADE
// choose method of solving matrix equations (if using Pade approximants)
#define USE_LU_DECOMP
//#define USE_QR_DECOMP
// Subtract off main diagonal in expSolveTaylor()
//#define SUB_MEAN_EIGEN

#ifdef USE_PADE
  #ifdef USE_SPARSE_MATRICES
    #error Pade approximants incompatible with sparse matrices
  #endif
#endif

//////////////////////////////// Helper functions /////////////////////////////

#ifdef USE_SPARSE_MATRICES
inline double
induced1Norm(const Matrix_t & mat, const size_t & numStates)
{
  // return the induced 1-norm (maximum over columns of sum of abs value of
  //  each column)
  // Restrict to only return induced 1-norm of the homogenous part
  double maxSum = 0.0;
  for(int col = 0; (size_t)col < numStates; ++col) {
    double sum = 0.0;
    for(Matrix_t::InnerIterator matItr(mat, col); matItr; ++matItr) {
      if((size_t)matItr.row() >= numStates) {
        break;
      }
      sum += abs(matItr.value());
    }
    if(sum > maxSum) {
      maxSum = sum;
    }
  }
  return maxSum;
}
inline double
induced1Norm(const Matrix_t & mat)
{
  // return the induced 1-norm (maximum over columns of sum of abs value of
  //  each column) for the whole matrix mat
  double maxSum = 0.0;
  for(int col = 0; col < mat.outerSize(); ++col) {
    double sum = 0.0;
    for(Matrix_t::InnerIterator matItr(mat, col); matItr; ++matItr) {
      sum += abs(matItr.value());
    }
    if(sum > maxSum) {
      maxSum = sum;
    }
  }
  return maxSum;
}

inline double
inducedInfNorm(const Matrix_t & mat)
{
  // return the induced infinity-norm (maximum over rows of sum of abs value of
  //  each row) for the whole matrix mat
  ArrayXd rowSums;
  rowSums.resize(mat.rows());
  rowSums.fill(0.0);
  for(int col = 0; col < mat.outerSize(); ++col) {
    for(Matrix_t::InnerIterator matItr(mat, col); matItr; ++matItr) {
      rowSums[matItr.row()] += abs(matItr.value());
    }
  }
  return rowSums.maxCoeff();
}
inline double
inducedInfNorm(const Matrix_t & mat, const VectorXd & v)
{
  // return the induced infinity-norm (maximum over rows of sum of abs value of
  //  each row) for the whole matrix mat
  // return the induced infinity-norm (maximum over rows of sum of abs value of
  //  each row) for the whole matrix mat
  ArrayXd rowSums = v.cwiseAbs().array();
  for(int col = 0; col < mat.outerSize(); ++col) {
    for(Matrix_t::InnerIterator matItr(mat, col); matItr; ++matItr) {
      rowSums[matItr.row()] += abs(matItr.value());
    }
  }
  return rowSums.maxCoeff();
}

inline double
minSparseValue(const Matrix_t & mat)
{
  double const* valPtr = mat.valuePtr();
  double const* const stopPtr = valPtr + (size_t)mat.nonZeros();
  double minVal = *valPtr;
  
  for(++valPtr; valPtr != stopPtr; ++valPtr) {
    if(*valPtr < minVal) {
      minVal = *valPtr;
    }
  }
  return minVal;
}
inline bool
isNonNegative(const Matrix_t & mat)
{
  
  double const* const stopPtr = mat.valuePtr() + (size_t)mat.nonZeros();
  for(double const* valPtr = mat.valuePtr(); valPtr != stopPtr; ++valPtr) {
    if(*valPtr < 0) {
      return false;
    }
  }
  return true;
}

#else
inline double
induced1Norm(const Matrix_t & mat, const size_t & numStates)
{
  // return the induced 1-norm (maximum over columns of sum of abs value of
  //  each column)
  // Restrict to only return induced 1-norm of the homogenous part
  
  return mat.topLeftCorner(numStates, numStates).cwiseAbs().colwise().sum().
    maxCoeff();
}
inline double
induced1Norm(const Matrix_t & mat)
{
  // return the induced 1-norm (maximum over columns of sum of abs value of
  //  each column) for the whole matrix mat
  return mat.cwiseAbs().colwise().sum().maxCoeff();
}

inline double
inducedInfNorm(const Matrix_t & mat)
{
  // return the induced infinity-norm (maximum over rows of sum of abs value of
  //  each row) for the whole matrix mat
  return mat.cwiseAbs().rowwise().sum().maxCoeff();
}
inline double
inducedInfNorm(const Matrix_t & mat, const VectorXd & v)
{
  // return the induced infinity-norm (maximum over rows of sum of abs value of
  //  each row) for the whole matrix mat
  return (mat.cwiseAbs().rowwise().sum().array() + v.array()).maxCoeff();
}
#endif

inline double
lInfNorm(const VectorXd & v)
{
  return v.cwiseAbs().maxCoeff();
}


inline ArrayXi
mySign(const VectorXd & vec)
{
  return (vec.array() >= 0).cast<int>();
}


inline double
fastPow(double base, size_t exponent)
{
  // exponentiation by squaring
  if(exponent) {
    while(!(exponent & 1)) {
      // while exponent is even, square base and divide exponent by two
      base *= base;
      exponent >>= 1;
    }
    // exponent is odd, so set current result to base
    double result = base;
    
    // set up loop so that base is only squared if the exponent is non-zero
    
    // divide exponent by two
    exponent >>= 1;
    // loop while exponent is non-zero
    while(exponent) {
      // square base
      base *= base;
      if(exponent & 1)  
        // if exponent is odd, multiply result by current base
        result *= base;
      // divide exponent by two
      exponent >>= 1;
    }
    return result;
  }
  else
    return 1.0;
}



//////////////////////////// Public member functions //////////////////////////



void
DiffEqFunctor::initializeDiffEq(list<Trace> & assignTraces)
{
  // initialize differential equation by attaching any traces
  #if DEBUG_LEVEL >= 1
    cout << "initializeDiffEq(assignTraces, assignNumStates)\n";
    cout.flush();
  #endif

  numStates = stateNames.size();
  
  if(numStates == 0) {
    throw DifferentialEquationException("Must have > 0 states");
  }
  
  // states1 has to be allocated so that traces can look at it
  states1.resize(numStates);
  
  for(size_t n = 0; n < numStates; ++n) {
    addTraceTarget(states1[n], stateNames[n], stateUnits[n]);
  }
  
  attachTraces(assignTraces);

  #if DEBUG_LEVEL >= 1
    cout << "done initialize\n"; cout.flush();
  #endif
}



void
DiffEqFunctor::initializeDiffEq(list<Trace> & assignTraces,
                                size_t assignNumStates)
{
  // initialize differential equation by setting the number of states,
  //   and attaching any traces
  #if DEBUG_LEVEL >= 1
    cout << "initializeDiffEq(assignTraces, assignNumStates)\n";
    cout.flush();
  #endif

  numStates = assignNumStates;
  
  // states1 has to be allocated so that traces can look at it
  states1.resize(numStates);
  
  stateUnits.resize(numStates, "none");
  stateNames.resize(numStates);
  for(size_t n = 0; n < numStates; ++n) {
    stateNames[n] = numToString(n);
    addTraceTarget(states1[n], stateNames[n], stateUnits[n]);
  }
  
  attachTraces(assignTraces);

  #if DEBUG_LEVEL >= 1
    cout << "done initialize\n"; cout.flush();
  #endif
}



void
DiffEqFunctor::initializeDiffEq(list<Trace> & assignTraces,
                                const vector<string> & assignNames,
                                const vector<string> & assignUnits)
{
  // initialize differential equation by setting the number of states,
  //   and attaching any traces
  #if DEBUG_LEVEL >= 1
    cout << "initializeDiffEq(assignTraces, assignNames, assignUnits)\n";
    cout.flush();
  #endif

  stateNames = assignNames;
  stateUnits = assignUnits;
  numStates = stateNames.size();
  
  // states1 has to be allocated so that traces can look at it
  states1.resize(numStates);
  
  if(assignUnits.size() != stateNames.size()) {
    stateUnits.resize(numStates, "none");
  }
  else {
    stateUnits = assignUnits;
  }
  
  for(size_t n = 0; n < numStates; ++n) {
    addTraceTarget(states1[n], stateNames[n], stateUnits[n]);
  }
  
  attachTraces(assignTraces);

  #if DEBUG_LEVEL >= 1
    cout << "done initialize\n"; cout.flush();
  #endif
}



size_t
DiffEqFunctor::addStateValue(string stateName, string stateUnit)
{
  // add a state to the DiffEqFunctor (before calling initialize), and
  // return the index to that state
  size_t stateInd = stateNames.size();

  if(stateName.empty()) {
    stateNames.push_back(numToString(stateInd));
  }
  else {
    stateNames.push_back(std::move(stateName));
  }

  stateUnits.push_back(std::move(stateUnit));
  return stateInd;
}



void
DiffEqFunctor::solve(const vector<double> & initialStates)
{
  // solve the differential equation
  
  #if DEBUG_LEVEL >= 1
    cout << "solve()\n"; cout.flush();
  #endif

  // ensure initialStates is not empty
  if(initialStates.empty()) {
    throw DifferentialEquationException("initialStates is empty");
  }
  // save the initial states in case solver must restart
  saveInitialStates(initialStates);
  
  numSolveFailures = 0;
  do {
    initializeSolve();
    
    if(verbosity > 0) {
      cout << "Starting time steps.\n"; cout.flush();
    }

    // solve for the first time step
    solveFirstStep();
    #if DEBUG_LEVEL >= 1
      cout << "t = " << t << ", dT = " << deltaT << '\n'; cout.flush();
    #endif
    
    while(continueSteps()) {
      #if DEBUG_LEVEL >= 1
        cout << "t = " << t << ", dT = " << deltaT << '\n'; cout.flush();
      #endif
      // loop until solve is successful
      // Repeatedly attempt to solve for the next time step
      solveStep();
    }
  } while(restartSteps && numSolveFailures <= maxSolveFailures);
  
  if(numSolveFailures > maxSolveFailures) {
    throw DifferentialEquationException(
      "Failed to converge to requested accuracy at t = " + numToString(t));
  }
  
  if(verbosity > 0) {
    cout << numTimeSteps << " steps, minDT = "
         << minSolutionDeltaT << ", max failures = " << greatestNumStepFailures
         << ", total step calls = " << totalTimeSteps << '\n';
    cout << numDerivatives << " function evaluations, "
         << numMatrixExponentials << " matrix exponentials\n";
    cout << "Finished time steps.\n"; cout.flush();
    cout.flush();
  }
  #if DEBUG_LEVEL >= 1
    cout << "finished solve()\n"; cout.flush();
  #endif
}



vector<double>
DiffEqFunctor::getFinalStates(void) const
{
  // yield the states after solve() completes
  
  return vector<double> (states3.data(), states3.data() + numStates);
}



vector<double>
DiffEqFunctor::getInitialStates(void) const
{
  // yield the initial states
  
  return vector<double> (startStates.data(), startStates.data() + numStates);
}



////////////////////////// Protected member functions /////////////////////////

void
DiffEqFunctor::addTraceTarget(double & target, const string & targetName,
                              const string & targetUnits)
{
  // let traces target a variable for recording/clamping/fitting
  traceTargetMap[targetName] = TraceTarget(target, targetUnits);
}



/////////////////////////// Private member functions //////////////////////////

void
DiffEqFunctor::attachTraces(list<Trace> & assignTraces)
{
  // attach traces to differential equation object, so that recording, etc
  //   can be carried out
  
  #if DEBUG_LEVEL >= 1
    cout << "\tattachTraces()\n"; cout.flush();
  #endif
  
  traces = assignTraces;
  shortestTraceTime = Inf;
  shortestTracePtr = NULL;
  for(auto & trace : traces) {
    const double traceTime = trace.deltaT * (double)(trace.trace.size() - 1);
    if(traceTime < shortestTraceTime) {
      shortestTraceTime = traceTime;
      shortestTracePtr = &trace;
    }
    
    size_t numMatches = 0;
    for(auto & target : traceTargetMap) {
      if(target.first == trace.targetName) {
        // found correct state name, so set trace target
        #if DEBUG_LEVEL >= 1
          cout << "\t\tsetting target " << trace.targetName; cout.flush();
        #endif
        trace.setTarget(target.second);
        ++numMatches;
        #if DEBUG_LEVEL >= 1
          cout << " ok\n"; cout.flush();
        #endif
      }
    }
    if(numMatches != 1) {
      cout << "Trace targets:\n";
      for(auto & target : traceTargetMap) {
        cout << "\t" << target.first << '\n';
      }
      cout.flush();
      if(numMatches == 0) {
        throw DifferentialEquationException("Invalid trace target: " +
                                            trace.targetName);
      }
      else {
        throw DifferentialEquationException(
          "Multiple matches to trace target: " + trace.targetName);
      }
    }
  }
  #if DEBUG_LEVEL >= 1
    cout << "\tDone attaching " << traces.size() << " traces\n"; cout.flush();
  #endif
}



void
DiffEqFunctor::saveInitialStates(const vector<double> & initialStates)
{
  // save the initial states in case solver must restart
  
  numStates = initialStates.size();
  startStates = Map<const VectorXd> (&initialStates[0], numStates, 1);
}



void
DiffEqFunctor::initializeSolve(void)
{
  // initialize bookkeeping for solve
  #if DEBUG_LEVEL >= 1
    cout << "\tinitializeSolve()\n"; cout.flush();
  #endif

  numTimeSteps = 0;
  totalTimeSteps = 0;
  greatestNumStepFailures = 0;
  numDerivatives = 0;
  numMatrixExponentials = 0;
  minSolutionDeltaT = Inf;
  tRemaining = tFinal - tInitial;
  discontinuous1 = false;
  restartSteps = false;

  if(numSolveFailures == 0) {
    // allocate space for vectors, matrices, etc
    if(fStates.size() != numStates) {
      allocateObjects();
    }
    // set errPerT to default value
    errPerT = accuracyGoal / sqrt((double)tRemaining) / stepSizeGain;
    // set error pow
    errorPow = 1.0 / (algorithmOrder - 0.5);
    // set terminalTracePtr if the shortest trace is responsible for stopping
    //   solve:
    terminalTracePtr = shortestTraceTime <= tRemaining ? shortestTracePtr :
                                                         NULL;
  }
  else {
    if(verbosity > 0) {
      cout << "Failed to converge to required accuracy at t = " << t
           << ". Restarting time steps.\n"; cout.flush();
    }
    // make stepSizeGain more conservative (hopefully this will avoid future
    // failure)
    stepSizeGain *= 0.1;
  }

  t = tInitial;  
  states1 = startStates;
  
  // execute traces before first step
  bool dummy;
  for(auto & trace : traces) {
    trace.initialize();
    if(trace.traceAction != Trace::Record) {
      trace.doTraceAction(dummy);
    }
  }
  
  getLinearizedDerivatives(states1, tInitial, dStates,
                           derivativesMat1, derivativesVec1);

  for(auto & trace : traces) {
    if(trace.traceAction == Trace::Record) {
      trace.doTraceAction(dummy);
    }
  }

  #if DEBUG_LEVEL >= 1
    cout << "\tfinished initializeSolve()\n"; cout.flush();
  #endif
}



void
DiffEqFunctor::allocateObjects(void)
{
  // allocate space for vectors, matrices, etc that will be accessed by index
  //  before being assigned (i.e. things that *must* be allocated)
  
  dStates.resize(numStates);
  dStates2.resize(numStates);
  
  derivativesMatM1.resize(numStates, numStates);
  derivativesVecM1.resize(numStates);
  derivativesMat0.resize(numStates, numStates);
  derivativesVec0.resize(numStates);
  derivativesMat1Left.resize(numStates, numStates);
  derivativesVec1Left.resize(numStates);  
  derivativesMat1.resize(numStates, numStates);
  derivativesVec1.resize(numStates);
  derivativesMat2.resize(numStates, numStates);
  derivativesVec2.resize(numStates);
  derivativesMat3.resize(numStates, numStates);
  derivativesVec3.resize(numStates);

  
  #ifdef USE_SPARSE_MATRICES
    magnusMat.resize(numStates, numStates);
    magnusVec.resize(numStates);
  #else
    derivativesMatM1.fill(0.0);
    derivativesMat0.fill(0.0);
    derivativesMat1Left.fill(0.0);
    derivativesMat1.fill(0.0);
    derivativesMat2.fill(0.0);
    derivativesMat3.fill(0.0);
    magnus.resize(numStates + 1, numStates + 1);
    magnus.fill(0.0);
  #endif
  
  dStates.resize(numStates);
  dStates2.resize(numStates);
  
  // prepare automatic differentiation
  fStates.resize(numStates);
  dFStates.resize(numStates);
  size_t n = 0;
  for(auto & fState : fStates) {
    fState.setDerivativeIndex(n);
    ++n;
  }
  
  // these could be accessed if a magnus matrix has non-finite elements
  // so they must be resized:
  states2.resize(numStates);
  states3.resize(numStates);
  predictedStates.resize(numStates);
}



bool
DiffEqFunctor::continueSteps(void) const
{
  // return true if simulation is not complete, false otherwise
  if(terminalTracePtr == NULL) {
    // simulation isn't terminated by a trace, solve while time remaining
    return tRemaining > 0 && !restartSteps;
  }
  else {
    // simulation terminates when the shortest traces reaches its end
    return !(terminalTracePtr->reachedEnd() || restartSteps);
  }
  
}



void
DiffEqFunctor::solveFirstStep(void)
{
  // Repeatedly attempt to solve for the first time step
  #if DEBUG_LEVEL >= 1
    cout << "\tsolveFirstStep()\n"; cout.flush();
  #endif

  numStepFailures = 0;
  
  // choose deltaT and make first attempt at time step
  setFirstDeltaT();
  if(!(deltaT > 0.0)) {
    ++numSolveFailures;
    restartSteps = true;
    return;
  }
  
  #if DEBUG_LEVEL >= 1
    cout << "\t\tPredicting first step\n"; cout.flush();
  #endif
  // predict states of first step
  predictFirstStep();
  
  // loop until solution is accurate enough or maxStepFailures is reached
  while(numStepFailures < maxStepFailures) {
    if(magnusError < maxErr) {
      // refine solution
      #if DEBUG_LEVEL >= 1
        cout << "\t\tRefining step with dT = " << deltaT << '\n';
        cout.flush();
      #endif
      timeStep();

      if(predictionError < maxErr) {
        // step is solved to sufficient accuracy
        break;
      }
      else {
        ++numStepFailures;
      }
    }
    else {
      // prediction failed or got NaN results from time step
      ++numStepFailures;

      // error is too large, so shrink deltaT
      setDeltaT();
      if(!(deltaT > 0.0)) {
        ++numSolveFailures;
        restartSteps = true;
        return;
      }

      #if DEBUG_LEVEL >= 1
        cout << "\t\tPredicting first step\n"; cout.flush();
      #endif
      predictFirstStep();
    }
  }

  if(numStepFailures >= maxStepFailures) {
    // was never able to make a reasonable estimate
    ++numSolveFailures;
    restartSteps = true;
    return;
  }
  #if DEBUG_LEVEL >= 1
    cout << "\t\tFinished refining\n";
    cout << "\t\tUpdating after successful step\n";
    cout.flush();
  #endif
  
  // update states, traces, bookkeeping
  updateAfterSuccessfulTimestep();
  #if DEBUG_LEVEL >= 1
    cout << "\tfinished solveFirstStep()\n"; cout.flush();
  #endif

}



void
DiffEqFunctor::solveStep(void)
{
  // Repeatedly attempt to solve for the next time step
  numStepFailures = 0;

  // choose deltaT and make first attempt at time step
  setDeltaT();
  if(!(deltaT > 0.0)) {
    ++numSolveFailures;
    restartSteps = true;
    return;
  }
  
  // predict states of next step
  predictStep();
  
  // loop until solution is accurate enough or maxStepFailures is reached
  while(numStepFailures < maxStepFailures) {
    if(magnusError < maxErr) {
      // refine solution
      timeStep();

      if(predictionError < maxErr) {
        // step is solved to sufficient accuracy
        break;
      }
      else {
        ++numStepFailures;
      }
    }
    else {
      // prediction failed or got NaN results from time step
      ++numStepFailures;

      // error is too large, so shrink deltaT
      setDeltaT();
      if(!(deltaT > 0.0)) {
        ++numSolveFailures;
        restartSteps = true;
        return;
      }
      predictStep();
    }
  }
  
  if(numStepFailures >= maxStepFailures) {
    // was never able to make a reasonable estimate
    ++numSolveFailures;
    restartSteps = true;
    return;
  }
  
  // update states, traces, bookkeeping
  updateAfterSuccessfulTimestep();
}



void
DiffEqFunctor::updateAfterSuccessfulTimestep(void)
{
  // update states, traces, bookkeeping

  // update states
  states1 = states3;
  
  // update oldest derivatives
  derivativesMatM1.swap(derivativesMat1);
  derivativesVecM1.swap(derivativesVec1);
  derivativesMat0.swap(derivativesMat2);
  derivativesVec0.swap(derivativesVec2);
  
    // advance model time
  t += (long double)deltaT;
  
  // reduce remaining simulation time
  if(deltaT >= (double)tRemaining) {
    tRemaining = 0;
  }
  else {
    tRemaining -= (long double)deltaT;
  }
  
  deltaTOld = deltaT;

  // advance the traces and carry out their action for the next time step
  discontinuous1 = false;
  for(auto & trace : traces) {
    trace.advance(deltaT);
    if(trace.traceAction != Trace::Record) {
      trace.doTraceAction(discontinuous1);
    }
  }
  
  // update newest derivatives
  if(discontinuous1) {
    derivativesMat1Left.swap(derivativesMat3);
    derivativesVec1Left.swap(derivativesVec3);
    // need to recompute derivatives due to action of trace
    getLinearizedDerivatives(states1, (double)t,
                             dStates, derivativesMat1, derivativesVec1);
  }
  else {
    derivativesMat1.swap(derivativesMat3);
    derivativesVec1.swap(derivativesVec3);
  }
  
  bool dummy;
  for(auto & trace : traces) {
    if(trace.traceAction == Trace::Record) {
      trace.doTraceAction(dummy);
    }
  }

  
  // bookkeeping
  ++numTimeSteps;
  if(numStepFailures > greatestNumStepFailures) {
    greatestNumStepFailures = numStepFailures;
  }
  if(deltaT < minSolutionDeltaT) {
    minSolutionDeltaT = deltaT;
  }
}



void
DiffEqFunctor::setFirstDeltaT(void)
{
  // set the value of deltaT for the first attempted time step
  
  double maxDeltaT;
  
  maxDeltaT = getMaxTraceDeltaT();
  if(terminalTracePtr == NULL) {
    maxDeltaT = std::min((double)tRemaining, maxDeltaT);
  }
  
  // set conservative guess for first value of deltaT
  const double errorRatio = 0.01;
  const double maxDeriv = dStates.cwiseAbs().maxCoeff();
  deltaT = pow(errPerT / pow(errorRatio * maxDeriv, algorithmOrder),
               errorPow);
  
  // adjust to not overstep maxDeltaT
  if(deltaT > maxDeltaT) {
    // desired deltaT is too large, take maximum allowed step
    deltaT = maxDeltaT;
  }
  else {
    if(deltaT < minDeltaT) {
      deltaT = minDeltaT;
    }

    if(!traces.empty()) {
      // desired deltaT is small, divide time step into integer number of smaller
      // steps, consistent with desired deltaT
      const double maxTraceDeltaT = getMaxTraceDeltaT();
      const double numDivs = ceil(maxTraceDeltaT / deltaT);
      deltaT = 1.01 * maxTraceDeltaT / numDivs;
    }
  }

  // compute maximum allowed error
  maxErr = errPerT * sqrt(deltaT);
}



void
DiffEqFunctor::setDeltaT(void)
{
  // set the value of deltaT for the next call to timestep()
  
  double maxDeltaT;
  
  if(numStepFailures == 0) {
    // first attempt at this step, try increasing step size
    maxDeltaT = getMaxTraceDeltaT();
    if(terminalTracePtr == NULL) {
      maxDeltaT = std::min((double)tRemaining, maxDeltaT);
    }
    // scale deltaT
    scaleDeltaT();
  }
  else {
    // failed previous step, try choosing next step size by error scaling
    
    // don't expand to longer than previous step
    maxDeltaT = deltaT;
    // scale deltaT
    scaleDeltaT();
  }
  
  // adjust to not overstep maxDeltaT
  if(deltaT > maxDeltaT) {
    // desired deltaT is too large, take maximum allowed step
    deltaT = maxDeltaT;
  }
  else {
    if(deltaT < minDeltaT) {
      deltaT = minDeltaT;
    }
    
    if(!traces.empty()) {
      // desired deltaT is small, divide time step into integer number of smaller
      // steps, consistent with desired deltaT
      const double maxTraceDeltaT = getMaxTraceDeltaT();
      const double numDivs = ceil(maxTraceDeltaT / deltaT);
      deltaT = 1.01 * maxTraceDeltaT / numDivs;
    }
  }

  // compute maximum allowed error
  maxErr = errPerT * sqrt(deltaT);
}



inline void
DiffEqFunctor::scaleDeltaT(void)
{
  // multiply deltaT by scale factor chosen based upon past error
  const double scaleFact = pow(stepSizeGain * maxErr / magnusError, errorPow);
  if(scaleFact > 0.01) {
    deltaT *= scaleFact;
  }
  else {
    deltaT *= 0.01;
  }
}



void
DiffEqFunctor::predictFirstStep(void)
{
  // predict states and derivatives of initial step
  
  // predict states2
  //   -fill magnus
  const double halfDT = 0.5 * deltaT;
  fillMagnus2ndOrder(derivativesMat1, derivativesVec1, halfDT);

  //   -solve states2 = exp(magnus) * states1
  #ifdef USE_PADE
    expSolvePade(states2, states1);
  #else
    expSolveTaylor(states2, states1);
  #endif
  
  // get derivativesMat/derivativesVec at estimated states2
  getLinearizedDerivatives(states2, (double)t + halfDT,
                           dStates, derivativesMat2, derivativesVec2);

  // predict derivatives of states3
  derivativesMat3 = 2 * derivativesMat2 - derivativesMat1;
  derivativesVec3 = 2 * derivativesVec2 - derivativesVec1;
  
  // compute states3
  //   -fill magnus
  fillMagnus4thOrder(derivativesMat1, derivativesVec1,
                     derivativesMat2, derivativesVec2,
                     derivativesMat3, derivativesVec3,
                     deltaT);

  //   -solve states3 = exp(magnus) * states2
  #ifdef USE_PADE
    expSolvePade(states3, states1);
  #else
    expSolveTaylor(states3, states1);
  #endif

  // estimate error in magnus expansion
  magnusError = getMagnusError();
}



void
DiffEqFunctor::predictStep(void)
{
  // predict state derivatives of next time step

  // predict derivatives of states2, states3
  const double coef1 = deltaT / deltaTOld;
  const double coef2 = coef1 * coef1;
  
  if(discontinuous1) {
    derivativesMat2 = derivativesMat1
      + (0.5 * coef1) *
         (3.0 * derivativesMat1Left + derivativesMatM1 - 4.0 * derivativesMat0)
      + (0.5 * coef2) *
         (derivativesMat1Left + derivativesMatM1 - 2.0 * derivativesMat0);
    derivativesVec2 = derivativesVec1
      + (0.5 * coef1) *
         (3.0 * derivativesVec1Left + derivativesVecM1 - 4.0 * derivativesVec0)
      + (0.5 * coef2) *
         (derivativesVec1Left + derivativesVecM1 - 2.0 * derivativesVec0);
    
    derivativesMat3 = derivativesMat1
      + (2 * coef1) *
         (3.0 * derivativesMat1Left + derivativesMatM1 - 4.0 * derivativesMat0)
      + coef2 *
         (derivativesMat1Left + derivativesMatM1 - 2.0 * derivativesMat0);
    derivativesVec3 = derivativesVec1
      + (2 * coef1) *
         (3.0 * derivativesVec1Left + derivativesVecM1 - 4.0 * derivativesVec0)
      + coef2 *
         (derivativesVec1Left + derivativesVecM1 - 2.0 * derivativesVec0);
  }
  else {
    derivativesMat2 = derivativesMat1
      + (0.5 * coef1) *
         (3.0 * derivativesMat1 + derivativesMatM1 - 4.0 * derivativesMat0)
      + (0.5 * coef2) *
         (derivativesMat1 + derivativesMatM1 - 2.0 * derivativesMat0);
    derivativesVec2 = derivativesVec1
      + (0.5 * coef1) *
         (3.0 * derivativesVec1 + derivativesVecM1 - 4.0 * derivativesVec0)
      + (0.5 * coef2) *
         (derivativesVec1 + derivativesVecM1 - 2.0 * derivativesVec0);
    
    derivativesMat3 = derivativesMat1
      + (2 * coef1) *
         (3.0 * derivativesMat1 + derivativesMatM1 - 4.0 * derivativesMat0)
      + coef2 *
         (derivativesMat1 + derivativesMatM1 - 2.0 * derivativesMat0);
    derivativesVec3 = derivativesVec1
      + (2 * coef1) *
         (3.0 * derivativesVec1 + derivativesVecM1 - 4.0 * derivativesVec0)
      + coef2 *
         (derivativesVec1 + derivativesVecM1 - 2.0 * derivativesVec0);
  }
  
  // compute states3
  //   -fill magnus
  fillMagnus4thOrder(derivativesMat1, derivativesVec1,
                     derivativesMat2, derivativesVec2,
                     derivativesMat3, derivativesVec3,
                     deltaT);
  //   -solve states3 = exp(magnus) * states1
  #ifdef USE_PADE
    expSolvePade(states3, states1);
  #else
    expSolveTaylor(states3, states1);
  #endif

  // get the difference between phony 6th-0order and true 4th-order states
  magnusError = getMagnusError();
}



void
DiffEqFunctor::timeStep(void)
{
  // advance the solution in time
  // get derivativesMat/derivativesVec at estimated states3
  getLinearizedDerivatives(states3, (double)t + deltaT,
                           dStates, derivativesMat3, derivativesVec3);
  
  // compute states2
  //   -fill magnus
  fillMagnusHalf4thOrder(derivativesMat1, derivativesVec1,
                         derivativesMat2, derivativesVec2,
                         derivativesMat3, derivativesVec3,
                         deltaT);
  //   -solve states2 = exp(magnus) * states1
  #ifdef USE_PADE
    expSolvePade(states2, states1);
  #else
    expSolveTaylor(states2, states1);
  #endif
    
  // get derivativesMat/derivativesVec at estimated states2
  getLinearizedDerivatives(states2, (double)t + 0.5 * deltaT,
                           dStates2, derivativesMat2, derivativesVec2);
  // save older states3 as the prediction from previous step
  predictedStates = states3;
  // compute states3
  //   -fill magnus
  fillMagnus4thOrder(derivativesMat1, derivativesVec1,
                     derivativesMat2, derivativesVec2,
                     derivativesMat3, derivativesVec3,
                     deltaT);
  
  //   -solve states3 = exp(magnus) * states1
  #ifdef USE_PADE
    expSolvePade(states3, states1);
  #else
    expSolveTaylor(states3, states1);
  #endif
  
  // compute the prediction error
  predictionError = getPredictionError();
  
  // increment number of attempted time steps
  ++totalTimeSteps;
}



inline double
DiffEqFunctor::getMaxTraceDeltaT(void) const
{
  // get the largest time step that won't step past a trace action
  
  double maxDeltaT = Inf;
  
  for(const auto & trace : traces) {
    const double nextTraceT = trace.getNextT();
    if(nextTraceT < maxDeltaT)
      maxDeltaT = nextTraceT;
  }
  return maxDeltaT;
}


void
DiffEqFunctor::getLinearizedDerivatives(const VectorXd & stateVec,
                                        const double & time,
                                        VectorXd & dStateVec,
                                        Matrix_t & derivativesMat,
                                        VectorXd & derivativesVec)
{
  // get derivativesMat and derivativesVec, specifying instantaneous
  // linear approximation to differential equation:
  //    d/dt states = derivativesMat * states + derivativesVec
  
  // fill fStates
  for(size_t n = 0; n < numStates; ++n) {
    fStates[n].value = stateVec[n];
  }
  
  // compute the derivatives of each state
  getDerivatives(fStates, time, dFStates);
  ++numDerivatives;
  
  #ifdef USE_SPARSE_MATRICES
    size_t n = 0;
    // for now, always set from triplets
    if(true || derivativesMat.nonZeros() == 0) {
      std::vector<Triplet_t> derivTriplets;
      for(const auto & dFState : dFStates) {
        dStateVec[n] = dFState.value;
        for(const auto & derivObj : dFState.derivatives) {
          derivTriplets.push_back(
            Triplet_t(n, derivObj.first, derivObj.second)
          );
        }
        ++n;
      }
      derivativesMat.setFromTriplets(derivTriplets.begin(),
                                     derivTriplets.end());
    }
    else {
      // maybe improve speed by iterating over sparse matrix?
      //   -might require switching Matrix_t from Column-major order to
      //    Row-major order
      for(const auto & dFState : dFStates) {
        dStateVec[n] = dFState.value;
        for(const auto & derivObj : dFState.derivatives) {
          derivativesMat.coeffRef(n,derivObj.first) = derivObj.second;
        }
        ++n;
      }
    }
  #else
    size_t n = 0;
    for(const auto & dFState : dFStates) {
      dStateVec[n] = dFState.value;
      for(const auto & derivObj : dFState.derivatives) {
        derivativesMat(n,derivObj.first) = derivObj.second;
      }
      ++n;
    }
  #endif  
  // compute derivativesVec
  derivativesVec.noalias() = dStateVec - derivativesMat * stateVec;
}



inline double
DiffEqFunctor::getPredictionError(void)
{
  // compute error (difference between predicted and actual states)
  if(states3.allFinite() && predictedStates.allFinite()) {
    return (predictedStates - states3).cwiseAbs().maxCoeff();
  }
  else {
    // Serious corruption of computed/predicted states
    magnusError = Inf;
    return Inf;
  }
}



inline void
DiffEqFunctor::fillMagnus2ndOrder(const Matrix_t & dMat, const VectorXd & dVec,
                                  const double & dT)
{
  // compute augmented magnus matrix to solve time endpoint at 2nd-order
  #ifdef USE_SPARSE_MATRICES
    magnusMat = dT * dMat;
    magnusVec = dT * dVec;
  #else
    magnus.topLeftCorner(numStates, numStates) = dT * dMat;
    magnus.rightCols<1>().head(numStates) = dT * dVec;
  #endif
}



inline void
DiffEqFunctor::fillMagnusHalf4thOrder(const Matrix_t & dMat1,
                                      const VectorXd & dVec1,
                                      const Matrix_t & dMat2,
                                      const VectorXd & dVec2,
                                      const Matrix_t & dMat3,
                                      const VectorXd & dVec3,
                                      const double & dT)
{
  // compute augmented magnus matrix to solve time midpoint at 4th-order
  
  // compute helper matrices/vectors
  //   (note, "mean" derivatives coefficients sum to 1/2 deltaT, because this
  //    expansion is intended to integrate over half the interval)
  const double coef = dT / 24.0;
  meanDerivMat = 5 * coef * dMat1
               + 8 * coef * dMat2
               - coef     * dMat3;
  meanDerivVec = 5 * coef * dVec1
               + 8 * coef * dVec2
               - coef     * dVec3;
  deltaDerivMat = coef * (derivativesMat2 - derivativesMat1);
  deltaDerivVec = coef * (derivativesVec2 - derivativesVec1);
  
  // fill magnus
  #ifdef USE_SPARSE_MATRICES
    magnusMat = meanDerivMat
      + deltaDerivMat * meanDerivMat - meanDerivMat * deltaDerivMat;
    magnusVec = meanDerivVec
      + deltaDerivMat * meanDerivVec - meanDerivMat * deltaDerivVec;
  #else
    magnus.topLeftCorner(numStates, numStates) = meanDerivMat;
    magnus.topLeftCorner(numStates, numStates).noalias() +=
      deltaDerivMat * meanDerivMat - meanDerivMat * deltaDerivMat;
    
    magnus.rightCols<1>().head(numStates) = meanDerivVec;
    magnus.rightCols<1>().head(numStates).noalias() +=
      deltaDerivMat * meanDerivVec - meanDerivMat * deltaDerivVec;
  #endif
}



inline void
DiffEqFunctor::fillMagnus4thOrder(const Matrix_t & dMat1,
                                  const VectorXd & dVec1,
                                  const Matrix_t & dMat2,
                                  const VectorXd & dVec2,
                                  const Matrix_t & dMat3,
                                  const VectorXd & dVec3,
                                  const double & dT)
{
  // compute augmented magnus matrix to solve time endpoint at 4th-order
  
  // compute helper matrices/vectors
  const double coef = dT / 6.0;
  meanDerivMat = coef       * dMat1
               + (4 * coef) * dMat2
               + coef       * dMat3;
  meanDerivVec = coef       * dVec1
               + (4 * coef) * dVec2
               + coef       * dVec3;
  deltaDerivMat = (0.5 * coef) * (dMat3 - dMat1);
  deltaDerivVec = (0.5 * coef) * (dVec3 - dVec1);
    
  // fill magnus
  #ifdef USE_SPARSE_MATRICES
    magnusMat = meanDerivMat
      + deltaDerivMat * meanDerivMat - meanDerivMat * deltaDerivMat;
    
    magnusVec = meanDerivVec
      + deltaDerivMat * meanDerivVec - meanDerivMat * deltaDerivVec;
  #else
    magnus.topLeftCorner(numStates, numStates) = meanDerivMat;
    magnus.topLeftCorner(numStates, numStates).noalias() +=
      deltaDerivMat * meanDerivMat - meanDerivMat * deltaDerivMat;
    
    magnus.rightCols<1>().head(numStates) = meanDerivVec;
    magnus.rightCols<1>().head(numStates).noalias() +=
      deltaDerivMat * meanDerivVec - meanDerivMat * deltaDerivVec;
  #endif
}



inline void
DiffEqFunctor::fillMagnusPhony6thOrder(const Matrix_t & dMat1,
                                       const VectorXd & dVec1,
                                       const Matrix_t & dMat2,
                                       const VectorXd & dVec2,
                                       const Matrix_t & dMat3,
                                       const VectorXd & dVec3,
                                       const double & dT)
{
  // compute augmented magnus matrix to solve time endpoint at 6th-order
  
  // compute helper matrices/vectors
  // note: assume mean and delta matrices/vectors are computed by previous call
  //       to fillMagnus4thOrder
  
  quadraticDerivMat = (2 * dT) * (dMat1 + dMat3) - (4 * dT) * dMat2;
  quadraticDerivVec = (2 * dT) * (dVec1 + dVec3) - (4 * dT) * dVec2;
  
  // compute commutators. note in the following:
  //   alpha1 = dT * d/dt states2
  //   alpha2 = dT * (d/dt states3 - d/dt states1)
  //   alpha3 = 2 * dT * (d/dt states1 + d/dt states3)
  //          - 4 * dT * d/dt states2
  
  //   comm1 = [alpha1, alpha2]
  #ifdef USE_SPARSE_MATRICES
    commMat1 = (12.0 * dT) *
      (dMat2 * deltaDerivMat - deltaDerivMat * dMat2);
    commVec1 = (12.0 * dT) *
      (dMat2 * deltaDerivVec - deltaDerivMat * dVec2);
  #else
    commMat1.noalias() = (12.0 * dT) *
      (dMat2 * deltaDerivMat - deltaDerivMat * dMat2);
    commVec1.noalias() = (12.0 * dT) *
      (dMat2 * deltaDerivVec - deltaDerivMat * dVec2);
  #endif
  
  // comm2 = -1/60 * [alpha1, comm1 + 2 * alpha3]
  //   first add comm1 += 2 * alpha3
  commMat1 += 2 * quadraticDerivMat;
  commVec1 += 2 * quadraticDerivVec;
  #ifdef USE_SPARSE_MATRICES
    commMat2 = (-dT/60.0) * (dMat2 * commMat1 - commMat1 * dMat2);
    commVec2 = (-dT/60.0) * (dMat2 * commVec1 - commMat1 * dVec2);
  #else
    commMat2.noalias() = (-dT/60.0) * (dMat2 * commMat1 - commMat1 * dMat2);
    commVec2.noalias() = (-dT/60.0) * (dMat2 * commVec1 - commMat1 * dVec2);
  #endif
  
  // comm3 = 1/240 [comm1 + 2/3 alpha3 - 20 meanDeriv, alpha3 + comm2]
  // first add comm2 += alpha3
  commMat2 += quadraticDerivMat;
  commVec2 += quadraticDerivVec;
  // next add comm1 += -4/3 alpha3 -20 meanDeriv
  //   (remember comm1 had 2 alpha3 added to it to compute comm2)
  commMat1 += (-4.0/3.0) * quadraticDerivMat - 20 * meanDerivMat;
  commVec1 += (-4.0/3.0) * quadraticDerivVec - 20 * meanDerivVec;
  
  // compute comm3 directly into magnus
  #ifdef USE_SPARSE_MATRICES
    magnusMat = meanDerivMat + (1.0/240.0) *
      (commMat1 * commMat2 - commMat2 * commMat1);
    
    magnusVec = meanDerivVec + (1.0/240.0) *
      (commMat1 * commVec2 - commMat2 * commVec1);
  #else
    magnus.topLeftCorner(numStates, numStates) = meanDerivMat;
    magnus.topLeftCorner(numStates, numStates).noalias() += (1.0/240.0) *
      (commMat1 * commMat2 - commMat2 * commMat1);
    
    magnus.rightCols<1>().head(numStates) = meanDerivVec;
    magnus.rightCols<1>().head(numStates).noalias() += (1.0/240.0) *
      (commMat1 * commVec2 - commMat2 * commVec1);
  #endif
}



inline double
DiffEqFunctor::getMagnusError(void)
{
  // compute phony higher-order magnus to estimate magnus expansion error
  
  #ifdef USE_PHONY_6TH
    fillMagnusPhony6thOrder(derivativesMat1, derivativesVec1,
                            derivativesMat2, derivativesVec2,
                            derivativesMat3, derivativesVec3,
                            deltaT);
  #else
    // this would be the 2nd-order version of this test:
    #ifdef USE_SPARSE_MATRICES
      magnusMat = meanDerivMat;
      magnusVec = meanDerivVec;
    #else
      magnus.topLeftCorner(numStates, numStates) = meanDerivMat;
      magnus.rightCols<1>().head(numStates) = meanDerivVec;
    #endif
  #endif  
  
  #ifdef USE_PADE
    expSolvePade(predictedStates, states1);
  #else
    expSolveTaylor(predictedStates, states1);
  #endif
  
  // get the difference between phony 6th-order and true 4th-order states
  return getPredictionError();
}


#ifndef USE_SPARSE_MATRICES
void
DiffEqFunctor::expSolvePade(VectorXd & statesFinal,
                            const Eigen::VectorXd & statesInitial)
{
  // solve inhomogenous equation with augmented matrix:
  //   [statesFinal, 1] = exp(magnus) * [statesInitial, 1]
  // via Padé approximant
  // gist:  approximate exp(magnus) ~ denominator^-1 * numerator
  //    then denominator * statesFinal = numerator * statesInitial
  //    use LU decomposition to solve this equation for statesFinal
  // compute norm of mat
  const double norm = induced1Norm(magnus, numStates);
  if(!isfinite(norm)) {
    // Won't be able to solve this equation, magnus has non-finite elements
    for(size_t n = 0; n < numStates; ++n) {
      statesFinal[n] = NaN;
    }
    return;
  }
  
  #ifdef LOW_ACCURACY_MATRIX_EXP
    const double lowAccuracyScale = 1.0e6;
    const double maxNorm3 = 1.495585217958292e-2 * lowAccuracyScale;
    const double maxNorm5 = 2.539398330063230e-1 * lowAccuracyScale;
    const double maxNorm7 = 9.504178996162932e-1 * lowAccuracyScale;
    const double maxNorm9 = 2.097847962157068e0 * lowAccuracyScale;
    const double maxNorm13 = 5.371920351148152e0 * lowAccuracyScale;
  #else
    // these are threshold norms for various orders of the Padé approximant
    const double maxNorm3 = 1.495585217958292e-2;
    const double maxNorm5 = 2.539398330063230e-1;
    const double maxNorm7 = 9.504178996162932e-1;
    const double maxNorm9 = 2.097847962157068e0;
    const double maxNorm13 = 5.371920351148152e0;
  #endif

  if(norm < maxNorm3) {
    //const size_t numSquares = scaleMagnus(norm, maxNorm3);
    const size_t numSquares = 0;
    //cout << "Order  3: numSquares = " << numSquares << '\n';
    pade3Solve(statesFinal, statesInitial, numSquares);
  }
  else if(norm < maxNorm5) {
    //const size_t numSquares = scaleMagnus(norm, maxNorm5);
    const size_t numSquares = 0;
    //cout << "Order  5: numSquares = " << numSquares << '\n';
    pade5Solve(statesFinal, statesInitial, numSquares);
  }
  else if(norm < maxNorm7) {
    //const size_t numSquares = scaleMagnus(norm, maxNorm7);
    const size_t numSquares = 0;
    //cout << "Order  7: numSquares = " << numSquares << '\n';
    pade7Solve(statesFinal, statesInitial, numSquares);
  }
  else if(norm < maxNorm9) {
    //const size_t numSquares = scaleMagnus(norm, maxNorm9);
    const size_t numSquares = 0;
    //cout << "Order  9: numSquares = " << numSquares << '\n';
    pade9Solve(statesFinal, statesInitial, numSquares);
  }
  else {
    const size_t numSquares = scaleMagnus(norm, maxNorm13);
    //cout << "Order 13: numSquares = " << numSquares << '\n'; cout.flush();
    pade13Solve(statesFinal, statesInitial, numSquares);
  }
  
  ++numMatrixExponentials;
}


inline int
DiffEqFunctor::scaleMagnus(const double & norm, const double & maxNorm)
{
  const double normRatio = norm / maxNorm;
  if(normRatio <= 1.0) {
    return 0;
  }
  
  const int numSquares = (int) ceil(log2(normRatio) );
  const double scale = fastPow(0.5, numSquares);
  magnus *= scale;
  return numSquares;
}


inline void
DiffEqFunctor::pade3Solve(VectorXd & statesFinal,
                          const Eigen::VectorXd & statesInitial,
                          size_t numSquares)
{
  // solve statesFinal = exp(magnus) * states1
  //   via Padé approximant of order 3/3
  
  // pN_j = arbitrary constant * (2N-j)! N! / ((2N)! (N-j)! j!)
  const double p3_0 = 120.0 / 120.0;
  const double p3_1 = 60.0 / 120.0;
  const double p3_2 = 12.0 / 120.0;
  const double p3_3 = 1.0 / 120.0;

  // compute necessary powers of magnus
  magnus2.noalias() = magnus * magnus;
  
  // store even terms
  evenTerms = p3_2 * magnus2;
  evenTerms.diagonal().array() += p3_0;
  
  // store m^-1 * odd terms
  oddTerms = p3_3 * magnus2;
  oddTerms.diagonal().array() += p3_1;
  
  if(numSquares == 0) {
    // With even and odd terms of approximation to exp(magnus) computed,
    //  solve for statesFinal
    padeSolveSolution(statesFinal, statesInitial);
  }
  else {
    // Have even and odd terms of approximation to root of exp(magnus).
    // Construct full approximation by squaring, then solve for statesFinal
    padeSolveSolutionWithSquares(statesFinal, statesInitial, numSquares);
  }
}



inline void
DiffEqFunctor::pade5Solve(VectorXd & statesFinal,
                          const Eigen::VectorXd & statesInitial,
                          size_t numSquares)
{
  // solve statesFinal = exp(magnus) * states1
  //   via Padé approximant of order 5/5
  
  // pN_j = arbitrary constant * (2N-j)! N! / ((2N)! (N-j)! j!)
  const double p5_0 = 30240.0 / 30240.0;
  const double p5_1 = 15120.0 / 30240.0;
  const double p5_2 = 3360.0 / 30240.0;
  const double p5_3 = 420.0 / 30240.0;
  const double p5_4 = 30.0 / 30240.0;
  const double p5_5 = 1.0 / 30240.0;

  // compute necessary powers of magnus
  magnus2.noalias() = magnus * magnus;
  magnus4.noalias() = magnus2 * magnus2;
  
  // store even terms
  evenTerms = p5_4 * magnus4 + p5_2 * magnus2;
  evenTerms.diagonal().array() += p5_0;
  
  // store m^-1 * odd terms
  oddTerms = p5_5 * magnus4 + p5_3 * magnus2;
  oddTerms.diagonal().array() += p5_1;
  
  if(numSquares == 0) {
    // With even and odd terms of approximation to exp(magnus) computed,
    //  solve for statesFinal
    padeSolveSolution(statesFinal, statesInitial);
  }
  else {
    // Have even and odd terms of approximation to root of exp(magnus).
    // Construct full approximation by squaring, then solve for statesFinal
    padeSolveSolutionWithSquares(statesFinal, statesInitial, numSquares);
  }
}



inline void
DiffEqFunctor::pade7Solve(VectorXd & statesFinal,
                          const Eigen::VectorXd & statesInitial,
                          size_t numSquares)
{
  // solve statesFinal = exp(magnus) * states1
  //   via Padé approximant of order 7/7
  
  // pN_j = arbitrary constant * (2N-j)! N! / ((2N)! (N-j)! j!)
  const double p7_0 = 17297280.0 / 17297280.0;
  const double p7_1 = 8648640.0 / 17297280.0;
  const double p7_2 = 1995840.0 / 17297280.0;
  const double p7_3 = 277200.0 / 17297280.0;
  const double p7_4 = 25200.0 / 17297280.0;
  const double p7_5 = 1512.0 / 17297280.0;
  const double p7_6 = 56.0 / 17297280.0;
  const double p7_7 = 1.0 / 17297280.0;

  // compute necessary powers of magnus
  magnus2.noalias() = magnus * magnus;
  magnus4.noalias() = magnus2 * magnus2;
  magnus6.noalias() = magnus4 * magnus2;
  
  // store even terms
  evenTerms = p7_6 * magnus6 + p7_4 * magnus4 + p7_2 * magnus2;
  evenTerms.diagonal().array() += p7_0;
  
  // store m^-1 * odd terms
  oddTerms = p7_7 * magnus6 + p7_5 * magnus4 + p7_3 * magnus2;
  oddTerms.diagonal().array() += p7_1;
  
  if(numSquares == 0) {
    // With even and odd terms of approximation to exp(magnus) computed,
    //  solve for statesFinal
    padeSolveSolution(statesFinal, statesInitial);
  }
  else {
    // Have even and odd terms of approximation to root of exp(magnus).
    // Construct full approximation by squaring, then solve for statesFinal
    padeSolveSolutionWithSquares(statesFinal, statesInitial, numSquares);
  }
}



inline void
DiffEqFunctor::pade9Solve(VectorXd & statesFinal,
                          const Eigen::VectorXd & statesInitial,
                          size_t numSquares)
{
  // solve statesFinal = exp(magnus) * states1
  //   via Padé approximant of order 9/9
  
  // pN_j = arbitrary constant * (2N-j)! N! / ((2N)! (N-j)! j!)
  const double p9_0 = 17643225600.0 / 17643225600.0;
  const double p9_1 = 8821612800.0 / 17643225600.0;
  const double p9_2 = 2075673600.0 / 17643225600.0;
  const double p9_3 = 302702400.0 / 17643225600.0;
  const double p9_4 = 30270240.0 / 17643225600.0;
  const double p9_5 = 2162160.0 / 17643225600.0;
  const double p9_6 = 110880.0 / 17643225600.0;
  const double p9_7 = 3960.0 / 17643225600.0;
  const double p9_8 = 90.0 / 17643225600.0;
  const double p9_9 = 1.0 / 17643225600.0;

  // compute necessary powers of magnus
  magnus2.noalias() = magnus * magnus;
  magnus4.noalias() = magnus2 * magnus2;
  magnus6.noalias() = magnus4 * magnus2;
  magnus8.noalias() = magnus4 * magnus4;
  
  // store even terms
  evenTerms = p9_8 * magnus8 + p9_6 * magnus6 + p9_4 * magnus4 +
              p9_2 * magnus2;
  evenTerms.diagonal().array() += p9_0;
  
  // store m^-1 * odd terms
  oddTerms = p9_9 * magnus8 + p9_7 * magnus6 + p9_5 * magnus4 +
             p9_3 * magnus2;
  oddTerms.diagonal().array() += p9_1;
  
  if(numSquares == 0) {
    // With even and odd terms of approximation to exp(magnus) computed,
    //  solve for statesFinal
    padeSolveSolution(statesFinal, statesInitial);
  }
  else {
    // Have even and odd terms of approximation to root of exp(magnus).
    // Construct full approximation by squaring, then solve for statesFinal
    padeSolveSolutionWithSquares(statesFinal, statesInitial, numSquares);
  }
}



inline void
DiffEqFunctor::pade13Solve(VectorXd & statesFinal,
                           const Eigen::VectorXd & statesInitial,
                           size_t numSquares)
{
  // solve statesFinal = exp(magnus) * states1
  //   via Padé approximant of order 13/13
  
  // pN_j = arbitrary constant * (2N-j)! N! / ((2N)! (N-j)! j!)
  const double p13_0 = 64764752532480008.0 / 64764752532480008.0;
  const double p13_1 = 32382376266240004.0 / 64764752532480008.0;
  const double p13_2 = 7771770303897600.0 / 64764752532480008.0;
  const double p13_3 = 1187353796428800.0 / 64764752532480008.0;
  const double p13_4 = 129060195264000.0 / 64764752532480008.0;
  const double p13_5 = 10559470521600.0 / 64764752532480008.0;
  const double p13_6 = 670442572800.0 / 64764752532480008.0;
  const double p13_7 = 33522128640.0 / 64764752532480008.0;
  const double p13_8 = 1323241920.0 / 64764752532480008.0;
  const double p13_9 = 40840800.0 / 64764752532480008.0;
  const double p13_10 = 960960.0 / 64764752532480008.0;
  const double p13_11 = 16380.0 / 64764752532480008.0;
  const double p13_12 = 182.0 / 64764752532480008.0;
  const double p13_13 = 1.0 / 64764752532480008.0;

  // compute necessary powers of m
  magnus2.noalias() = magnus * magnus;
  magnus4.noalias() = magnus2 * magnus2;
  magnus6.noalias() = magnus4 * magnus2;
  
  // temporarily store m^-6 * high-power even terms in numerator
  numerator.noalias() = p13_12 * magnus6 + p13_10 * magnus4 + p13_8 * magnus2;
  // temporarily store m^-6 * high-power odd terms in denominator
  denominator.noalias()= p13_13 * magnus6 + p13_11 * magnus4 + p13_9 * magnus2;
  
  // store even terms
  evenTerms.noalias() = magnus6 * numerator
                      + p13_6 * magnus6 + p13_4 * magnus4 + p13_2 * magnus2;
  evenTerms.diagonal().array() += p13_0;
  
  // store m^-1 * odd terms
  oddTerms.noalias() = magnus6 * denominator;
                     + p13_7 * magnus6 + p13_5 * magnus4 + p13_3 * magnus2;
  oddTerms.diagonal().array() += p13_1;
  
  if(numSquares == 0) {
    // With even and odd terms of approximation to exp(magnus) computed,
    //  solve for statesFinal
    padeSolveSolution(statesFinal, statesInitial);
  }
  else {
    // Have even and odd terms of approximation to root of exp(magnus).
    // Construct full approximation by squaring, then solve for statesFinal
    padeSolveSolutionWithSquares(statesFinal, statesInitial, numSquares);
  }
}



inline void
DiffEqFunctor::padeSolveSolution(VectorXd & statesFinal,
                                 const Eigen::VectorXd & statesInitial)
{
  // With even and odd terms of approximation to exp(magnus) computed,
  //  solve for statesFinal
  
  // compute full odd terms, temporarily store in numerator
  numerator.noalias() = oddTerms * magnus;
  
  // compute denominator by subtracting even - odd
  denominator = evenTerms - numerator;
  // computer numerator by adding even + odd
  numerator += evenTerms;
  
  // compute numeratorTimesStates1 = numerator * states1
  //   (note, extra fiddling required because magnus is an augmented matrix,
  //    but states1 and statesFinal don't have extra 1.0 appended to them)
  numeratorTimesStates1 =
     numerator.topLeftCorner(numStates, numStates) * statesInitial
   + numerator.col(numStates).head(numStates)
   - denominator.col(numStates).head(numStates);
  
  // solve denominator * statesFinal = numeratorTimesStates1 for statesFinal
  #ifdef USE_LU_DECOMP
    statesFinal = denominator.topLeftCorner(numStates, numStates).
                    partialPivLu().solve(numeratorTimesStates1);
  #elif USE_QR_DECOMP
    statesFinal = denominator.topLeftCorner(numStates, numStates).
                    colPivHouseholderQr().solve(numeratorTimesStates1);
  #else
    #error Must define either USE_LU_DECOMP or USE_QR_DECOMP
  #endif
}



inline void
DiffEqFunctor::padeSolveSolutionWithSquares(VectorXd & statesFinal,
                                          const Eigen::VectorXd & statesInitial,
                                          size_t numSquares)
{
  // Have even and odd terms of approximation to root of exp(magnus).
  // Construct full approximation by squaring, then solve for statesFinal
  
    // compute full odd terms, temporarily store in numerator
  numerator.noalias() = oddTerms * magnus;
  
  // compute denominator by subtracting even - odd
  denominator = evenTerms - numerator;
  // computer numerator by adding even + odd
  numerator += evenTerms;
  
  //if(numSquares <= numStates) {
  if(false) {
    // few squares to perform, so square denominator and numerator separately,
    // and solve as usual
    
    // square numerator and denominator the requisite number of times
    while(numSquares != 0) {
      magnus2.noalias() = denominator * denominator;
      denominator.swap(magnus2);
      
      magnus2.noalias() = numerator * numerator;
      numerator.swap(magnus2);
      --numSquares;
    }
    
    // compute numeratorTimesStates1 = numerator * states1
    //   (note, extra fiddling required because magnus is an augmented matrix,
    //    but states1 and statesFinal don't have extra 1.0 appended to them)
    numeratorTimesStates1 =
       numerator.topLeftCorner(numStates, numStates) * statesInitial
     + numerator.rightCols<1>().head(numStates)
     - denominator.rightCols<1>().head(numStates);
    
    // solve denominator * statesFinal = numeratorTimesStates1 for statesFinal
    #ifdef USE_LU_DECOMP
      statesFinal = denominator.topLeftCorner(numStates, numStates).
                      partialPivLu().solve(numeratorTimesStates1);
    #elif USE_QR_DECOMP
      statesFinal = denominator.topLeftCorner(numStates, numStates).
                      colPivHouseholderQr().solve(numeratorTimesStates1);
    #else
      #error Must define either USE_LU_DECOMP or USE_QR_DECOMP
    #endif
  }
  else {
    // many squares to perform, so construct full matrix exponential and
    // repeatedly square it
    
    // solve for exp(magnus), store result in magnus2
    #ifdef USE_LU_DECOMP
      magnus2 = denominator.partialPivLu().solve(numerator);
    #elif USE_QR_DECOMP
      magnus2 = denominator.colPivHouseholderQr().solve(numerator);
    #else
      #error Must define either USE_LU_DECOMP or USE_QR_DECOMP
    #endif
    
    // square the matrix exponential the requisite number of times
    while(numSquares != 0) {
      magnus4.noalias() = magnus2 * magnus2;
      magnus2.swap(magnus4);
      --numSquares;
    }
    
    // compute statesFinal by matrix multiplication
    statesFinal = magnus2.topLeftCorner(numStates, numStates) * statesInitial
                + magnus2.rightCols<1>().head(numStates);
  }
}
#endif

#ifdef USE_SPARSE_MATRICES
void
DiffEqFunctor::expSolveTaylor(VectorXd & statesFinal,
                              const VectorXd & statesInitial)
{
  // subtract off the average eigenvector (may aid in numerical convergence)
  #ifdef SUB_MEAN_EIGEN
    const double meanEigenvalue
      = magnusMat.diagonal().sum() / (double)(numStates + 1);
    magnusMat.diagonal().array() -= meanEigenvalue;
    // compute the L1 norm of magnus
    double normMagnus = std::max(induced1Norm(magnusMat),
                             magnusVec.lpNorm<1>() + abs(meanEigenvalue));
  #else
    // compute the L1 norm of magnus
    double normMagnus = std::max(induced1Norm(magnusMat),
                                 magnusVec.lpNorm<1>());
  #endif

  if(normMagnus == 0) {
    // if the L1 norm is zero, statesFinal = statesInitial, so we're done
    return;
  }
  // put upper bound on largest eigenvalue
  #ifdef SUB_MEAN_EIGEN
    magnusBR = -meanEigenValue;
    const double infNormMagnus = std::max(abs(magnusBR),
                                         inducedInfNorm(magnusMat, magnusVec));
  #else
    magnusBR = 0;
    const double infNormMagnus = inducedInfNorm(magnusMat, magnusVec);
  #endif
  const double maxEigen = sqrt(normMagnus * infNormMagnus);
  
  // determine the number of intervals (e.g. number of matrix exponentials)
  const double maxAllowedEigenvalue = 30.0;
  const size_t numIntervals = (size_t)ceil(maxEigen / maxAllowedEigenvalue);
  if(numIntervals > 1) {
    magnusMat /= (double)numIntervals;
    magnusVec /= (double)numIntervals;
    #ifdef SUB_MEAN_EIGEN
      magnusBR /= (double)numIntervals;
    #endif
    // Why is the following line wrong???
    //normMagnus /= (double)numIntervals;
  }

  #ifdef LOW_ACCURACY_MATRIX_EXP
    double stopTol = 0.1 * maxErr / (double)numIntervals;
  #else
    double stopTol = 1.1e-16;
  #endif  

  // choose the optimal number of steps and degree of truncated Taylor series
  //   to solve the problem
  size_t numStepsPerInterval, taylorDegree;
  tie(numStepsPerInterval, taylorDegree)
    = getTruncatedParameters(normMagnus, stopTol);
  if(taylorDegree < 1) {
    // expansion is trivial
    return;
  }

  statesFinal = statesInitial;
  for(size_t interval = 0; interval < numIntervals; ++interval) {
    expSolveTaylorInner(statesFinal, stopTol, numStepsPerInterval,
                        taylorDegree);
  }
}
#else
void
DiffEqFunctor::expSolveTaylor(VectorXd & statesFinal,
                              const VectorXd & statesInitial)
{
  // subtract off the average eigenvector (may aid in numerical convergence)
  #ifdef SUB_MEAN_EIGEN
    const double meanEigenvalue = magnus.diagonal().mean();
    magnus.diagonal().head(numStates).array() -= meanEigenvalue;
    magnusBR = -meanEigenvalue;
  #endif

  // compute the L1 norm of magnus
  double normMagnus = induced1Norm(magnus);
  
  if(normMagnus == 0) {
    // if the L1 norm is zero, statesFinal = statesInitial, so we're done
    return;
  }
  // put upper bound on largest eigenvalue
  const double maxEigen = sqrt(normMagnus * inducedInfNorm(magnus));
  
  // determine the number of intervals (e.g. number of matrix exponentials)
  const double maxAllowedEigenValue = 30.0;
  const size_t numIntervals = (size_t)ceil(maxEigen / maxAllowedEigenValue);
  
  if(numIntervals > 1) {
    magnus /= (double)numIntervals;
    magnusBR /= (double)numIntervals;
    // Why is the following line wrong???
    //normMagnus /= (double)numIntervals;
  }

  #ifdef LOW_ACCURACY_MATRIX_EXP
    double stopTol = 0.1 * maxErr / (double)numIntervals;
  #else
    double stopTol = 1.1e-16;
  #endif  

  // choose the optimal number of steps and degree of truncated Taylor series
  //   to solve the problem
  size_t numStepsPerInterval, taylorDegree;
  tie(numStepsPerInterval, taylorDegree)
    = getTruncatedParameters(normMagnus, stopTol);
  if(taylorDegree < 1) {
    // expansion is trivial
    return;
  }

  statesFinal = statesInitial;
  for(size_t interval = 0; interval < numIntervals; ++interval) {
    expSolveTaylorInner(statesFinal, stopTol, numStepsPerInterval,
                        taylorDegree);
  }
}
#endif

#ifdef USE_SPARSE_MATRICES
void
DiffEqFunctor::expSolveTaylorInner(VectorXd & statesFinal,
                                   const double & stopTol,
                                   const size_t & numSteps,
                                   const size_t & taylorDegree)
{
  // solve inhomogenous equation with augmented matrix:
  //   [statesFinal, 1] = exp(magnus) * [statesInitial, 1]
  // via Truncated Taylor Series
  
  for(size_t step = 0; step < numSteps; ++step) {
    // Compute statesFinal *= exp(magnus/numSteps) for each step
    
    // Update working states to the value of each term in the series,
    //  keeping track of the norms of the last two terms in the expansion
    taylorTerm = statesFinal;  // first term is multiplied by identity
    double c1 = max(taylorTerm.lpNorm<Infinity>(), 1.0); // norm of term
    
    taylorTerm = (magnusMat * taylorTerm + magnusVec) / (double)(numSteps);

    #ifdef SUB_MEAN_EIGEN
      // now the inhomogenous (ih) term in magnus is non-zero, and must be
      // kept track of
      double ihFact = 1.0;
      double ihTerm = magnusBR / (double)(numSteps);
      ihFact += ihTerm;
      double c2 = max(taylorTerm.lpNorm<Infinity>(),
                      abs(ihTerm)); // norm of last term
    #else
      double c2 = taylorTerm.lpNorm<Infinity>(); // norm of last term
    #endif
    
    statesFinal += taylorTerm;
    
    for(size_t term = 2; term <= taylorDegree; ++term) {
      #ifdef SUB_MEAN_EIGEN
        if(c1 + c2 < stopTol * max(statesFinal.lpNorm<Infinity>(),
                                   abs(ihFact))) {
          // taylorTerm has shrunk enough that series has converged
          break;
        }
        const double factorialFact = (double)(term * numSteps);
        taylorTerm =
          ( magnusMat * taylorTerm + magnusVec * ihTerm).eval()
            / factorialFact;
        ihTerm *= magnusBR / factorialFact;
        ihFact += ihTerm;
        statesFinal += taylorTerm;
        c1 = c2;
        c2 = max(taylorTerm.lpNorm<Infinity>(), abs(ihTerm));
      #else
        if(c1 + c2 < stopTol * max(statesFinal.lpNorm<Infinity>(), 1.0)) {
          // taylorTerm has shrunk enough that series has converged
          break;
        }
        taylorTerm = ( magnusMat * taylorTerm ).eval()
            / (double) (term * numSteps);
        statesFinal += taylorTerm;
        c1 = c2;
        c2 = taylorTerm.lpNorm<Infinity>();
      #endif
    }
    #ifdef SUB_MEAN_EIGEN
      statesFinal /= ihFact;
    #endif
  }
  ++numMatrixExponentials;
}
#else
void
DiffEqFunctor::expSolveTaylorInner(VectorXd & statesFinal,
                                   const double & stopTol,
                                   const size_t & numSteps,
                                   const size_t & taylorDegree)
{
  // solve inhomogenous equation with augmented matrix:
  //   [statesFinal, 1] = exp(magnus) * [statesInitial, 1]
  // via Truncated Taylor Series
  
  for(size_t step = 0; step < numSteps; ++step) {
    // Compute statesFinal *= exp(magnus/numSteps) for each step
    
    // Update working states to the value of each term in the series,
    //  keeping track of the norms of the last two terms in the expansion
    taylorTerm = statesFinal;  // first term is multiplied by identity
    double c1 = max(taylorTerm.lpNorm<Infinity>(), 1.0); // norm of term
    
    taylorTerm = ( magnus.topLeftCorner(numStates, numStates) * taylorTerm
                 + magnus.rightCols<1>().head(numStates)).eval()
                   / (double)(numSteps);
    #ifdef SUB_MEAN_EIGEN
      // now the inhomogenous (ih) term in magnus is non-zero, and must be
      // kept track of
      double ihFact = 1.0;
      double ihTerm = magnusBR / (double)(numSteps);
      ihFact += ihTerm;
      double c2 = max(taylorTerm.lpNorm<Infinity>(),
                      abs(ihTerm)); // norm of last term
    #else
      double c2 = taylorTerm.lpNorm<Infinity>(); // norm of last term
    #endif
    
    statesFinal += taylorTerm;
    
    for(size_t term = 2; term <= taylorDegree; ++term) {
      #ifdef SUB_MEAN_EIGEN
        if(c1 + c2 < stopTol * max(statesFinal.lpNorm<Infinity>(),
                                   abs(ihFact))) {
          // taylorTerm has shrunk enough that series has converged
          break;
        }
        const double factorialFact = (double)(term * numSteps);
        taylorTerm =
          ( magnus.topLeftCorner(numStates, numStates) * taylorTerm
          + magnus.rightCols<1>().head(numStates) * ihTerm).eval()
            / factorialFact;
        ihTerm *= magnusBR / factorialFact;
        ihFact += ihTerm;
        statesFinal += taylorTerm;
        c1 = c2;
        c2 = max(taylorTerm.lpNorm<Infinity>(), abs(ihTerm));
      #else
        if(c1 + c2 < stopTol * max(statesFinal.lpNorm<Infinity>(), 1.0)) {
          // taylorTerm has shrunk enough that series has converged
          break;
        }
        taylorTerm =
          ( magnus.topLeftCorner(numStates, numStates) * taylorTerm ).eval()
            / (double) (term * numSteps);
        statesFinal += taylorTerm;
        c1 = c2;
        c2 = taylorTerm.lpNorm<Infinity>();
      #endif
    }
    #ifdef SUB_MEAN_EIGEN
      statesFinal /= ihFact;
    #endif
  }
  ++numMatrixExponentials;
}
#endif

tuple<size_t, size_t>
DiffEqFunctor::getTruncatedParameters(const double & normMagnus,
                                      const double & stopTol)
{
  // Compute optimal number of steps and degree of truncated Taylor series
  
  // Maximum power of magnus to consider norm of:
  const size_t pMax = 8;
  // Maximum length of truncated Taylor series
  const size_t mMax = 100;
  
  // halfVec lists the values of theta_m for tol = 2^-10, 1 <= m <= 100
  const double halfVec [] = {
    0.00195058, 0.0744366, 0.266455, 0.524205, 0.810269, 1.10823, 1.41082,
    1.71493, 2.01903, 2.32247, 2.62492, 2.9263, 3.22657, 3.52578, 3.82398,
    4.12123, 4.41759, 4.71314, 5.00792, 5.302, 5.59543, 5.88825, 6.18051,
    6.47224, 6.76348, 7.05427, 7.34464, 7.6346, 7.92419, 8.21342, 8.50232,
    8.79091, 9.07921, 9.36722, 9.65497, 9.94247, 10.2297, 10.5168, 10.8036,
    11.0902, 11.3766, 11.6629, 11.9489, 12.2348, 12.5205, 12.8061, 13.0915,
    13.3768, 13.6619, 13.9469, 14.2318, 14.5165, 14.8012, 15.0857, 15.3702,
    15.6547, 15.9391, 16.2235, 16.5078, 16.792, 17.0761, 17.3599, 17.6437,
    17.9273, 18.2108, 18.4943, 18.7776, 19.0609, 19.3442, 19.6273, 19.9104,
    20.1935, 20.4764, 20.7593, 21.0422, 21.325, 21.6077, 21.8904, 22.173,
    22.4556, 22.7381, 23.0206, 23.303, 23.5854, 23.8677, 24.15, 24.4322,
    24.7144, 24.9966, 25.2787, 25.5608, 25.8428, 26.1248, 26.4068, 26.6887,
    26.9706, 27.2524, 27.5342, 27.816, 28.0978
  };
  
  // doubleVec lists the values of theta_m for tol = 2^-53, 1 <= m <= 100
  const double doubleVec [] = {
    2.22045e-16, 2.58096e-08, 1.38635e-05, 0.000339717, 0.00240088, 0.00906566,
    0.0238446, 0.0499123, 0.0895776, 0.144183, 0.214236, 0.299616, 0.399778,
    0.513915, 0.641084, 0.780287, 0.930533, 1.09086, 1.26038, 1.43825, 1.62372,
    1.81608, 2.01471, 2.21905, 2.42858, 2.64285, 2.86145, 3.084, 3.31017,
    3.53967, 3.77221, 4.00756, 4.2455, 4.48582, 4.72835, 4.97292, 5.21938,
    5.46759, 5.71744, 5.9688, 6.22158, 6.47568, 6.73102, 6.9875, 7.24507,
    7.50365, 7.76317, 8.02359, 8.28485, 8.5469, 8.80969, 9.07319, 9.33734,
    9.60212, 9.8675, 10.1334, 10.3999, 10.6669, 10.9343, 11.2022, 11.4705,
    11.7392, 12.0084, 12.2778, 12.5477, 12.8178, 13.0883, 13.3591, 13.6302,
    13.9016, 14.1732, 14.4451, 14.7172, 14.9896, 15.2622, 15.535, 15.808,
    16.0812, 16.3546, 16.6281, 16.9019, 17.1758, 17.4498, 17.724, 17.9984,
    18.2729, 18.5476, 18.8223, 19.0972, 19.3723, 19.6474, 19.9227, 20.198,
    20.4735, 20.7491, 21.0248, 21.3005, 21.5764, 21.8523, 22.1284
  };
  
  double const* const theta = log2(stopTol) >= -10.0 ? halfVec : doubleVec;
  double const* const thetaEnd = theta + mMax;
  const double & thetaBack = theta[mMax - 1];
  
  /*
    // example of interpolation between doubleVec and halfVec
    const double logTol = log2(stopTol);
    const double interpDouble = (10.0 + logTol) / (10.0 - 53.0);
    const double interpHalf = 1.0 - interpDouble;
    double const* dVecPtr = doubleVec;
    double const* hVecPtr = halfVec;
    for(double & t : theta) {
      // t = *dVecPtr * interpDouble + *hVecPtr * interpHalf;
      t = exp(log(*dVecPtr) * interpDouble + log(*hVecPtr) * interpHalf);
      ++dVecPtr;
      ++hVecPtr;
    }
  */
  
  // Compute numSteps and taylorDegree
  size_t numSteps = -1, taylorDegree = -1;
  if( normMagnus <= 4.0 * thetaBack * (double)(pMax * (pMax + 3))
                    / (double)mMax ) {
    // If the norm is small enough, can compute numSteps and taylorDegree using
    // a simpler algorithm
    size_t m = 1;
    double bestCost = Inf;
    for(double const* thetaPtr = theta; thetaPtr != thetaEnd; ++thetaPtr) {
      const double steps_m = ceil(normMagnus / *thetaPtr);
      const double cost_m = steps_m * (double)m;
      if(cost_m < bestCost) {
        bestCost = cost_m;
        numSteps = (size_t)steps_m;
        taylorDegree = m;
      }
      ++m;
    }
  }
  else {
    // Must use more costly algorithm
    double bestCost = Inf;
    double prevNorm = estimatePNorm(2);
    for(size_t p = 2; p <= pMax; ++p) {
      // compute alpha_p, a norm on magnus that provides a bound on cost
      const double nextNorm = estimatePNorm(p + 1);
      const double alpha_p = max(prevNorm, nextNorm);
      prevNorm = nextNorm;
      
      const size_t mMin = p * (p-1) - 1;
      auto thetaPtr = theta + (mMin - 1);
      for(size_t m = mMin; m <= mMax; ++m) {
        const double steps_m = ceil(alpha_p / *thetaPtr);
        const double cost_m = steps_m * (double)m;
        if(cost_m < bestCost) {
          bestCost = cost_m;
          numSteps = (size_t) steps_m;
          taylorDegree = m;
        }
        ++thetaPtr;
      }
    }
  }
  
  return make_tuple(numSteps, taylorDegree);
}

#ifdef USE_SPARSE_MATRICES
inline void
assignColwiseSumTranspose(const Matrix_t & mat, VectorXd & v)
{
  // return the sum of mat's columns, as a transposed vector
  v.resize(mat.cols(), 1);
  for(int col = 0; col < mat.outerSize(); ++col) {
    double & sum = v[col];
    sum = 0.0;
    for(Matrix_t::InnerIterator matItr(mat, col); matItr; ++matItr) {
      sum += matItr.value();
    }
  }
}

inline double
DiffEqFunctor::estimatePNorm(const size_t & p)
{
  // Estimate the p-th root of the 1-norm of magnus^p, for an integer p
  // Assumes that:
  //   magnus is real
  //   only one column test-vector is used
  //   p is an integer >= 1
  // Reference:
  //  Nicholas J. Higham and Fran\c{c}oise Tisseur, "A Block Algorithm for
  //    Matrix 1-Norm Estimation with an Application to 1-Norm Pseudospectra",
  //    SIAM J. Matrix Anal. App. 21, 1185-1201, 2000. 
  // Adapted from normest1.m
  

  if(isNonNegative(magnusMat) && (magnusVec.array() >= 0).all()) {
    // The augmented matrix has all non-negative entries, so use simpler
    // estimate
    assignColwiseSumTranspose(magnusMat, testVector);
    #ifdef SUB_MEAN_EIGEN
      testVectorB = magnusVec.sum() + magnusBR;
      for(size_t j = 1; j < p; ++j) {
        testVector = (magnusMat * testVector).eval()
                   + magnusVec * testVectorB;
        testVectorB *= magnusBR;
      }
      return pow(std::max(testVector.lpNorm<Infinity>(),
                          abs(testVectorB)),
                 1.0/(double)p);
    #else
      for(size_t j = 1; j < p; ++j) {
        testVector = (magnusMat * testVector).eval();
      }
      return pow(testVector.lpNorm<Infinity>(), 1.0/(double)p);
    #endif
      
  }
  else if(magnusMat.rows() <= 3) {  // (should be 3 because of extra row)
    #ifdef SUB_MEAN_EIGEN
      throw DifferentialEquationException("Unimplemented");
    #endif
    /*
    testMatrix = magnus;
    for(size_t j = 1; j < p; ++j) {
      testMatrix = (magnus * testMatrix).eval();
    }
    return pow(induced1Norm(magnus), 1.0/(double)p);
    */
    return NaN;
  }
  
  // Not an easy case, will have to use block 1-norm power method
  // Take multiple passes; in each pass a different simple test vector is
  //   acted on by magnus^p and the 1-norm of the resulting vector is used
  //   as an estimate of the 1-norm of magnus^p
  // For first pass, use a ones(numRows, 1) / numRows as the test vector
  // For subsequent passes, use different unit vectors
  // After each pass, run a variety of tests to determine if it is time to
  //   terminate and choose the best current estimate
  
  const size_t maxPasses = 6; // maximum number of passes to make
  
  // First pass:
  size_t passNum = 0;
  double normEstimate = 0.0;
  double oldNormEstimate;
  size_t unitIndex;
  set<size_t> indexHistory;
  int nonNegBot;
  double pNonNegBot;
  // Since testVector is so simple, assign magnus * testVector
  assignColwiseSumTranspose(magnusMat, testVector);
  testVector /= (double)(numStates + 1);
  #ifdef SUB_MEAN_EIGEN
    testVectorB = (magnusVec.sum() + magnusBR) / (double)(numStates + 1);
  #else
    (void)nonNegBot;
    (void)pNonNegBot;
  #endif
  while(true) {
    oldNormEstimate = normEstimate;
    #ifdef SUB_MEAN_EIGEN
      for(size_t j = 1; j < p; ++j) {
        testVector = (magnusMat * testVector).eval()
                   + magnusVec * testVectorB;
        testVectorB *= magnusBR;
      }
      normEstimate = testVector.lpNorm<1>() + abs(testVectorB);
    #else
      for(size_t j = 1; j < p; ++j) {
        testVector = (magnusMat * testVector).eval();
      }
      normEstimate = testVector.lpNorm<1>();
    #endif

    // increment the number of passes and test for convergence
    ++passNum;
    if(passNum > 1 && normEstimate <= oldNormEstimate) {
      // if normEstimate decreases, take the larger, saved value and quit
      return oldNormEstimate;
    }
    
    if(passNum >= maxPasses) {
      // reached the maximum number of passes, choose the current estimate
      return normEstimate;
    }
    
    // Test changes to non-negativity of testVector
    if(passNum == 1) {
      nonNegVector = mySign(testVector);
      #ifdef SUB_MEAN_EIGEN
        nonNegBot = (int)(testVectorB >= 0);
      #endif
    }
    else {
      oldNonNegVec = nonNegVector;
      nonNegVector = mySign(testVector);
      #ifdef SUB_MEAN_EIGEN
        oldNonNegBot = nonNegBot;
        nonNegBot = (int)(testVectorB >= 0);
        if( abs(oldNonNegVec.dot(nonNegVector) + nonNegbot * oldNonNegBot)
            == magnusMat.rows() + 1 ) {
         // non-negativity of Y hasn't changed since last pass (or all changed)
          return normEstimate;
        }
      #else
        if( abs(oldNonNegVec.dot(nonNegVector)) == magnusMat.rows() ) {
         // non-negativity of Y hasn't changed since last pass (or all changed)
          return normEstimate;
        }
      #endif
    }

    // Set magnusPNonNeg = abs( magnus.transpose()^p * nonNegVector )
    #ifdef SUB_MEAN_EIGEN
      pNonNegBot = magnusVec.dot(nonNegVector.cast<double>())
                 + magnusBR * (double)nonNegBot;
      magnusPNonNeg = magnusMat.transpose() * nonNegVector.cast<double>();
      for(size_t j = 1; j < p; ++j) {
        pNonNegBot = magnusVec.dot(manusPNonNeg)) + magnusBR * pNonNegBot;
        magnusPNonNeg = (magnusMat.transpose() * magnusPNonNeg).eval();
      }
    #else
      pNonNegBot = magnusVec.dot(nonNegVector.cast<double>());
      magnusPNonNeg = magnusMat.transpose() * nonNegVector.cast<double>();
      for(size_t j = 1; j < p; ++j) {
        pNonNegBot = magnusVec.dot(magnusPNonNeg);
        magnusPNonNeg = (magnusMat.transpose() * magnusPNonNeg).eval();
      }
    #endif
    
    // Set the new unitIndex to point to the maximum value of magnusPNonNeg
    magnusPNonNeg.maxCoeff(&unitIndex);
    if(magnusPNonNeg[unitIndex] < pNonNegBot) {
      unitIndex = magnusPNonNeg.rows() + 1;
    }
    if(passNum > 1) {
      if(indexHistory.count(unitIndex)) {
        // unitIndex has repeated, so terminate
        return normEstimate;
      }
    }
    // Select the next testVector
    if(unitIndex < (size_t)magnusMat.cols()) {
      testVector = magnusMat.col(unitIndex);
    }
    else {
      testVector = magnusVec;
    }
    
    // Record that this unitIndex has been used
    indexHistory.insert(unitIndex);
  }
}
#else
inline double
DiffEqFunctor::estimatePNorm(const size_t & p)
{
  // Estimate the p-th root of the 1-norm of magnus^p, for an integer p
  // Assumes that:
  //   magnus is real
  //   only one column test-vector is used
  //   p is an integer >= 1
  // Reference:
  //  Nicholas J. Higham and Fran\c{c}oise Tisseur, "A Block Algorithm for
  //    Matrix 1-Norm Estimation with an Application to 1-Norm Pseudospectra",
  //    SIAM J. Matrix Anal. App. 21, 1185-1201, 2000. 
  // Adapted from normest1.m
  
  if((magnus.array() >= 0).all()) {
    // The matrix has all non-negative entries, so use simpler estimate
    testVector = magnus.colwise().sum().transpose();
    for(size_t j = 1; j < p; ++j) {
      testVector = (magnus * testVector).eval();
    }
    return pow(testVector.lpNorm<Infinity>(), 1.0/(double)p);
  }
  else if(magnus.rows() <= 4) {
    testMatrix = magnus;
    for(size_t j = 1; j < p; ++j) {
      testMatrix = (magnus * testMatrix).eval();
    }
    return pow(induced1Norm(magnus), 1.0/(double)p);
  }
  
  // Not an easy case, will have to use block 1-norm power method
  // Take multiple passes; in each pass a different simple test vector is
  //   acted on by magnus^p and the 1-norm of the resulting vector is used
  //   as an estimate of the 1-norm of magnus^p
  // For first pass, use a ones(numRows, 1) / numRows as the test vector
  // For subsequent passes, use different unit vectors
  // After each pass, run a variety of tests to determine if it is time to
  //   terminate and choose the best current estimate
  
  const size_t maxPasses = 6; // maximum number of passes to make
  
  // First pass:
  size_t passNum = 0;
  double normEstimate = 0.0;
  double oldNormEstimate;
  size_t unitIndex;
  set<size_t> indexHistory;
  // Since testVector is so simple, assign magnus * testVector
  testVector = magnus.colwise().mean().transpose();
  
  while(true) {
    // apply magnus^(p-1) to testVector (it already has magnus^1 multiplied in)
    for(size_t j = 1; j < p; ++j) {
      testVector = (magnus * testVector).eval();
    }
    // compute the 1-norm of testVector, saving previous estimate
    oldNormEstimate = normEstimate;
    normEstimate = testVector.lpNorm<1>();
    
    // increment the number of passes and test for convergence
    ++passNum;
    if(passNum > 1 && normEstimate <= oldNormEstimate) {
      // if normEstimate decreases, take the larger, saved value and quit
      return oldNormEstimate;
    }
    
    if(passNum >= maxPasses) {
      // reached the maximum number of passes, choose the current estimate
      return normEstimate;
    }
    
    // Test changes to non-negativity of testVector
    if(passNum == 1) {
      nonNegVector = mySign(testVector);
    }
    else {
      oldNonNegVec = nonNegVector;
      nonNegVector = mySign(testVector);
      if( abs(oldNonNegVec.dot(nonNegVector)) == magnus.rows() ) {
        // non-negativity of Y hasn't changed since last pass (or all changed)
        return normEstimate;
      }
    }

    // Set magnusPNonNeg = abs( magnus.transpose()^p * nonNegVector )
    magnusPNonNeg = magnus.transpose() * nonNegVector.cast<double>();
    for(size_t j = 1; j < p; ++j) {
      magnusPNonNeg = (magnus.transpose() * magnusPNonNeg).eval();
    }
    magnusPNonNeg = magnusPNonNeg.array().abs();
    
    // Set the new unitIndex to point to the maximum value of magnusPNonNeg
    magnusPNonNeg.maxCoeff(&unitIndex);
    if(passNum > 1) {
      if(indexHistory.count(unitIndex)) {
        // unitIndex has repeated, so terminate
        return normEstimate;
      }
    }
    // Record that this unitIndex has been used
    indexHistory.insert(unitIndex);
    
    // Select the next testVector
    testVector = magnus.col(unitIndex);
  }
}
#endif
