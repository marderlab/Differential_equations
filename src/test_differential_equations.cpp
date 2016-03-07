#include <iostream>
#include "../include/differential_equations.h"
#include "../include/Timer.h"


using namespace std;
using namespace Eigen;



/////////////////////////////// Helper functions //////////////////////////////

template <class T>
inline T
abs(const T & x)
{
  return x >= 0 ? x : -x;
}



vector<double>
fitPolynomial(const vector<double> & x, const vector<double> & y,
              size_t degree)
{
  // fit a polynomial of specified degree to observations of y and x, and
  // return coefficients of the fit
  
  // sanity check on sizes of x and y
  size_t numX = x.size();
  size_t numY = y.size();
  if(numX != numY) {
    throw runtime_error("Size of x and y must match");
  }
  const size_t numCoefs = degree + 1;
  
  // convert STL vectors into Eigen vectors
  VectorXd xVec = Map<const VectorXd> (&x[0], numX, 1);
  VectorXd yVec = Map<const VectorXd> (&y[0], numX, 1);
  
  
  // create a matrix to hold powers of x
  MatrixXd xMat (numX, numCoefs);
  
  // fill out the values of xMat: xMat_nm = x_n ^ (degree + 1 - m)
  xMat.col(degree).setOnes();
  xMat.col(degree - 1).array() = xVec.array();
  for(size_t colNum = degree - 1; colNum > 0; --colNum) {
    xMat.col(colNum - 1) = xMat.col(colNum).array() * xVec.array();
  }

  // solve xMat * coefVec = yVec for coefVec
  VectorXd coefVec = xMat.colPivHouseholderQr().solve(yVec);

  // return STL vector  
  return vector<double> (coefVec.data(), coefVec.data() + numCoefs);
}



void
fitErrors(const vector<double> & logX, const vector<double> & logErrors,
          const string xLabel)
{
  auto xItr = logX.cbegin();
  const double xFront = *xItr;
  bool constX = true;
  for(++xItr; xItr != logX.cend(); ++xItr) {
    if(*xItr != xFront) {
      constX = false;
      break;
    }
  }
  
  if(!constX) {
    vector<double> coefs = fitPolynomial(logX, logErrors, 1);
    double num = pow(10.0, -coefs[1]/coefs[0]);
    cout << "Err ~ (" << xLabel << '/' << num << ")^" << coefs[0] << '\n';
  }
  else {
    double errMean = 0.0;
    for(const double & err : logErrors) {
      errMean += pow(10.0, err);
    }
    errMean /= logErrors.size();
    cout << "Err ~ " << errMean << '\n';
  }
}



/////////////////////////////// TestDiffEq Class //////////////////////////////

class TestDiffEq : public DiffEqFunctor
{
  public:
    // return the true error in the final states at requested accuracy goal
    double getError(const double & accuracyGoal);
    // return the true error in the recorded traces at requested accuracy goal
    double getTraceError(const double & accuracyGoal);
    // run comprehensive test
    void test(double sampleDT=0.0);

    // this function must at least set the value of initStates
    virtual void initializeTest(void) = 0;
    // return the correct states at the requested time
    virtual vector<double> getCorrectStates(const double & time) const = 0;
    // return the correct states at the end of the trial
    vector<double> getCorrectStates(void) const {
      return getCorrectStates(tFinal);
    }
    
    vector<double> initStates;
    double externalForcing;
  protected:
    // compute d/dt of each state
    virtual void getDerivatives(const std::vector<FDouble> & states,
                                const double & time,
                                std::vector<FDouble> & dStates) const = 0;
    // return the integral of externalForcing
    double integralForcing(const double & time) const;
    
  private:
    // return Clamp Trace with synthetic data to be clamped into
    //   externalForcing
    Trace getForcingTrace(const double & sampleDT);
    Trace const* forcingTrace;
};


double
TestDiffEq::getError(const double & accGoal)
{
  // return the true error in the final states at requested accuracy goal
  accuracyGoal = accGoal;
  solve(initStates);
  
  const vector<double> finalStates = getFinalStates();
  const vector<double> correctStates = getCorrectStates();
  
  
  auto cItr = correctStates.cbegin();
  auto fItr = finalStates.cbegin();
  double maxErr = abs(*cItr - *fItr);
  for(; fItr != finalStates.cend(); ++fItr, ++cItr) {
    maxErr = std::max(maxErr, abs(*cItr - *fItr));
  }
  return maxErr;
}


double
TestDiffEq::getTraceError(const double & accGoal)
{
  // return the true error in the recorded traces at requested accuracy goal
  accuracyGoal = accGoal;
  solve(initStates);
  
  size_t numTraces = 0;
  double error = 0;
  for(const Trace & trace : traces) {
    if(trace.traceAction != Trace::Fit) {
      // this isn't a fit trace, ignore it
      continue;
    }
    ++numTraces;
    error = std::max(error, trace.getMaxError());
  }
  if(numTraces == 0) {
    return NaN;
  }
  else {
    return error;
  }
}



Trace
TestDiffEq::getForcingTrace(const double & sampleDT)
{
  // return Clamp Trace with synthetic data to be clamped into externalForcing
  Trace trace;
  for(double time = tInitial; time <= tFinal; time += sampleDT) {
    trace.trace.push_back(1.0 + sin(Pi * time));
  }
  trace.setTargetName("externalForcing");
  trace.traceAction = Trace::Clamp;
  trace.deltaT = sampleDT;
  trace.waitTime = 0.0;
  return trace;
}



void
TestDiffEq::test(double sampleDT)
{
  // run comprehensive test

  verbosity = 0; // set verbosity to zero so the test output isn't overwhelming
  // initialize test and ensure that initial state is okay
  initializeTest();
  const size_t numStates = initStates.size();
  if(numStates == 0) {
    throw runtime_error("TestDiffEq::initializeTest() did not set initStates");
  }
  
  // let traces target the externalForcing variable, which has no units
  externalForcing = 0.0;
  addTraceTarget(externalForcing, "externalForcing", "none");
  
  // if sampleDT == 0.0, don't create traces, otherwise create 1 for each state
  //   and 1 more for clamping external forcing
  list<Trace> initTraces (sampleDT <= 0.0 ? 0 : numStates + 1);
  
  if(sampleDT > 0.0) {
    // create forcing trace
    initTraces.back() = getForcingTrace(sampleDT);
    forcingTrace = &initTraces.back();
    
    // create a Fit trace for each state. These will compute error
    
    // fill out each trace with the correct values at each time point
    for(double time = tInitial; time <= tFinal; time += sampleDT) {
      const vector<double> correct = getCorrectStates(time);
      auto cItr = correct.cbegin();
      for(Trace & trace : initTraces) {
        trace.trace.push_back(*cItr);
        ++cItr;
      }
    }
    // resize the error trace for each trace, and set the trace as Fit type
    size_t n = 0;
    for(Trace & trace : initTraces) {
      trace.setTargetName(numToString(n));
      trace.traceAction = Trace::Fit;
      trace.errorTrace.resize(trace.trace.size());
      trace.deltaT = sampleDT;
      trace.fitTau = Inf;
      trace.waitTime = 0.0;
      trace.filterTau = 0.0;
      ++n;
      if(n == numStates) {
        break;
      }
    }
  }
  
  // initialize the DiffEqFunctor object
  initializeDiffEq(initTraces, numStates);
  
  // set pointer to forcingTrace
  forcingTrace = NULL;
  for(const Trace & trace : traces) {
    if(trace.traceAction == Trace::Clamp &&
       trace.targetName == "externalForcing") {
      forcingTrace = &trace;
      break;
    }
  }
  
  vector<double> accGoals =
    {1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10, 1.0e-11,
     1.0e-12};
  
  vector<double> logGoals (accGoals.size());
  vector<double> logErrors (accGoals.size());
  vector<double> logNumEvals (accGoals.size());
  vector<double> logNumMatrixExp (accGoals.size());
  for(size_t n = 0; n < accGoals.size(); ++n) {
    const double error_n = traces.empty() ?
                           getError(accGoals[n]) :
                           getTraceError(accGoals[n]);
    
    const size_t numEvals_n = getNumFunctionEvaluations();
    const size_t numMatrixExp_n = getNumMatrixExponentials();
    cout << "  acc Goal=" << accGoals[n] << " : err=" << error_n
         << ", numEvals=" << numEvals_n << ", numMatExp=" << numMatrixExp_n
         << '\n';
    cout.flush();
    
    logGoals[n] = log10(accGoals[n]);
    logErrors[n] = log10(error_n);
    logNumEvals[n] = log10((double)numEvals_n);
    logNumMatrixExp[n] = log10((double)numMatrixExp_n);
  }
  
  fitErrors(logNumEvals, logErrors, "numEvals");
  fitErrors(logNumEvals, logErrors, "numMatrixExp");
}



double
TestDiffEq::integralForcing(const double & time) const
{
  // return the integral of externalForcing
  if(forcingTrace == NULL || time <= tInitial) {
    return 0.0;
  }
  
  double integral = 0;
  double tRemain = time - tInitial;
  const double & dT = forcingTrace->deltaT;
  
  auto forceItr = forcingTrace->trace.cbegin();
  while(tRemain > dT) {
    integral += *forceItr * dT;
    ++forceItr;
    tRemain -= dT;
  }
  integral += *forceItr * tRemain;
  return integral;
}



 //////////////////////////// Linear TestDiffEq Class /////////////////////////

class LinearTestDiffEq : public TestDiffEq
{
  public:
    virtual void initializeTest(void) {
      // this function must at least set the value of initStates
      gain = 2 * Pi;
      tInitial = 0.0;
      tFinal = 10.25;
      initStates = {1.0, 0.0};
    }
    
    virtual vector<double> getCorrectStates(const double & time) const {
      // return the correct states at the end of the trial
      const double dT = gain * (time - tInitial);
      const double cosT = cos(dT);
      const double sinT = sin(dT);
      return {
        cosT * initStates[0] - sinT * initStates[1],
        sinT * initStates[0] + cosT * initStates[1]
      };
    }

    double gain;
    
  protected:
    virtual void getDerivatives(const vector<FDouble> & states,
                                const double & time,
                                vector<FDouble> & dStates) const {
      // compute d/dt of each state
      (void)time;
      dStates[0] = -gain * states[1];
      dStates[1] = gain * states[0];
    }
};



/////////////////////////// Nonlinear TestDiffEq Class ///////////////////////

class NonlinearTestDiffEq : public TestDiffEq
{
  public:
    virtual void initializeTest(void) {
      // this function must at least set the value of initStates
      gainLogistic = 1.0;
      gainDiff = 10.0;
      tInitial = 0.0;
      tFinal = 7.0;
      initStates = {10.0, -9.99};
    }
  
    virtual vector<double> getCorrectStates(const double & time) const {
      // return the correct states at the end of the trial
      const double dT = time - tInitial;
      const double u0 = 0.5 * (initStates[0] + initStates[1]);
      const double d0 = 0.5 * (initStates[0] - initStates[1]);
      const double u = u0 * exp(gainLogistic * dT)
                       / (1.0 - u0 + u0 * exp(gainLogistic * dT));
      const double d = d0 * exp(-gainDiff * dT - integralForcing(time));
      return {u + d, u - d};
    }


    double gainLogistic;
    double gainDiff;

  protected:
    virtual void getDerivatives(const vector<FDouble> & states,
                                const double & time,
                                vector<FDouble> & dStates) const {
      // compute d/dt of each state
      (void)time;
      const double gD = gainDiff + externalForcing;
      dStates[0] = 0.5 * (gainLogistic * (states[0] + states[1]) *
                                         (1.0 - 0.5 * (states[0] + states[1]))
                         -gD * (states[0] - states[1]));
      dStates[1] = 0.5 * (gainLogistic * (states[0] + states[1]) *
                                         (1.0 - 0.5 * (states[0] + states[1]))
                         +gD * (states[0] - states[1]));
    }
};



///////////////////////////// Stiff TestDiffEq Class //////////////////////////

class StiffTestDiffEq : public TestDiffEq
{
  public:
    virtual void initializeTest(void) {
      // this function must at least set the value of initStates
      
      // seems to work for positive integer powers, bombs otherwise
      tPower = 0.5;
      gainCross = 0.1;
      tInitial = 0.1;
      tFinal = tInitial + 10.0;
      const double ratio = 1000.0;  // seemingly a measure of stiffness
      
      // complete 1 full UV circle
      gainUV = 2.0 * Pi / (tFinal - tInitial);
      // complete ratio full XY circles
      const double intPow = 1.0 + tPower;
      gainXY = gainUV * ratio * intPow
             / (pow(tFinal, intPow) - pow(tInitial, intPow));

      initStates = {1.0, 0.0, 1.0, 0.0};
    }
  
    virtual vector<double> getCorrectStates(const double & time) const {
      // return the correct states at the end of the trial
      
      const double intPow = 1.0 + tPower;
      const double F = (gainXY / intPow)
                     * (pow(time, intPow) - pow(tInitial, intPow));
      const double cosF = cos(F);
      const double sinF = sin(F);
      const double x = initStates[0] * cosF - initStates[1] * sinF;
      const double y = initStates[0] * sinF + initStates[1] * cosF;
      const double G = gainUV * (time - tInitial)
                     + gainCross * (y - initStates[1]);
      const double cosG = cos(G);
      const double sinG = sin(G);
      const double u = initStates[2] * cosG - initStates[3] * sinG;
      const double v = initStates[2] * sinG + initStates[3] * cosG;
      
      return {x, y, u, v};
    }


    double gainXY;
    double tPower;
    double gainUV;
    double gainCross;

  protected:
    virtual void getDerivatives(const vector<FDouble> & states,
                                const double & time,
                                vector<FDouble> & dStates) const {
      // compute d/dt of each state
      const double f = tPower == 0.0 ? gainXY : gainXY * pow(time, tPower);
      dStates[0] = -f * states[1];
      dStates[1] = f * states[0];
      const auto g = gainUV + gainCross * dStates[1];
      dStates[2] = -g * states[3];
      dStates[3] = g * states[2];
    }
};



///////////////////////////////////// main ////////////////////////////////////

int main(void)
{
  Timer timer;
  cout << "Testing linear differential equation (no traces, final state error):\n";
  LinearTestDiffEq diffLinear;
  diffLinear.test();
  cout << "Testing linear differential equation (with trace error):\n";
  diffLinear.test(0.1);
  
  cout << "Testing non-linear differential equation (no traces, final state error):\n";
  NonlinearTestDiffEq diffNonlinear;
  diffNonlinear.test();
  cout << "Testing non-linear differential equation: (with trace error)\n";
  diffNonlinear.test(0.1);
  
  cout << "Testing stiff non-linear differential equation (no traces, final state error):\n";
  StiffTestDiffEq diffStiff;
  diffStiff.test();
  cout << "Testing stiff non-linear differential equation (with trace error):\n";
  diffStiff.test(0.1);
  
  cout << "\nTotal elapsed time: " << timer << '\n';
  return 0;
}
