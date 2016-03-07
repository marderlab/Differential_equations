#include "../include/Trace.h"
#include <cmath>
#include <iostream>


using namespace std;



//////////////////////////////  Helper functions  /////////////////////////////



inline double
square(const double & x)
{
  return x * x;
}

///////////////////////  Trace public member functions  ///////////////////////



void
Trace::setTargetName(const string & targetStr)
{
  targetName = targetStr;
  targetDescriptiveName = targetName;
}



void
Trace::setTarget(TraceTarget traceTarget)
{
  targetPtr = traceTarget.targetPtr;
  units = traceTarget.units;
}



bool
Trace::matches(const string & targetStr,
               const TraceAction & traceAct,
               const double & dT) const
{
  // check if the information provided matches this trace
  
  if(targetStr != targetName || traceAct != traceAction)
    return false;
  if(traceAction == Record && deltaT != dT)
    return false;
  
  return true;
}



void
Trace::initialize(void)
{
  // prepare trace to begin simulation/fit or read-in data
  //   note: targetPtr can be NULL if data is going to be read in

  traceItr = trace.begin();
  switch(traceAction) {
    case Record: {
      nextT = 0.0;
      break;
    }
    case Clamp: {
      nextT = deltaT;
      break;
    }
    case Fit: {
      nextT = 0.0;
      errorItr = errorTrace.begin();
      errorFact = 1.0 - exp(-deltaT / fitTau);
      filterTau = 0.0;
      countdownT = waitTime;
      
      if(filterTau > 0 && filter.empty()) {
        const size_t midLen = (size_t)(filterTau / deltaT);
        const size_t filterLen = 1 + 2 * midLen;
        filter.resize(filterLen);
        filter[midLen] = 1.0;
        size_t indPlus = midLen + 1;
        size_t indMinus = midLen - 1;
        const double dTFilter = deltaT / filterTau;
        double x = dTFilter;
        double total = 1.0;
        for(size_t n = 0; n < midLen; ++n) {
          filter[indPlus] = exp(-x*x);
          filter[indMinus] = filter[indPlus];
          total += 2.0 * filter[indPlus];
          ++indPlus;
          --indMinus;
          x -= dTFilter;
        }
        for(auto & filtVal : filter)
          filtVal /= total;
      }
      break;
      
    }
  }
}



void
Trace::readInValue(const double & value)
{
  *traceItr = value;
  ++traceItr;
}



void
Trace::doTraceAction(bool & alterModel)
{
  // carry out traceAction
  
  // nextT stores time until next trace action, so the trace will iterate
  // forward while nextT < 0
  
  switch(traceAction) {
    case Record: {
      if(nextT == 0.0) {
        // record target
        *traceItr = *targetPtr;
        // increment counter by the trace recording period minus dT
        nextT = deltaT;
        // advance trace
        ++traceItr;
      }
      break;
    }
    case Clamp: {
      if(nextT == 0.0) {
        // signal that the model will be altered
        alterModel = true;
        // increment counter by the trace recording period minus dT
        nextT = deltaT;
        // advance trace
        ++traceItr;
      }

      // clamp target to trace value
      *targetPtr = *traceItr;
      break;
    }
    
    case Fit: {
      if(nextT == 0.0) {
        // signal that the model will be altered
        alterModel = true;
        // error is the difference between target value and trace value
        *errorItr = *targetPtr - *traceItr;
        // partially correct the error
        if(countdownT <= 0.0) {
          *targetPtr -= (*errorItr * errorFact);
        }
        else {
          countdownT -= deltaT;
        }
        // increment counter by the trace recording period minus dT
        nextT = deltaT;
        // advance traces
        ++traceItr;
        ++errorItr;
      }
    }
  }
}



void
Trace::advance(const double & dT)
{
  // advance time, but do not perform trace action

  nextT -= dT;
}



double
Trace::getRMSError(void) const
{
  // return the L2 error of the trace
  if(traceAction != Fit) {
    return 0.0;
  }

  if(filter.empty()) {
    double error = 0.0;
    for(const auto & err : errorTrace)
      error += err * err;
    return sqrt(error / (double)(errorTrace.size()));
  }
  else {
    double error = 0.0;
    auto eStartItr = errorTrace.cbegin();
    auto eStopItr = eStartItr + filter.size();
    while(eStopItr != errorTrace.cend()) {
      double filtErr = 0.0;
      auto filtItr = filter.cbegin();
      for(auto eItr = eStartItr; eItr != eStopItr; ++eItr) {
        filtErr += *eItr * *filtItr;
      }
      error += filtErr * filtErr;
      ++eStartItr;
      ++eStopItr;
    }
    return sqrt(error / (double)(errorTrace.size() - filter.size()));
  }
}


double
Trace::getMaxError(void) const
{
  // return the max square error of the trace
  if(traceAction != Fit) {
    return 0.0;
  }

  if(filter.empty()) {
    double error = 0.0;
    for(const auto & err : errorTrace)
      if(abs(err) > error)
        error = abs(err);
    return error;
  }
  else {
    double error = 0.0;
    auto eStartItr = errorTrace.cbegin();
    auto eStopItr = eStartItr + filter.size();
    while(eStopItr != errorTrace.cend()) {
      double filtErr = 0.0;
      auto filtItr = filter.cbegin();
      for(auto eItr = eStartItr; eItr != eStopItr; ++eItr) {
        filtErr += *eItr * *filtItr;
      }
      if(abs(filtErr) > error)
        error = abs(filtErr);
      ++eStartItr;
      ++eStopItr;
    }
    return error;
  }
}

double
Trace::getL4Error(void) const
{
  // return the L4 (root mean 4th-power) error of the trace
  if(traceAction != Fit) {
    return 0.0;
  }

  if(filter.empty()) {
    double error = 0.0;
    for(const auto & err : errorTrace)
      error += square(err * err);
    return pow(error / (double)(errorTrace.size()), 0.25);
  }
  else {
    double error = 0.0;
    auto eStartItr = errorTrace.cbegin();
    auto eStopItr = eStartItr + filter.size();
    while(eStopItr != errorTrace.cend()) {
      double filtErr = 0.0;
      auto filtItr = filter.cbegin();
      for(auto eItr = eStartItr; eItr != eStopItr; ++eItr) {
        filtErr += *eItr * *filtItr;
      }
      error += square(filtErr * filtErr);
      ++eStartItr;
      ++eStopItr;
    }
    return pow(error / (double)(errorTrace.size() - filter.size()),
                0.25);
  }
}


double
Trace::getL1Error(void) const
{
  // return the L1 (mean) error of the trace
  if(traceAction != Fit) {
    return 0.0;
  }
  
  if(filter.empty()) {
    double error = 0.0;
    for(const auto & err : errorTrace)
      error += abs(err);
    return (error / (double)(errorTrace.size()));
  }
  else {
    double error = 0.0;
    auto eStartItr = errorTrace.cbegin();
    auto eStopItr = eStartItr + filter.size();
    while(eStopItr != errorTrace.cend()) {
      double filtErr = 0.0;
      auto filtItr = filter.cbegin();
      for(auto eItr = eStartItr; eItr != eStopItr; ++eItr) {
        filtErr += *eItr * *filtItr;
      }
      error += abs(filtErr);
      ++eStartItr;
      ++eStopItr;
    }
    return (error / (double)(errorTrace.size() - filter.size()));
  }
}

