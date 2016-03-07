#ifndef Trace_h
#define Trace_h



#include "io_functions.h"

#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>  // for sqrt in root mean square error


//// TraceTarget class provides a convenient way to group Trace target info ///

struct TraceTarget
{
  TraceTarget() {}
  TraceTarget(double & target, const std::string & targetUnits)
    : targetPtr (&target), units (targetUnits) {}

  double* targetPtr;
  std::string units;
};



// Trace class handles recording, clamping and fitting data to the inner
//   workings of a differential equation
class Trace
{
  public:
    virtual ~Trace() {}
    enum TraceAction { Record, Clamp, Fit };
    
    Trace(void) { targetPtr = NULL; }
    // parse target string to extract components of target information
    void setTargetName(const std::string & targetStr);
    // point electrophys trace (recording/injecting/etc) at a target variable
    //   and compartment
    void setTarget(TraceTarget traceTarget);
    // check if trace has target specified by targetStr and traceAct
    bool matches(const std::string & targetStr,
                 const TraceAction & traceAct,
                 const double & dT = 0.0) const;
    // prepare trace to begin simulation/fit or read-in data
    void initialize(void);
    // read a value into the trace (typically from disk)
    void readInValue(const double & value);
    // carry out traceAction
    void doTraceAction(bool & alterModel);
    // advance time, but do not perform trace action
    void advance(const double & dT);
    // return the L2 error of the trace
    double getRMSError(void) const;
    // return the max square error of the trace
    double getMaxError(void) const;
    // return the L4 (root mean 4th-power) error of the trace
    double getL4Error(void) const;
    // return the L1 (mean) error of the trace
    double getL1Error(void) const;
    // return the time until the next trace action is taken
    inline const double & getNextT(void) const
      { return nextT == 0.0 ? deltaT : nextT; }
    // return true if the trace has iterated to its end, false otherwise
    inline bool reachedEnd(void) const { return traceItr == trace.end(); }

    // string specifying the Trace target:
    std::string targetName;
    // extra name string in case targetName is user-unfriendly. Unless
    // explicitly modified, it's set to targetName. Should only affect
    // output descriptions of the target
    std::string targetDescriptiveName;
    // string describing the units of the Trace
    std::string units;

    TraceAction traceAction;

    std::vector<double> trace;
    std::vector<double> errorTrace;
    double deltaT;
    double fitTau;
    double waitTime;
    double filterTau;
    
  private:
    double* targetPtr;
    double nextT;
    double errorFact;
    double countdownT;
    double nextScaleDeltaT;
    double shortestTraceTime;
    
    std::vector<double> filter;
    
    std::vector<double>::iterator traceItr;
    std::vector<double>::iterator errorItr;
};



//////////////// Exception class for errors related to Traces /////////////////

class TraceException : public std::runtime_error
{
  public:
    TraceException() :
      std::runtime_error("TraceException") {}
    TraceException(const std::string & str) :
      std::runtime_error(str) {}
};
 
#endif
