#ifndef _TIMER_H_
#define _TIMER_H_



#include <sys/time.h>
#include <iostream>
#include <stdexcept>



class Timer
{
  public:
    Timer();
    Timer(timeval & startT);
    
    // start the timer (called automatically by constructor)
    void start(void);
    // return elapsed time (seconds) as a double
    double elapsedTime(void) const;
    // return current time as a timeval
    timeval currentTime(void) const;
    // for future calls to elapsed time, add constant offset (in seconds)
    // (cumulative over successive calls)
    void addElapsedTime(const double & addTime);
    // output elapsed time to stream, as human readable string
    friend std::ostream & operator<<(std::ostream & out, const Timer & timer);
    // input elapsed time from stram, as human readable string
    friend std::istream & operator>>(std::istream & in, Timer & timer);

  protected:
    mutable timeval now;
    timeval startTime;
};



class TimerException : public std::runtime_error
{
  public:
    TimerException() : std::runtime_error("TimerException") { }
    TimerException(const std::string & str) : std::runtime_error(str) { }
};



#endif
