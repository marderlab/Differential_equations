#include "../include/Timer.h"
#include "../include/io_functions.h"



using namespace std;



///////////////////////// Timer public member functions ///////////////////////
Timer::Timer()
{
  gettimeofday(&startTime, NULL);
}



Timer::Timer(timeval & startT)
{
  startTime = startT;
}



void
Timer::start(void)
{
  gettimeofday(&startTime, NULL);
}



double
Timer::elapsedTime(void) const
{
  gettimeofday(&now, NULL);
  return (double)(now.tv_sec - startTime.tv_sec) +
                   (double)(now.tv_usec - startTime.tv_usec) / 1000000.0;
}



timeval
Timer::currentTime(void) const
{
  gettimeofday(&now, NULL);
  return now;
}



void
Timer::addElapsedTime(const double & addTime)
{
  double currentStart = (double)startTime.tv_sec +
                        1.0e-6 * (double)startTime.tv_usec;
  currentStart -= addTime;

  startTime.tv_sec = (__time_t)currentStart;
  startTime.tv_usec = (__suseconds_t)(1.0e6
			      * (currentStart - (double)startTime.tv_sec));
}



//////////////////////////// Timer friend functions ///////////////////////////
ostream &
operator<<(ostream & out, const Timer & timer)
{
  // output elapsed time to stream, as human readable string
  
  double seconds = timer.elapsedTime();
  
  size_t minutes = (size_t) (seconds / 60);
  seconds -= (double)(minutes * 60);
  
  size_t hours = minutes / 60;
  minutes -= hours * 60;
  
  size_t days = hours / 24;
  hours -= days * 24;

  // record format flags
  const ios_base::fmtflags formatFlags = out.flags();
  // record width;
  const streamsize width = out.width();
  // record precision
  const streamsize precision = out.precision();
  
  if(days > 0)
    out << left << days << "d ";
  if(hours > 0)
    out << left << hours << "h ";
  if(minutes > 0)
    out << left << minutes << "m ";

  out.precision(6);
  out.setf(ios::fixed);
  streamsize secWidth = seconds >= 10.0 ? 8 : 7;
  out.width(secWidth);
  out << left << seconds << 's';

  // restore flags
  out.setf(formatFlags);
  // restore width
  out.width(width);
  // restore precision
  out.precision(precision);
  return out;
}



istream &
operator>>(istream & in, Timer & timer)
{
  // input elapsed time from stram, as human readable string
  
  double addTime = 0;
  string word;
  
  bool reading = true;

  while(reading) {
    //read in the next word
    in >> word;
    
    double wordTime;
    switch(*(word.rbegin())) {
      // adjust time value by the unit appended
      case 'd':
        wordTime = 86400.0;
        break;
      case 'h':
        wordTime = 3600.0;
        break;
      case 'm':
        wordTime = 60.0;
        break;
      case 's':
        reading = false;
        wordTime = 1.0;
        break;
      default:
        throw TimerException("Error reading Timer from input stream");
    }
      
    // remove the unit from the word
    word.erase(word.end() - 1);

    // multiply the value of the word into wordTime
    wordTime *= stringToDouble(word);
    
    // add the value of the word to addTime
    addTime += wordTime;
  }
  
  // add addTime to elapsed time
  timer.addElapsedTime(addTime);
  
  return in;
}
