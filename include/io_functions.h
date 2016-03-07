#ifndef _IO_FUNCTIONS_H_
#define _IO_FUNCTIONS_H_



#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <map>



bool fileExists(const char* fileName);
bool fileExists(const std::string & fileName);



// create a typedef to hold a map of argument names to values
typedef std::map<std::string, std::string> ArgMap;

// create a typedef to hold a map of argument names to boolean (overridden or
// not)
typedef std::map<std::string, bool> OverrideMap;

// parse input arguments, looking for overrides of default arguments
// the override map always contains the key "help", and will be set to true
// if one of the arguments is "-help"
OverrideMap parseArgs(int numArgs, char** arguments, ArgMap & argMap);



// search string for '~' and replace it with the user's home directory
void userExpand(std::string & str);
std::string userExpand(const std::string & constStr);



// remove all the white space from a string
void
removeWhiteSpace(std::string & str);
// strip leading and trailing white space from a string
std::string
strip(const std::string & str);
// (left) strip leading white space from a string
std::string
lstrip(const std::string & str);
// (right) strip trailing white space from a string
std::string
rstrip(const std::string & str);



// get the next non-empty line from a stream, removing comments
std::string
getNextLine(std::istream & in, size_t & lineNum,
            const char & commentChar = '\0');



// define a type for lines that have been split into words
// (currently implemented as vector of strings)
typedef std::vector<std::string> SplitLine;

// define a string that contains white space delimiters
const std::string whitespace (" \t\n\v\f\r");

// split up a line into a list of words
//   words are separated by delimiters in string 'delimiter'
//     (default = whitespace)
//   perform maxSplits splits
//   return a Splitline
SplitLine
split(const std::string & line, const std::string & delimiter = whitespace,
      const size_t & maxSplits =(size_t)-1);
// (right) split up a line into a list of words, starting from the end
//   words are separated by delimiters in string 'delimiter'
//     (default = whitespace)
//   perform maxSplits splits
//   return a Splitline
SplitLine
rsplit(const std::string & line, const std::string & delimiter = whitespace,
       const size_t & maxSplits =(size_t)-1);



// convert string to numbers
double
stringToDouble(const std::string & word);
long
stringToInt(const std::string & word);
size_t
stringToSize(const std::string & word);
/*
// problems with stringToNum()
//    -reading nan, inf, -inf
//    -slow
template <class T>
void stringToNum(const std::string & word, T & num);
*/



// convert number to string
template <class T>
std::string
numToString(const T & num);



// get numbers from stream
double getDouble(std::istream & in);
long getInt(std::istream & in);
size_t getSize(std::istream & in);



//########### Template function definitions #####################
template <class T>
std::string numToString(const T & num)
{
  std::stringstream ss;
  ss << num;
  return ss.str();
}


/*
template <class T>
void stringToNum(const std::string & word, T & num)
{
  std::stringstream streamWord (word);
  streamWord >> num;
}
*/



class IOException : public std::runtime_error
{
  public:
    IOException() :
      std::runtime_error("IOException") {}
    IOException(const std::string & str) :
      std::runtime_error(str) {}
 };


#endif
