#ifndef _IO_FUNCTIONS_CPP_
#define _IO_FUNCTIONS_CPP_



#include "../include/io_functions.h"
#include "../include/constants.h"

#include <algorithm>



using namespace std;



bool
fileExists(const char* fileName)
{
  return ifstream(fileName) != NULL;
}



bool
fileExists(const string & fileName)
{
  return ifstream(fileName.c_str()) != NULL;
}



void
userExpand(string & path)
{
  // search path string for '~' and replace it with the user's home directory
  
  size_t tildePos = path.find('~');
  if(tildePos != std::string::npos) {
    // found at least one tilde (presumably the only one, but we'll check)
    
    // get the home directory
    const std::string homeDir = getenv("HOME");
    // replace the tilde
    path.replace(tildePos, 1, homeDir);
    // search for the next one
    tildePos = path.find('~', tildePos + homeDir.size());
    
    while(tildePos != std::string::npos) {
      // while we find new tildes, replace and search for next
      path.replace(tildePos, 1, homeDir);
      tildePos = path.find('~', tildePos + homeDir.size());
    }
  }
}



OverrideMap parseArgs(int numArgs, char** arguments, ArgMap & argMap)
{
  // parse input arguments, looking for overrides of default arguments
  
  // create a map to determine if arguments are override
  OverrideMap override;
  ArgMap::iterator argMapItr;
  for(argMapItr = argMap.begin(); argMapItr != argMap.end(); ++argMapItr)
    // initialize all to false (not override)
    override[argMapItr->first] = false;
  // add "help" to override
  override["help"] = false;
  
  int argNum = 1;
  while(argNum < numArgs) {
    std::string arg = arguments[argNum];
    if(arg == "help" || arg == "-help" || arg == "--help") {
      override["help"] = true;
      ++argNum;
    }
    else if(arg[0] == '-') {
      // keyword, next argument should be a value
      // remove the dash
      arg.erase(0, 1);
      // inc argNum
      ++argNum;
      if(argNum == numArgs)
        throw IOException("Argument passed with keyword but no value");
      // set the next value
      argMap[arg] = arguments[argNum];
      // mark as override
      override[arg] = true;
      // inc argNum
      ++argNum;
    }
    else {
      // override the next non-overriden argument
      OverrideMap::iterator overrideItr = override.begin();
      for(; overrideItr != override.end(); ++overrideItr)
        if(!overrideItr->second) {
          argMap[overrideItr->first] = arg;
          overrideItr->second = true;
          ++argNum;
          break;
        }
      if(overrideItr == override.end())
        throw IOException("Too many arguments");
    }
  }
  return override;
}



std::string
userExpand(const string & constStr)
{
  // search string for '~' and replace it with the user's home directory
  
  string path = constStr;
  userExpand(path);
  return path;
}



void
removeWhiteSpace(string & str)
{
  // remove all the white space from a string
  str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());
}



string
strip(const string & str)
{
  // strip leading and trailing white space from a string
  
  // find an iterator to the first non-white-space
  string::const_iterator start;
  string::const_iterator stop = str.end();
  for(start = str.begin(); start != stop; ++start)
    if(!::isspace(*start))
      break;

  if(start == stop)
    // we've reached the end, return an empty string
    return string();
  
  // find an iterator past the last non-white-space
  for(--stop; (stop != start) && (::isspace(*stop)); --stop)
    ;
  // inc stop
  ++stop;
  
  // return a new string object
  return string(start, stop);
}



string
lstrip(const string & str)
{
  // strip leading white space from a string
  
  // find an iterator to the first non-white-space
  string::const_iterator start;
  string::const_iterator stop = str.end();
  for(start = str.begin(); start != stop; ++start)
    if(!::isspace(*start))
      break;

  if(start == stop)
    // we've reached the end, return an empty string
    return string();
  
  // return a new string object
  return string(start, stop);
}



string
rstrip(const string & str)
{
  // strip trailing white space from a string
  
  // find an iterator past the last non-white-space
  string::const_iterator start = str.begin();
  string::const_iterator stop = str.end();
  for(--stop; (stop != start) && (::isspace(*stop)); --stop)
    ;

  if(start == stop && ::isspace(*stop))
    // we've reached the beginning, return an empty string
    return string();
  // inc stop
  ++stop;
  
  // return a new string object
  return string(start, stop);
}



std::string
getNextLine(istream & in, size_t & lineNum, const char & commentChar)
{
  // get the next non-empty line from a stream, removing comments
  
  while(in.good()) {
    // loop until done
    
    // inc lineNum
    ++lineNum;
    
    // get the next line
    string line;
    getline(in, line);
    in.peek();  //peek helps trip eof flag
    
    if(commentChar != '\0') {
      // look for comment, remove if present
      size_t commentInd = line.find(commentChar);
      if(commentInd != string::npos) {
        // remove everything after comment
        line.erase(commentInd);
      }
    }
    
    // remove leading and trailing white space
    line = strip(line);
    
    if(!line.empty())
      // line is non-empty, return it
      return line;
  }
  // It's probably better to return an empty string than to throw an error
  // throw IOException("Reached end of stream without finding next line");
  return string ();
}



SplitLine
split(const string & line, const string & delimiter, const size_t & maxSplits)
{
  // split up a line into a list of words
  //   words are separated by delimiters in string 'delimiter'
  //     (if delimiter is empty, split by whitespace)
  //   perform maxSplits splits
  //   return a Splitline
  
  if(delimiter.empty())
    // assume user wants whitespace
    return split(line, whitespace, maxSplits);
  
  if(line.empty())
    // return empty list for empty line
    return SplitLine (1, string());
  
  size_t numSplits = 0;
  SplitLine splitLine;
    
  // find the first char that IS NOT a delimiter
  size_t ind1 = line.find_first_not_of(delimiter, 0);
  
  while(ind1 != string::npos) {
    // loop until the last word is found
    
    size_t ind2;
    if(numSplits == maxSplits) {
      // put the rest as the last word of splitLine
      
      // find index to last non-delimiter char
      ind2 = line.find_last_not_of(delimiter);
      // point to just after last non-delimiter
      ++ind2;
      // resize splitLine
      splitLine.resize(numSplits + 1);
      // put the word into the last string in splitLine
      splitLine.back() = line.substr(ind1, ind2 - ind1);
      break;
    }
    else
      // find the first char after ind1 that IS a delimiter
      ind2 = line.find_first_of(delimiter, ind1 + 1);
    
    // inc numSplits
    ++numSplits;
    // resize splitLine
    splitLine.resize(numSplits);
    // put the word into the last string in splitLine
    if(ind2 == string::npos) {
      splitLine.back() = line.substr(ind1, string::npos);
      ind1 = string::npos;  // got the last word
    }
    else {
      splitLine.back() = line.substr(ind1, ind2 - ind1);
      // advance to next char after ind2 that IS NOT a delimiter
      ind1 = line.find_first_not_of(delimiter, ind2 + 1);
    }
  }

  return splitLine;
}



inline string
reverse(const string & str)
{
  return string (str.rbegin(), str.rend());
}

inline SplitLine
reverse(const SplitLine & splitLine)
{
  SplitLine flipped;
  flipped.resize(splitLine.size());
  
  size_t m = splitLine.size() - 1;
  for(size_t n = 0; n < splitLine.size(); ++n, --m)
    flipped[n] = reverse(splitLine[m]);
  
  return flipped;
}

SplitLine
rsplit(const std::string & line, const std::string & delimiter,
       const size_t & maxSplits)
{
  return reverse(split(reverse(line), delimiter, maxSplits));
}



double
stringToDouble(const string & word)
{
  return atof(word.c_str());
}



long
stringToInt(const string & word)
{
  return atol(word.c_str());
}



size_t
stringToSize(const string & word)
{
  return (size_t) strtoul(word.c_str(), NULL, 0);
}



double
getDouble(istream & in)
{
  string word;
  in >> word;
  return atof(word.c_str());
}



long
getInt(istream & in)
{
  string word;
  in >> word;
  return atol(word.c_str());
}



size_t
getSize(istream & in)
{
  string word;
  in >> word;
  return (size_t) strtoul(word.c_str(), NULL, 0);
}



bool
getNumber(istream & fileIn, double & x, string* priorStr)
{
  int numFound;
  if(priorStr != NULL)
    *priorStr = "";
  string temp;
  do{
    fileIn >> temp;
    numFound = sscanf(temp.c_str(), "%lg", &x);
    if(priorStr != NULL && numFound == 0)
      *priorStr = temp;
  } while(numFound == 0 && !fileIn.fail() && !fileIn.eof());
  
  if(numFound == 1)
    return true;
  else
    return false;
}



bool
getNumber(istream & fileIn, unsigned long & x, string* priorStr)
{
  int numFound;
  if(priorStr != NULL)
    *priorStr = "";
  string temp;
  do{
    fileIn >> temp;
    numFound = sscanf(temp.c_str(), "%lu", &x);
    if(priorStr != NULL && numFound == 0)
      *priorStr = temp;
  } while(numFound == 0 && !fileIn.fail() && !fileIn.eof());

  if(numFound == 1)
    return true;
  else
    return false;
}
#endif
