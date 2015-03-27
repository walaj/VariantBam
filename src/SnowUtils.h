#ifndef SNOWUTILS_H
#define SNOWUTILS_H

#include <string>
#include <time.h>
#include <ctime>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>

namespace SnowUtils {

template <typename T> 
std::string AddCommas(T data) { 
  std::stringstream ss; ss << data; std::string s = ss.str();
  if (s.length() > 3)
     for (int i = s.length()-3; i > 0; i -= 3)
       s.insert(i,",");
   return s;
}

inline bool exist_test (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

inline bool read_access_test (const std::string& name) {
  return (access (name.c_str(), R_OK) == 0); 
}



inline void displayRuntime(const timespec start) {

  struct timespec finish;
  clock_gettime(CLOCK_MONOTONIC, &finish);
  double elapsed = (finish.tv_sec - start.tv_sec);
  int t = clock()/CLOCKS_PER_SEC;
  int min = (int)floor(elapsed / 60.0);
  int sec = (int)(elapsed-min*60);
  char buffer[100];
  sprintf (buffer, "CPU: %4dm%02ds Wall: %4dm%02ds", 
            (int)floor( ((double)t) /60.0), t % 60, min, sec);
  printf ("%s",buffer);
}

inline void rcomplement(std::string &a) {

  std::reverse(&a[0], &a[a.size()]);
  std::string::iterator it = a.begin();
  for (; it != a.end(); it++)
    if (*it == 'A')
      *it = 'T';
    else if (*it == 'T')
      *it = 'A';
    else if (*it == 'C')
      *it = 'G';
    else
      *it = 'C';
}

  // calculate the percentage
 template <typename T> inline int percentCalc(T numer, T denom) {
   if (denom <= 0)
     return 0;
   int perc  = static_cast<int>(floor((float)numer / (float)denom * 100.0));
   return perc;
 }

 // remove the last character from a string
 inline std::string cutLastChar(std::string in) {
   if (in.length() == 0)
     return in;
   else 
     return in.substr(0, in.length() - 1);
 }

 // remove substrings from a string
 inline std::string scrubString(std::string toscrub, std::string toremove) {
   std::string::size_type i = toscrub.find(toremove);
   while (i != std::string::npos) {
     toscrub.erase(i, toremove.length());
     i = toscrub.find(toremove);
   }
   return toscrub;
 }


}

#endif
