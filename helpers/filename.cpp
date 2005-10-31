#include "filename.h"
#include <sstream>
#include <iostream>

void SimFileName::setReplace(std::string id, int val)
{
  replaceMap[id] = val;
}

std::string SimFileName::filename()
{
  std::cerr << "Creating file name from " << *this << std::endl;
  std::string parsed = *this;
   
  std::map<std::string, int>::iterator it;
  
  for (it=replaceMap.begin(); it!=replaceMap.end(); ++it)
  {
  
    std::ostringstream repstr;
    repstr << it->second;
    std::string rep = repstr.str();
    int length = it->first.length();
   int pos = parsed.find(it->first);
   while (pos != std::string::npos)
   {
     parsed.replace(pos,length,rep);
     pos = parsed.find(it->first);
   }
  }

  std::cerr << "  -- " << parsed << std::endl;
  
  return parsed;
}
