#ifndef FILENAME_H
#define FILENAME_H


#include <string>
#include <map>

class SimFileName : public std::string
{
  private:
    std::map<std::string, int> replaceMap;
  public:
    SimFileName() 
      : std::string() {}
      
    SimFileName(const char *fname) 
      : std::string(fname) {}
      
    SimFileName(const std::string &fname) 
      : std::string(fname) {}
      
    SimFileName(const SimFileName &fname) 
      : std::string(fname), replaceMap(fname.replaceMap) {}
    
    void setReplace(std::string id, int val);
    std::string filename();
};

#endif
