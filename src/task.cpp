// -*- C++ -*-
// $Id$

#include "task.h"
#include <algorithm>


/** @brief Register the parameters needed for the Task
 *
 *  A derived class should overwrite this whenever it needs additional 
 *  parameters from the setup file. It should then ALWAYS call the 
 *  MakeParamMap of its superclass.
 */
PARAMETERMAP* Task::MakeParamMap (PARAMETERMAP* pm) {
  if (NULL == pm) pm = new PARAMETERMAP;
  return pm;
}

std::string Task::Rebuild (std::istream& in) {
  PARAMETERMAP* pm = MakeParamMap();
  std::string strToken;
  in >> strToken;
  while ((strToken != "}") && (!in.eof())) {
    while (strToken == "//") { // Kommentar
      char ch;
      
      do {
		    in.get(ch);
      } while ((ch != '\n') && (!in.eof()));
      
      if (!in.eof()) in >> strToken;
    }
    if ("}" != strToken) {
      if (0 == pm->count(strToken)) { // Token nicht vorgesehen
        std::string strParam = strToken;
				in >> strToken;
				if (strToken == "{") { // ist ein Task, sollte also uebersprungen werden
				  std::cerr << "Unknown Task " << strParam << ". Skipping... " << std::endl;
				  int nLevel = 0;
				  do {
				    in >> strToken;
				    if (strToken == "{") nLevel++;
				    if (strToken == "}") nLevel--;
				  } while (((strToken != "}") || (nLevel >= 0)) && (!in.eof()));
				  
          if (in.eof()) 
				    std::cerr << "Unexpected end of file." << std::endl;
				  else
				    in >> strToken;
				}
				else // nur den Namen des Parameters ausgeben, es koennte ein richtiger kommen
				  std::cerr << "Unknown Parameter " << strParam << std::endl;
      }
      else
        strToken = (((Parameter*)(*pm)[strToken])->Rebuild(in));
    }
  }
  if (strToken == "}") in >> strToken;
  while (strToken == "//") { // Kommentar
    char ch;
    do {
      in.get(ch);
    } while ((ch != '\n') && (!in.eof()));
    if (!in.eof()) 
      in >> strToken;
    else
      break;
  }
  
  // hier muss noch pm aufgeraeumt werden
  while (!pm->empty()) {
    PARAMETERMAP::iterator iter = pm->begin();
    delete ((Parameter*)(*iter).second);
    pm->erase(iter);
  }
  delete pm;
  return strToken;
}



