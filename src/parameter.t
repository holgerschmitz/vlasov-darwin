// -*- C++ -*-
// $Id: parameter.t,v 1.1.1.1 2006/05/11 10:02:00 th Exp $

//the rebuild method template. called by Rebuildable objects
//gets *pValue from stream returns next token
template<class Type>
std::string ParameterValue<Type>::Rebuild (std::istream& in) {
	  in >> *pValue;
	  std::string strToken;
	  in >> strToken;
	return strToken;
}


template<class Type, class BaseType>
std::string ParameterRebuild<Type,BaseType>::Rebuild (std::istream& in) {
	  // in die Liste der Subtasks des Elterntasks eintragen
	  (*value) = NewInstance();
  
	  // Task aus der Datei wiederherstellen
	  std::string strToken;
	  in >> strToken;
	  if ((strToken != "{") || in.eof()) return strToken;
	  std::string nextToken = (*value)->Rebuild(in);
	  (*value)->finalize();
  return nextToken;
}
