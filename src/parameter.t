// -*- C++ -*-
// $Id:$

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
  return (*value)->Rebuild(in);
}
