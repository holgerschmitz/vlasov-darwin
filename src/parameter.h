/***************************************************************************
                          parameter.h  -  description
                             -------------------
    begin                : Fri Jun 2 2000
    copyright            : (C) 2000 by Holger Schmitz
    email                : 
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/** @file parameter.h
 *  The classes in this file are used for reading parameters from
 *  an input file. The parameters are stored in a parameter map
 */

#include <map>
#include "stlpwrapper.h"

#ifndef PARAMETER_H
#define PARAMETER_H

class Task;

class Parameter;
typedef PtrWrapper<Parameter> WParameter;
typedef map<string,WParameter> PARAMETERMAP;

/** @brief An abstract base class for parameters.
 *
 *  There is only one abstract method: Rebuild. This method has to rebuild 
 *  the data from the input file and return the next token of th input file
 *  in a string.
 *  @todo Is the < operator really needed?
 */
class Parameter {
    public:
        /// Default constructor
        Parameter () {};
        // Destructor
        virtual ~Parameter () {};

        /** @brief This method has to rebuild 
         *  the data from the input file and return the next token of th input file
         *  in a string.
         */
        virtual string Rebuild (istream& in) = 0; // abstract class Parameter

        /// Always returns true since sorting is irrelevant
        friend bool operator< (const Parameter& left, const Parameter& right) {
            return true; 
        };
}; // Parameter


/** @brief Single value paramter of some type.
 *
 *  This reads in a single value of some type TYPE from the input file
 *  and stores it under a given address. If the parameter is not encountered
 *  in the input, a default value is set.
 */
template<class TYPE>
class ParameterValue : public Parameter {
    private:
        TYPE* pValue;   ///< Pointer to the variable to set
        TYPE Default;   ///< A default value
    public:
        /// Constructor takes a pointer to the variable and a default value
        ParameterValue (TYPE* _pValue, TYPE _Default) : pValue(_pValue) { 
            SetDefault(_Default);
        };
        /// Destructor
        virtual ~ParameterValue () {};
  
        /// Sets the default value
        void SetDefault (TYPE _Default) { 
            Default = _Default;
            *pValue = Default; 
        };
        
        /// Reads the value from the stream and returns the next token
        virtual string Rebuild (istream& in);
}; // ParameterValue

template<class TYPE>
string ParameterValue<TYPE>::Rebuild (istream& in) {
  in >> *pValue;
  string strToken;
  in >> strToken;
  return strToken;
}

/** @brief Base class for reading a subclass of the Task class from the input. 
 *
 *  A Task class is read from the input file by invoking the creating the
 *  object and calling its Task::Rebuild method. The parameters of the newly 
 *  built Task have to be enclosed in curly brackets.
 *  The new Task is then attached to its 'parent'.
 *
 *  The abstract method New should actually create the object. ParameterTask 
 *  and DiagnosticTask overwrite this method.
 */
template<class TYPE>
class ParameterTask_BASE : public Parameter {
    protected:
        Task* pParent;    ///< Pointer to the parent of the Task
    public:
        /// Constructor takes a pointer to the parent Task
        ParameterTask_BASE (Task* _pParent) : pParent(_pParent) {}
        /// Destructor
        virtual ~ParameterTask_BASE () {}

        /** @brief Creates a new Task, calls its Task::Rebuild method 
         *  and attaches it to its parent Task
         */
        virtual string Rebuild (istream& in);
        /// Abstract method that should return a pointer to a new task
        virtual TYPE* New ()=0;
}; // ParameterTask_BASE

template<class TYPE>
string ParameterTask_BASE<TYPE>::Rebuild (istream& in) {
  // in die Liste der Subtasks des Elterntasks eintragen
  TYPE* pTask = New();
  pParent->AttachTask(pTask);
  // Task aus der Datei wiederherstellen
  string strToken;
  in >> strToken;
  if ((strToken != "{") || in.eof()) return strToken;
  return pTask->Rebuild(in);
}

/** @brief An implementation of ParameterTask_BASE for reading in most Tasks.
 *
 *  The New method is overwritten to create the Task object. For normal Tasks
 *  the default constructor ca be called.
 */
template<class TYPE>
class ParameterTask : public ParameterTask_BASE<TYPE> {
    public:
        /// Constructor takes a pointer to the parent Task
        ParameterTask (Task* _pParent) : ParameterTask_BASE<TYPE>(_pParent) {}
        /// Destructor
        virtual ~ParameterTask () {}

        /// Creates a new Task by calling the new operator
        virtual TYPE* New () { return new TYPE(); }
}; // ParameterTask


#endif // PARAMETER_H

