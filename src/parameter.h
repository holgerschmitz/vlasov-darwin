// -*- C++ -*-
// $Id$

/** @file parameter.h
 *  The classes in this file are used for reading parameters from
 *  an input file. The parameters are stored in a parameter map
 */

#include <map>
#include <string>
#include <iostream>
#include "stlpwrapper.h"

#ifndef PARAMETER_H
#define PARAMETER_H

class Task;

class Parameter;
typedef PtrWrapper<Parameter> WParameter;
typedef std::map<std::string,WParameter> PARAMETERMAP;

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
        virtual std::string Rebuild (std::istream& in) = 0; // abstract class Parameter

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
        virtual std::string Rebuild(std::istream& in);
}; // ParameterValue


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
template<class Type, class BaseType>
class ParameterTask : public Parameter {
	public:
    protected:
        BaseType **value; ///< Pointer to the pointer of the object
    public:
        /// Constructor takes a pointer to the parent Task
        ParameterTask (BaseType **value_) : value(value_){}
        /// Destructor
        virtual ~ParameterTask () {}

        /** @brief Creates a new Task, calls its Task::Rebuild method 
         *  and attaches it to its parent Task
         */
        virtual std::string Rebuild (std::istream& in);
        /// Abstract method that should return a pointer to a new task
        virtual BaseType* NewInstance ()  { return new Type; }
}; // ParameterTask_BASE


#include "parameter.t"

#endif // PARAMETER_H

