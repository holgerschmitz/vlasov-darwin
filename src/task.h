// -*- C++ -*-
// $Id$

/** @file task.h
 *  Contains the declarations of the fundamental Task class and the 
 *  also important Process class
 */

#include <list>
#include <vector>
#include <string>
#include "parameter.h"

#ifndef TASK_H
#define TASK_H


/** @brief The Task class is the fundamental class from which all the 
 *  classes that carry out simulation calculations are derived.
 *
 *  The Task contains information of when it is to be executed.
 *  Also every object instantiating Task has a name.
 *  A Task can also include subtasks which are executed whenever the 'parent' Task
 *  is executed.
 *  @param name The name of the Task
 *  @param start The timestep when the Task starts execution
 *  @param end The timestep when the Task ends execution
 *  @param step The interval between execution in timesteps
 *  @param Norm The name of the referenced Normalization object
 */
class Task {
    public:
        Task () {}
        virtual ~Task () {} ///< virtual destructor

    protected:
        virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
    public:
        /** @brief Rebuilds the task from the setup. This normally does not need 
         *  to be overwritten.
         *
         *  Rebuild always returns the next token that does not belong to the
         *  object setup.
         */
        std::string Rebuild (std::istream& in);
}; // Task

#endif // TASK_H

