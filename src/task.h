// -*- C++ -*-
// $Id$

/** @file task.h
 *  Contains the declarations of the fundamental Rebuildable class and the 
 *  also important Process class
 */

#include <list>
#include <vector>
#include <string>
#include "parameter.h"

#ifndef TASK_H
#define TASK_H


/** @brief 
 */
class Rebuildable {
    public:
        Rebuildable () {}
        virtual ~Rebuildable () {} ///< virtual destructor

    protected:
        virtual PARAMETERMAP* MakeParamMap (PARAMETERMAP* pm = NULL);
    public:
        /** @brief Rebuilds the task from the setup. This normally does not need 
         *  to be overwritten.
         *
         *  Rebuild always returns the next token that does not belong to the
         *  object setup.
         */
        virtual std::string Rebuild (std::istream& in);
        virtual void finalize() {};
}; // Rebuildable


#endif // TASK_H

