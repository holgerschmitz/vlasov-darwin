/***************************************************************************
                          stlpwrapper.h  -  description
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
/** @file stlpwrapper.h
 *  Contains a pointer wrapper derived from Mhumit Khan stl.html
 *  @todo implement a proper smart pointer
 */


// the user still has to delete the objects from the heap

#ifndef STLPWRAPPER_H
#define STLPWRAPPER_H

#ifdef _DEBUG
#define PtrWrapper XP
#endif

/** @brief A simple pointer wrapper
 *
 *  Wraps a pointer, so that it can be stored in an STL container.
 *  The current implementation has however little advantages
 */
template <class TYPE>
class PtrWrapper {
    private:
        TYPE* m_pObj;   ///< The pointer
    public:
        /// Constructor witha pointer to the TYPE
        PtrWrapper(TYPE* pObj = 0) : m_pObj(pObj) {};
        /// Copy constructor
        PtrWrapper(const PtrWrapper<TYPE>& wrapper) { 
            m_pObj = wrapper.m_pObj;
        };
        /// Assignment operator
        PtrWrapper<TYPE>& operator= (const PtrWrapper<TYPE>& wrapper) {
            m_pObj = wrapper.m_pObj;
            return *this;
        }  
        /// Destructor doesn't free pointer
        ~PtrWrapper() {};
        
        //@{
        operator const TYPE* () const { return m_pObj; };
        operator TYPE* () { return m_pObj; };
        ///< Typecast to the original pointer type
        //@}

        /// Explicit method for the typecast
        TYPE* pObj () const { return m_pObj; };

        /// Compares the two objects
        friend bool operator== (const PtrWrapper<TYPE>& left, const PtrWrapper<TYPE>& right) {
            return (left.m_pObj == right.m_pObj);
        };

        /// Compares the two objects
        friend bool operator!= (const PtrWrapper<TYPE>& left, const PtrWrapper<TYPE>& right) {
            return (left.m_pObj != right.m_pObj);
        };

        /// Compares the two objects
//          friend bool operator< (const PtrWrapper<TYPE>& left, const PtrWrapper<TYPE>& right) {
//              return (*(left.m_pObj) < *(right.m_pObj));
//          };
};

#endif // STLPWRAPPER_H
