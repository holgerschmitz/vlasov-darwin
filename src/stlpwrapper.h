// -*- C++ -*-
// $Id$


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
        /// The pointer
        TYPE* m_pObj;   
    public:
        /// Constructor with a pointer to the TYPE
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
        operator const TYPE* () const { return m_pObj; }
        operator TYPE* () { return m_pObj; }
        TYPE* operator->() { return m_pObj; }
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
