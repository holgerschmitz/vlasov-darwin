// -*- C++ -*-
// $Id$

#ifndef INDEX_H
#define INDEX_H

#include <assert.h>
#include <iostream>

// using namespace std;

/** @brief A class template to implement fixed length arrays of various
 *  types and standard operations on them
 */
template<class T, int rank>
class FixedArray {
  private:
      /// The internal data array
      T data[rank];
  public:
      /// Return a reference to the i-th data element
      T &operator[](int i) { return data[i]; }
      /// Return a const reference to the i-th data element
      const T &operator[](int i) const { return data[i]; }

      /// Default constructor
      FixedArray() {}
      
      /** @brief Single value constructor.
       *  The whole array is filled with a single value
       */
      FixedArray(const T& a) {
          for (int i=0; i<rank; ++i) data[i] =a;
      }

      /// Constructor for a rank 2 array
      FixedArray(const T& a1, const T& a2) {
          assert(rank==2);
          data[0]=a1;
          data[1]=a2;
      }

      /// Constructor for a rank 3 array
      FixedArray(const T& a1, const T& a2, 
                 const T& a3) {
          assert(rank==3);
          data[0]=a1;
          data[1]=a2;
          data[2]=a3;
      }

      /// Constructor for a rank 4 array
      FixedArray(const T& a1, const T& a2, 
                 const T& a3, const T& a4) {
          assert(rank==4);
          data[0]=a1;
          data[1]=a2;
          data[2]=a3;
          data[3]=a4;
      }

      /// Constructor for a rank 5 array
      FixedArray(const T& a1, const T& a2, 
                 const T& a3, const T& a4,
                 const T& a5) {
          assert(rank==5);
          data[0]=a1;
          data[1]=a2;
          data[2]=a3;
          data[3]=a4;
          data[4]=a5;
      }

      /// Copy constructor
      FixedArray(const FixedArray& A) {
          for (int i=0; i<rank; ++i) data[i] = A[i];
      }

      /// Assignement operator, copies the values
      FixedArray &operator=(const FixedArray& A) {
          for (int i=0; i<rank; ++i) data[i] = A[i];
          return *this;
      }

      /** @brief Addition operator.
       *  Adds the elements and returns a new array
       */
      template<class T2>
      FixedArray operator+(const FixedArray<T2,rank> & A) {
          FixedArray<T,rank> result;
          for (int i=0; i<rank; ++i) result[i] = data[i] + A[i];
          return result;
      }

      /** @brief Subtraction operator.
       *  Subtracts the elements and returns a new array
       */
      template<class T2>
      FixedArray operator-(const FixedArray<T2,rank>& A) {
          FixedArray<T,rank> result;
          for (int i=0; i<rank; ++i) result[i] = data[i] - A[i];
          return result;
      }

      /** @brief Multiplication operator.
       *  Multiplies element by element and returns a new array
       */
      template<class T2>
      FixedArray operator*(const FixedArray<T2,rank>& A) {
          FixedArray<T,rank> result;
          for (int i=0; i<rank; ++i) result[i] = data[i] * A[i];
          return result;
      }

      /** @brief Division operator.
       *  Divides element by element and returns a new array
       */
      template<class T2>
      FixedArray operator/(const FixedArray<T2,rank>& A) {
          FixedArray<T,rank> result;
          for (int i=0; i<rank; ++i) result[i] = data[i] / A[i];
          return result;
      }

      /** @brief Scalar multiplication operator.
       *  Multiplies all the element with the scalar and returns 
       *  a new array
       */
      FixedArray operator*(T A) {
          FixedArray<T,rank> result;
          for (int i=0; i<rank; ++i) result[i] = data[i]*A;
          return result;
      }

      /** @brief Scalar division operator.
       *  Divides all the element by the scalar and returns 
       *  a new array
       */
      FixedArray operator/(T A) {
          FixedArray<T,rank> result;
          for (int i=0; i<rank; ++i) result[i] = data[i]/A;
          return result;
      }

      /** @brief Scalar addition operator.
       *  Adds a scalar to all the element and returns 
       *  a new array
       */
      FixedArray operator+(T A) {
          FixedArray<T,rank> result;
          for (int i=0; i<rank; ++i) result[i] = data[i]+A;
          return result;
      }

      /** @brief Scalar subtraction operator.
       *  Subtracts a scalar from all the element and returns 
       *  a new array
       */
      FixedArray operator-(T A) {
          FixedArray<T,rank> result;
          for (int i=0; i<rank; ++i) result[i] = data[i]-A;
          return result;
      }

      /** @brief Returns the product of all the elements
       */
      T product() {
          T Result=1;
          for (int i=0; i<rank; ++i) Result *= data[i];
          return Result;
      }

      /** @brief Construct using a C-array
       *  @todo This is not safe and should be replaced by something 
       *  better
       */
      FixedArray(const T *A) {
          for (int i=0; i<rank; ++i) data[i] = A[i];
      }
      /** @brief Assign a C-array to a FixedArray
       *  @todo This is not safe and should be replaced by something 
       *  better
       */
     FixedArray &operator=(T *A) {
          for (int i=0; i<rank; ++i) data[i] = A[i];
          return *this;
      }
      
      /** @brief Get a pointer to the internal data array
       *  @todo This is not very nice and only for compatibility with the 
       *  Matrix class
       */
      T *Data() { return data; }

      /** @brief Cast to a pointer to the internal data array
       *  @todo This is not very nice and only for compatibility with the 
       *  Matrix class
       */
      operator T*()  { return data; }

      /// Return the static length of the array
      int length() { return rank; }
};

/** @brief Write a FixedArray to the output stream
 */
template<class T, int rank>
ostream &operator<<(ostream &O, FixedArray<T,rank> A) {
    O << A[0];
    for (int i=1; i<rank; ++i) O << ", " << A[i];
    return O;
}

#endif
