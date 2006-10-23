#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <assert.h>

/** @file array.h
   * @brief elementary array classes
   *
   * Defines the elementary array class templates Array and NumArray.
   */

using namespace std;

// -------------------------------------------------------------
// declarations

template<class T>
class Array;
/** @brief An elementary array class template*/
template<class T>
class Array {
protected:
  /** pointer to data element */
  T* arr;
  /** length of array */
  int len;
public:

//--------------------------------------------------------------
// constructors & destructor
/**   @brief constructor
  *
  *	constructs an Array of length len and copies the
  *   contents of arr into it if provided
  */
  Array(int len, const T* arr = 0);

/** @brief copy constructor */
  Array(const Array<T>&);

/** @brief named constructor 
  *
  * constructs Array of length 1 containing {t0}
  */
  static Array set (T t0);

/** @brief named constructor
  *
  * constructs Array of length 2 containing {t0,t1}
  */
  static Array set (T t0, T t1);

/** @brief named constructor
  *
  * constructs Array of length 3 containing {t0,t1,t2}  
  */
  static Array set (T t0, T t1, T t2);

/** @brief named constructor
  *
  *constructs Array of length 4 containing {t0,t1,t2,t3} 
  */
  static Array set (T t0, T t1, T t2, T t3);

/** @brief named constructor
  *
  * constructs Array of length 5 containing {t0,t1,t2,t3,t4} 
  */
  static Array set (T t0, T t1, T t2, T t3, T t4);

/** destructor */
  ~Array();
//--------------------------------------------------------------
/** returns lenght of the Array */
  int getLength() const;

/** returns pointer to the Array, modifying forbidden */
  const T* data() const;

/** returns pointer to the Array, modifying allowed */
  T* data();

/** cast into C Array, modifying forbidden */
  operator const T*() const;
  
 /** cast into C Array, modifying allowed */
  operator T*();

/** index operator for writing */
  T& operator[](int pos);
  
 /** index operator for reading */
  T operator[](int pos) const;

/** assign */
  Array<T>& operator=(const Array<T>&);

/** comparison */
  int operator==(const Array<T>&) const;

};

/**  puts reasonable output for Array in ostream

  Example:
  \verbatim
  Array a.set(2,5);
  cout << a;
  \endverbatim
  */

template<class T>
ostream& operator<< (ostream&, const Array<T>&);

// -------------------------------------------------------------
//NumArray
/** @brief Array of fixed size which can do elementary mathematical operations */
template<class T>
class NumArray : public Array<T> {
public:
/** @brief constructor
  *
  * constructs an NumArray of length len and copies the
  * contents of arr into it if provided 
  */
  NumArray(int len, const T* arr = 0);

/** copy constructor */
  NumArray(const Array<T>&);

/** @brief named constructor
  *
  * constructs NumArray of length 1 containing {t0}
  */
  static NumArray set (T);

/** @brief named constructor
  *
  * constructs NumArray of length 2 containing {t0,t1}
  */
  static NumArray set (T, T);

/** @brief named constructor
  *
  * constructs NumArray of length 3 containing {t0,t1,t2} 
  */
  static NumArray set (T, T, T);

/** @brief named constructor
  *
  * constructs NumArray of length 4 containing {t0,t1,t2,t3} 
  */
  static NumArray set (T t0, T t1, T t2, T t3);

/** @brief named constructor
  *
  * constructs NumArray of length 5 containing {t0,t1,t2,t3,t4} 
  */
  static NumArray set (T t0, T t1, T t2, T t3, T t4);

  /** set every element to zero */
  void clear();

  /** sum & assign for constant NumArrays*/
  NumArray<T>& operator += (const NumArray<T>&);
  /** multiply & assign (by element) for constant NumArrays*/
  NumArray<T>& operator -= (const NumArray<T>&);
  /** sum & assign (by element)*/
  NumArray<T>& operator += (T);
  /** substract & assign (by element)*/
  NumArray<T>& operator -= (T);
  /** multiply & assign (by element)*/
  NumArray<T>& operator *= (T);
  /** divide & assign (by element)*/
  NumArray<T>& operator /= (T);
/** addition operator, using the += operator for constant NumArrays*/
  NumArray<T> operator+ (const NumArray<T>&) const;
/** subtraction operator, using the -= operator for constant NumArrays*/
  NumArray<T> operator- (const NumArray<T>&) const;
/** addition operator, using the += operator */
  NumArray<T> operator+ (T) const;
/** substraction operator, using the -= operator */
  NumArray<T> operator- (T) const;
/** multiplication operator, using the *= operator */
  NumArray<T> operator* (T) const;
/** division operator, using the /= operator */
  NumArray<T> operator/ (T) const;
};

// -------------------------------------------------------------
// inline functions

/** index operator, read & write*/
template<class T>
inline T& Array<T>::operator[](int pos)
{
  return arr[pos];
}
/** index operator, read only */
template<class T>
inline T Array<T>::operator[](int pos) const
{
  return arr[pos];
}

/** returns point to first element of data, read only*/
template<class T>
inline const T* Array<T>::data() const
{
  return arr;
}

/** returns point to first element of data, read & write*/
template<class T>
inline T* Array<T>::data()
{
  return arr;
}

/** "pointer to" operator, read only*/
template<class T>
inline Array<T>::operator const T*() const
{
  return data();
}
/** "pointer to" operator, read & write */
template<class T>
inline Array<T>::operator T*()
{
  return data();
}

/** @brief Integer NumArray
  *
  * typedef NumArray<int> Index
  */
typedef NumArray<int> Index;

//--------------------------------------------------------------
//including template member functions from array.t

#include "array.t"

#endif
