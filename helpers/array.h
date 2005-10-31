#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <assert.h>

using namespace std;


// ---------------------------------------------------------
// declarations

template<class T>
class Array;

/** Array of fixed size */
template<class T>
class Array {
protected:
  T* arr;
  int len;

public:
  /** constructs an Array of length len and copies the
      contents of arr into it if provided
      @memo constructor */
  Array(int len, const T* arr = 0);

  /** copy constructor */
  Array(const Array<T>&);

  /** constructs Array of length 1 containing {t0}
      @memo named constructor */
  static Array set (T t0);

  /** constructs Array of length 2 containing {t0,t1}
      @memo named constructor */
  static Array set (T t0, T t1);

  /** constructs Array of length 3 containing {t0,t1,t2} 
      @memo named constructor */
  static Array set (T t0, T t1, T t2);

  /** constructs Array of length 4 containing {t0,t1,t2,t3} 
      @memo named constructor */
  static Array set (T t0, T t1, T t2, T t3);

  /** constructs Array of length 5 containing {t0,t1,t2,t3,t4} 
      @memo named constructor */
  static Array set (T t0, T t1, T t2, T t3, T t4);

  /** destructor */
  ~Array();

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

  /** index operator for writing
      @memo writing */
  T& operator[](int pos);

  /** index operator for reading
      @memo reading */
  T operator[](int pos) const;

  /** assign */
  Array<T>& operator=(const Array<T>&);

  /** comparison */
  int operator==(const Array<T>&) const;

};

/**  puts reasonable output for Array in ostream

  Example:
  \begin{verbatim}
  Array a.set(2,5);
  cout << a;
  \end{verbatim}
  @memo output for array */

template<class T>
ostream& operator<< (ostream&, const Array<T>&);

// --------------------------------------------

/** Array of fixed size which can do elementary mathematical operations */
template<class T>
class NumArray : public Array<T> {
public:
  /** constructs an NumArray of length len and copies the
      contents of arr into it if provided 
      @memo constructor */
  NumArray(int len, const T* arr = 0);

  /** copy constructor */
  NumArray(const Array<T>&);

  /** constructs NumArray of length 1 containing {t0}
      @memo named constructor */
  static NumArray set (T);

  /** constructs NumArray of length 2 containing {t0,t1}
      @memo named constructor */
  static NumArray set (T, T);

  /** constructs NumArray of length 3 containing {t0,t1,t2} 
      @memo named constructor */
  static NumArray set (T, T, T);

  /** constructs NumArray of length 4 containing {t0,t1,t2,t3} 
      @memo named constructor */
  static NumArray set (T t0, T t1, T t2, T t3);

  /** constructs NumArray of length 5 containing {t0,t1,t2,t3,t4} 
      @memo named constructor */
  static NumArray set (T t0, T t1, T t2, T t3, T t4);

  /** set every element to zero
      @memo clear the NumArray */
  void clear();

  /** */
  NumArray<T>& operator += (const NumArray<T>&);
  /** */
  NumArray<T>& operator -= (const NumArray<T>&);
  /** */
  NumArray<T>& operator += (T);
  /** */
  NumArray<T>& operator -= (T);
  /** */
  NumArray<T>& operator *= (T);
  /** */
  NumArray<T>& operator /= (T);
  /** */
  NumArray<T> operator+ (const NumArray<T>&) const;
  /** */
  NumArray<T> operator- (const NumArray<T>&) const;
  /** */
  NumArray<T> operator+ (T) const;
  /** */
  NumArray<T> operator- (T) const;
  /** */
  NumArray<T> operator* (T) const;
  /** */
  NumArray<T> operator/ (T) const;
};

// -------------------------------------------------------------
// inline functions

template<class T>
inline T& Array<T>::operator[](int pos)
{
  return arr[pos];
}

template<class T>
inline T Array<T>::operator[](int pos) const
{
  return arr[pos];
}

template<class T>
inline const T* Array<T>::data() const
{
  return arr;
}

template<class T>
inline T* Array<T>::data()
{
  return arr;
}

template<class T>
inline Array<T>::operator const T*() const
{
  return data();
}

template<class T>
inline Array<T>::operator T*()
{
  return data();
}

/** typedef NumArray<int> Index
    @memo Integer NumArray */
typedef NumArray<int> Index;

#include "array.t"

#endif
