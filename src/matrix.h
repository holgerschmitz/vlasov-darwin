#ifndef MATRIX_H
#define MATRIX_H

#include <assert.h>
#include <iostream>
#include "util.h"
#include "array.h"

using namespace std;
template<class T, int rank>
class Matrix;

template<class T, int rank>
ostream& operator<< (ostream&, const Matrix<T, rank>&);


template<class T>
class ValueChecker {
  protected:
    T& checkValue(T &val) const {
//        assert( (val<1) || (val>0) );
        return val;
    }
};

// ---------------------------------------------------------
// declarations

/** An elementary matrix class */
template<class T, int rank>
class Matrix : public ValueChecker<T> {
  protected:
    T* matr;
    T* matr_fast;
    void* index;
    int size;
    int lo[rank];
    int hi[rank];
    int dims[rank];

  public:
    /** default constructur */
    Matrix();

    /** constructor, which builds Matrix of size d[0] x ... x d[rank-1]

        Example: 
        \begin{verbatim}
        int size[2]={512,512};
        Matrix<double,2> m(size);
        \end{verbatim}
        @memo constructor */
    Matrix(const int d[rank]);

    /** constructor, which builds Matrix with lower indices l[0],...,l[rank-1]
        and upper indices h[0],...,h[rank-1]

        Example: 
        \begin{verbatim}
        int l[2]={-5,-10};
        int h[2]={15,36}
        Matrix<double,2> m(l,h);
        \end{verbatim}
        @memo constructor */
    Matrix(const int l[rank], const int h[rank]);

    /** copy constructor */
    Matrix(const Matrix<T, rank>&);

    /** destructor */
    ~Matrix();

    /** */
    const int* getLow() const;
    /** */
    const int* getHigh() const;
    /** */
    const int* getDims() const;
    /** */
    int getLow(int) const;
    /** */
    int getHigh(int) const;
    /** */
    int getDims(int) const;
    /** */
    int getSize() const;
    /** return pointer to the first element holding the data, 
        modifying allowed */
    T* data();
    /** return pointer to the first element holding the data, 
        modifying forbidden */
    const T* data() const;
    /** return pointer to the first element holding the data, 
        modifying allowed */
    operator T*();
    /** return pointer to the first element holding the data, 
        modifying forbidden */
    operator const T*() const;

    /** resizes matrix to size d[0] x ... x d[rank-1] */
    void resize(const int* d);
    /** resizes to matrix with lower indices l[0],...,l[rank-1]
        and upper indices h[0],...,h[rank-1] */
    void resize(const int* l, const int* h);
    /** resizes to values of matrix m */
    void resize(const Matrix<T, rank>& m);

    /** index operator, writing */
    T& operator[](const int* pos); // write
    /** index operator, reading */
    T  operator[](const int* pos) const; // read
    /** index operator, writing */
    T& operator()(int i);
    /** index operator, reading */
    T  operator()(int i) const;
    /** index operator, writing */
    T& operator()(int i, int j);
    /** index operator, reading */
    T  operator()(int i, int j) const;
    /** index operator, writing */
    T& operator()(int i, int j, int k);
    /** index operator, reading */
    T  operator()(int i, int j, int k) const;
    /** index operator, writing */
    T& operator()(int i, int j, int k, int l);
    /** index operator, reading */
    T  operator()(int i, int j, int k, int l) const;
    /** index operator, writing */
    T& operator()(int i, int j, int k, int l, int m);
    /** index operator, reading */
    T  operator()(int i, int j, int k, int l, int m) const;

    /** assign */
    Matrix<T, rank>& operator=(const Matrix<T, rank>&);
    /** comparison */
    int operator==(const Matrix<T, rank>&) const;
    /** comparison */
    int operator!=(const Matrix<T, rank>&) const;

    /** output to ostream */
    friend ostream& operator<< <>(ostream&, const Matrix<T, rank>&);

  protected:
    /** */
    void deleteData();
    /** */
    void newData(const int* l, const int* h);
};


// ---------------------------------------------------------

/** An matrix class which can do 
    elementary mathematical operations */
template<class T, int rank>
class NumMatrix : public Matrix<T, rank> {
public:
  /** default constructur */
  NumMatrix();

  /** constructor, which builds NumMatrix of size d[0] x ... x d[rank-1]
      @memo constructor */
  NumMatrix(const int d[rank]);

    /** constructor, which builds NumMatrix with lower 
	indices l[0],...,l[rank-1]
	and upper indices h[0],...,h[rank-1]
	@memo constructor */
  NumMatrix(const int l[rank], const int h[rank]);

  /** copy constructor */
  NumMatrix(const Matrix<T, rank>&);

  /** setting all elements of NumMatrix to zero */
  void clear();
  /** maximum norm of all elements of NumMatrix */
  double max_norm();

  /** mean value of direct neighborhood */
  T mean(const int *pos) const;
  /** mean value of direct neighborhood */
  T mean(int i) const;
  /** mean value of direct neighborhood */
  T mean(int i, int j) const;
  /** mean value of direct neighborhood */
  T mean(int i, int j, int k) const;

  /** interpolates value at pos+diff/r (linear)
      @memo linear interpolation */
  T interpol(const int* pos, const int* diff, int r) const;
  /** interpolates value at pos+diff/r (cubic)
      @memo cubic interpolation */
  T interpol3(const int* pos, const int* diff, int r) const;

  /** */
  NumMatrix<T, rank>& operator += (const NumMatrix<T, rank>&);
  /** */
  NumMatrix<T, rank>& operator -= (const NumMatrix<T, rank>&);
  /** */
  NumMatrix<T, rank>& operator *= (T);
  /** */
  NumMatrix<T, rank>& operator /= (T);

  /** */
  NumMatrix<T, rank> operator+(const NumMatrix<T, rank>&) const;
  /** */
  NumMatrix<T, rank> operator-(const NumMatrix<T, rank>&) const;
  /** */
  NumMatrix<T, rank> operator*(T) const;
  /** */
  NumMatrix<T, rank> operator/(T) const;
};

// -------------------------------------------------------------
// inline functions

template<class T, int rank>
inline T& Matrix<T, rank>::operator[](const int* pos)
{
  assert((pos[rank-1] >= lo[rank-1]) && (pos[rank-1] <= hi[rank-1]));
  int p = pos[rank-1]-lo[rank-1];
      
  for(int i = rank-2; i >= 0; i--) {
    assert((pos[i] >= lo[i]) && (pos[i] <= hi[i]));
    p = p*dims[i] + (pos[i]-lo[i]);
  }

  assert((p >= 0) && (p < size));

  return matr[p];
}

template<class T, int rank>
inline T Matrix<T, rank>::operator[](const int* pos) const
{
  assert((pos[rank-1] >= lo[rank-1]) && (pos[rank-1] <= hi[rank-1]));
  int p = pos[rank-1];
      
  for(int i = rank-2; i >= 0; i--) {
    assert((pos[i] >= lo[i]) && (pos[i] <= hi[i]));
    p = p*dims[i] + pos[i];
  }

  return matr_fast[p];
}

template<class T, int rank> 
inline T* Matrix<T, rank>::data()
{
  return matr;
}

template<class T, int rank> 
inline const T* Matrix<T, rank>::data() const
{
  return matr;
}

template<class T, int rank>
inline Matrix<T, rank>::operator T*()
{
  return data();
}

template<class T, int rank>
inline Matrix<T, rank>::operator const T*() const
{
  return data();
}

template<class T, int rank>
inline  T& Matrix<T, rank>::operator ()(int i)
{
  assert(rank == 1);
  assert((i >= lo[0]) && (i <= hi[0]));
  return checkValue(matr_fast[i]);
}

template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i) const
{
  assert(rank == 1);
  assert((i >= lo[0]) && (i <= hi[0]));
  return checkValue(matr_fast[i]);
}

template<class T, int rank>
inline  T& Matrix<T, rank>::operator ()(int i, int j)
{
  assert(rank == 2);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  return checkValue(matr_fast[i + dims[0] * j]);
}

template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i, int j) const
{
  assert(rank == 2);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  return checkValue(matr_fast[i + dims[0] * j]);
}

template<class T, int rank>
inline  T& Matrix<T, rank>::operator ()(int i, int j, int k) 
{
  assert(rank == 3);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  assert((k >= lo[2]) && (k <= hi[2]));
  return checkValue(matr_fast[i + dims[0]*(j + dims[1]*k)]);
}

template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i, int j, int k) const
{
  assert(rank == 3);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  assert((k >= lo[2]) && (k <= hi[2]));
  return checkValue(matr_fast[i + dims[0]*(j + dims[1]*k)]);
}

template<class T, int rank>
inline  T& Matrix<T, rank>::operator ()(int i, int j, int k, int l) 
{
  assert(rank == 4);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  assert((k >= lo[2]) && (k <= hi[2]));
  assert((l >= lo[3]) && (l <= hi[3]));
  return checkValue(matr_fast[i + dims[0]*(j + dims[1]*(k+dims[2]*l))]);
}

template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i, int j, int k, int l) const
{
  assert(rank == 4);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  assert((k >= lo[2]) && (k <= hi[2]));
  assert((l >= lo[3]) && (l <= hi[3]));
  return checkValue(matr_fast[i + dims[0]*(j + dims[1]*(k+dims[2]*l))]);
}

template<class T, int rank>
inline  T& Matrix<T, rank>::operator ()(int i, int j, int k, int l, int m) 
{
  assert(rank == 5);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  assert((k >= lo[2]) && (k <= hi[2]));
  assert((l >= lo[3]) && (l <= hi[3]));
  assert((m >= lo[4]) && (m <= hi[4]));
  return checkValue(matr_fast[i + dims[0]*(j + dims[1]*(k+dims[2]*(l+dims[3]*m)))]);
}

template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i, int j, int k, int l, int m) const
{
  assert(rank == 5);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  assert((k >= lo[2]) && (k <= hi[2]));
  assert((m >= lo[4]) && (m <= hi[4]));
  return checkValue(matr_fast[i + dims[0]*(j + dims[1]*(k+dims[2]*(l+dims[3]*m)))]);
}

#include "matrix.t"

#endif
