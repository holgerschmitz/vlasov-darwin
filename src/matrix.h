#ifndef MATRIX_H
#define MATRIX_H
	
#include <assert.h>
#include <iostream>
#include "util.h"
#include "array.h"
/** @file matrix.h
   * @brief elementary matrix class templates
   *
   * Defines elementary matrix templates, Matrix and NumMatrix. NumMatrix can do basic mathematical operations.
   */
using namespace std;

//--------------------------------------------------------
/** @brief an elementary matrix class template
 */
template<class T, int rank>
class Matrix;

template<class T, int rank>
ostream& operator<< (ostream&, const Matrix<T, rank>&);

/// @brief class for testing on "nan".
template<class T>
class ValueChecker {
  protected:
  ///can be used to check T, for example on "nan", does nothing now
    T& checkValue(T &val) const {
//        assert( (val<1) || (val>0) );
        return val;
    }
};

// ---------------------------------------------------------
// declarations
/** @brief an elementary matrix class template
 */
template<class T, int rank>
class Matrix : public ValueChecker<T> {
  protected:
    /**   pointer to elements */
    T* matr;
    /** pointer to matrix elements, but starts at index 0. 
    (no need to compute local element from global indices)*/
    T* matr_fast;
    /** index pointer */
    void* index;
    /**  number of elements*/
    int size;
    /**  lower boundaries, position in a global matrix */
    int lo[rank];
    /**  upper boundaries, position in a global matrix */
    int hi[rank];
    /**   dimension in each direction*/
    int dims[rank];

  public:
//--------------------------------------------------------------
//constructors & destructor
    /** default constructur */
    Matrix();

    /** constructor, which builds Matrix of size d[0] x ... x d[rank-1]

        Example: 
        \verbatim
        int size[2]={512,512};
        Matrix<double,2> m(size);
        \endverbatim
        */
    Matrix(const int d[rank]);

    /** constructor, which builds Matrix with lower indices l[0],...,l[rank-1]
        and upper indices h[0],...,h[rank-1]
	  for multi-process runs, dividing a matrix into sub-matrices with the boundaries l[] and h[]
        Example: 
        \verbatim
        int l[2]={-5,-10};
        int h[2]={15,36}
        Matrix<double,2> m(l,h);
        \endverbatim
        */
    Matrix(const int l[rank], const int h[rank]);

    /** copy constructor */
    Matrix(const Matrix<T, rank>&);

    /** destructor */
    ~Matrix();
//--------------------------------------------------------------
    /**   swaps matrix */
    void swap(Matrix<T, rank>&);
    /**   returns pointer to first element of lower idex array, read only*/
    const int* getLow() const;
     /**  return pointer to first element of upper index-array, read only*/
    const int* getHigh() const;
    /**   returns (pointer to first element of) dimension array, read only*/
    const int* getDims() const;
    /**   returns lower boundary index, read only*/
    int getLow(int) const;
    /**  returns lower index i,read only*/
    int getHigh(int) const;
    /**  returns upper index i, read only*/
    int getDims(int) const;
    /**   returns number of elements, read only*/
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

//     index operator, writing 
//    T& operator[](const int* pos); // write
//    index operator, reading 
//    T  operator[](const int* pos) const; // read
    ///index operator read & write
    T& operator()(int i);
	/**index operator, read only*/
    T  operator()(int i) const;
	/**index operator read & write*/
    T& operator()(int i, int j);
	/**index operator, read only*/
    T  operator()(int i, int j) const;
	/**index operator read & write*/
    T& operator()(int i, int j, int k);
	/**index operator, read only*/
    T  operator()(int i, int j, int k) const;
	/**index operator read & write*/
    T& operator()(int i, int j, int k, int l);
	/**index operator, read only*/
    T  operator()(int i, int j, int k, int l) const;
	/**index operator read & write*/	
    T& operator()(int i, int j, int k, int l, int m);
	/**index operator, read only*/
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
    /** deletes matrix elements*/
    void deleteData();
    /**   allocates new matrix data, upper and lower boundaries set to l[] and h[]*/
    void newData(const int* l, const int* h);
};


// ---------------------------------------------------------
//NumMatrix declarations
/** @brief An matrix class which can do 
    elementary mathematical operations */
/** An matrix class which can do 
    elementary mathematical operations */
template<class T, int rank>
class NumMatrix : public Matrix<T, rank> {
public:
  /** default constructor */
  NumMatrix();

  /** constructor, which builds NumMatrix of size d[0] x ... x d[rank-1]
      */
  NumMatrix(const int d[rank]);

    /** constructor, which builds NumMatrix with lower 
	indices l[0],...,l[rank-1]
	and upper indices h[0],...,h[rank-1]
	*/
  NumMatrix(const int l[rank], const int h[rank]);

  /** copy constructor */
  NumMatrix(const Matrix<T, rank>&);

  /** setting all elements of NumMatrix to zero */
  void clear();

  /** mean value of direct neighborhood */
  T mean(const int *pos) const;
  /** mean value of direct neighborhood */
  T mean(int i) const;
  /** mean value of direct neighborhood */
  T mean(int i, int j) const;
  /** mean value of direct neighborhood */
  T mean(int i, int j, int k) const;

  /** interpolates value at pos+diff/r (linear)
      */
  T interpol(const int* pos, const int* diff, int r) const;
  /** interpolates value at pos+diff/r (cubic)
      */
  T interpol3(const int* pos, const int* diff, int r) const;

  /**   assign&add operator for NumMatrix*/
  NumMatrix<T, rank>& operator += (const NumMatrix<T, rank>&);
  /**   assign&substract operator for NumMatrix */
  NumMatrix<T, rank>& operator -= (const NumMatrix<T, rank>&);
  /**   assign&multiply operator for NumMatrix (by element)*/
  NumMatrix<T, rank>& operator *= (T);
  /**   assign&divide operator for NumMatrix (by element)*/
  NumMatrix<T, rank>& operator /= (T);
  /** addition operator for NumMatrix*/
  NumMatrix<T, rank> operator+(const NumMatrix<T, rank>&) const;
  /**  substraction operator for NumMatrix*/
  NumMatrix<T, rank> operator-(const NumMatrix<T, rank>&) const;
  /**  multiplication operator for NumMatrix (by element)*/
  NumMatrix<T, rank> operator*(T) const;
  /**  division operator for NumMatrix (by element)*/
  NumMatrix<T, rank> operator/(T) const;
};

// -------------------------------------------------------------
// inline functions

// template<class T, int rank>
// inline T& Matrix<T, rank>::operator[](const int* pos)
// {
//   assert((pos[rank-1] >= lo[rank-1]) && (pos[rank-1] <= hi[rank-1]));
//   int p = pos[0]-lo[0];
//       
//   for(int i = 1; i < rank; i--) {
//     assert((pos[i] >= lo[i]) && (pos[i] <= hi[i]));
//     p = p*dims[i] + (pos[i]-lo[i]);
//   }
// 
//   assert((p >= 0) && (p < size));
// 
//   return matr[p];
// }
// 
// template<class T, int rank>
// inline T Matrix<T, rank>::operator[](const int* pos) const
// {
//   assert((pos[rank-1] >= lo[rank-1]) && (pos[rank-1] <= hi[rank-1]));
//   int p = pos[rank-1];
//       
//   for(int i = rank-2; i >= 0; i--) {
//     assert((pos[i] >= lo[i]) && (pos[i] <= hi[i]));
//     p = p*dims[i] + pos[i];
//   }
// 
//   return matr_fast[p];
// }
/** returns pointer to elements of the matrix, read & write*/
template<class T, int rank> 
inline T* Matrix<T, rank>::data()
{
  return matr;
}
/** returns pointer to elements of the matrix, read only */
template<class T, int rank> 
inline const T* Matrix<T, rank>::data() const
{
  return matr;
}

/**  "pointer to" operator, read & write */
template<class T, int rank>
inline Matrix<T, rank>::operator T*()
{
  return data();
}
/**  "pointer to" operator, read only */
template<class T, int rank>
inline Matrix<T, rank>::operator const T*() const
{
  return data();
}

/*--------------------------------------------------------------*/
/* checking the index limits for various dimensions (debug code)*/
/* disable assert () with -DNDEBUG or #define NDEBUG*/
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
  return checkValue(matr_fast[j + dims[1] * i]);
}



template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i, int j) const
{
  assert(rank == 2);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  return checkValue(matr_fast[j + dims[1] * i]);
}

template<class T, int rank>
inline  T& Matrix<T, rank>::operator ()(int i, int j, int k) 
{
  assert(rank == 3);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  assert((k >= lo[2]) && (k <= hi[2]));
  return checkValue(matr_fast[k + dims[2]*(j + dims[1]*i)]);
}

template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i, int j, int k) const
{
  assert(rank == 3);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  assert((k >= lo[2]) && (k <= hi[2]));
  return checkValue(matr_fast[k + dims[2]*(j + dims[1]*i)]);
}

template<class T, int rank>
inline  T& Matrix<T, rank>::operator ()(int i, int j, int k, int l) 
{
  assert(rank == 4);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  assert((k >= lo[2]) && (k <= hi[2]));
  assert((l >= lo[3]) && (l <= hi[3]));
  return checkValue(matr_fast[l + dims[3]*(k + dims[2]*(j+dims[1]*i))]);
}

template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i, int j, int k, int l) const
{
  assert(rank == 4);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  assert((k >= lo[2]) && (k <= hi[2]));
  assert((l >= lo[3]) && (l <= hi[3]));
  return checkValue(matr_fast[l + dims[3]*(k + dims[2]*(j+dims[1]*i))]);
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
  return checkValue(matr_fast[m + dims[4]*(l + dims[3]*(k+dims[2]*(j+dims[1]*i)))]);
}

template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i, int j, int k, int l, int m) const
{
  assert(rank == 5);
  assert((i >= lo[0]) && (i <= hi[0]));
  assert((j >= lo[1]) && (j <= hi[1]));
  assert((k >= lo[2]) && (k <= hi[2]));
  assert((m >= lo[4]) && (m <= hi[4]));
  return checkValue(matr_fast[m + dims[4]*(l + dims[3]*(k+dims[2]*(j+dims[1]*i)))]);
}

/*-------------------------------------------------------------*/
/*include template member functions from matrix.t */
#include "matrix.t"

#endif
