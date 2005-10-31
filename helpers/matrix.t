#include <algorithm>
#include <math.h>
#include <stdlib.h>

// --------------------------------------------------------------
// implementation

template<class T, int rank>
Matrix<T, rank>::Matrix()
{
  matr = 0;
  size = 0;
  index = 0;
}

template<class T, int rank>
Matrix<T, rank>::Matrix(const int* d)
{
  matr = 0;
  size = 0;
  
  Index l(rank);
  Index h(rank);
  l.clear();
  for (int i = 0; i < rank; i++) 
    h[i] = d[i] - 1;
  
  newData(l, h);
}

template<class T, int rank>
Matrix<T, rank>::Matrix(const int* l, const int* h)
{
  matr = 0;
  size = 0;
  newData(l, h);
}

template<class T, int rank>
Matrix<T, rank>::Matrix(const Matrix<T, rank>& _matr)
{
  matr = 0;
  size = 0;

  newData(_matr.lo, _matr.hi);

  for (int i = 0; i < size; i++)    
    matr[i] = checkValue(_matr.matr[i]);
}
 
template<class T, int rank> 
Matrix<T, rank>::~Matrix()
{
  deleteData();
}

template<class T, int rank>
void Matrix<T, rank>::swap(Matrix<T, rank>& M)
{
  T* matr_tmp=matr;             matr=M.matr;              M.matr=matr_tmp;
  T* matr_fast_tmp=matr_fast;   matr_fast=M.matr_fast;    M.matr_fast=matr_fast_tmp;
  void* index_tmp=index;        index=M.index;            M.index=index_tmp;
  int size_tmp=size;            size=M.size;              M.size=size_tmp;
  
  for(int i=0; i<rank; ++i)
  {
    int tmp;
    tmp = lo[i];      lo[i]=M.lo[i];      M.lo[i]=tmp;
    tmp = hi[i];      hi[i]=M.hi[i];      M.hi[i]=tmp;
    tmp = dims[i];    dims[i]=M.dims[i];  M.dims[i]=tmp;
  }
   
}

template<class T, int rank> 
const int* Matrix<T, rank>::getLow() const
{
  return lo;
}

template<class T, int rank> 
const int* Matrix<T, rank>::getHigh() const
{
  return hi;
}

template<class T, int rank> 
const int* Matrix<T, rank>::getDims() const
{
  return dims;
}

template<class T, int rank>
int Matrix<T, rank>::getLow(int i) const
{
  return lo[i];
}

template<class T, int rank>
int Matrix<T, rank>::getHigh(int i) const
{
  return hi[i];
}

template<class T, int rank>
int Matrix<T, rank>::getDims(int i) const
{
  return dims[i];
}

template<class T, int rank>
int Matrix<T, rank>::getSize() const
{
  return size;
}

template<class T, int rank>
Matrix<T, rank>& Matrix<T, rank>::operator=(const Matrix<T, rank>& _matr)
{
  for (int d = 0; d < rank; d++) {
    if (_matr.lo[d] != lo[d] ||
	    _matr.hi[d] != hi[d]) {
      deleteData();
      newData(_matr.lo, _matr.hi);
      break;
    }
  }

  for (int i = 0; i < size; i++)
    matr[i] = checkValue(_matr.matr[i]);
    
  return *this;
}
 
template<class T, int rank>
int Matrix<T, rank>::operator==(const Matrix<T, rank>& _matr) const
{
  int i;
  for (i = 0; i < rank; i++) {
    if (lo[i] != _matr.lo[i])
      return 0;
    if (hi[i] != _matr.hi[i])
      return 0;
  }

  for (i = 0; i < size; i++)
    if (matr[i] != _matr.matr[i])
      return 0;
  return 1;
}

template<class T, int rank>
int Matrix<T, rank>::operator!=(const Matrix<T, rank>& _matr) const
{
  return !operator==(_matr);
}

template<class T, int rank>
void Matrix<T, rank>::resize(const int* d)
{
  Index l(rank);
  Index h(rank);
  l.clear();
  for (int i = 0; i < rank; i++) 
    h[i] = d[i] - 1;

  resize(l, h);
}

template<class T, int rank>
void Matrix<T, rank>::resize(const int* l, const int* h)
{
  bool needresize = false;
  for (int i = 0; i < rank; i++)
  {
    needresize = needresize || (l[i] != lo[i]) || (h[i] != hi[i]);
  }
  
  if (needresize)
  {
    deleteData();
    newData(l, h);
  }
}

template<class T, int rank>
void Matrix<T, rank>::resize(const Matrix<T, rank>& _matr)
{
  resize(_matr.lo, _matr.hi);
}

template<class T, int rank>
void Matrix<T, rank>::deleteData()
{
  if (matr)
    delete[] matr;
  matr = 0;
  size = 0;
}

template<class T, int rank>
void Matrix<T, rank>::newData(const int* l, const int* h)
{
  size = 1;
  int d;
  for (d = 0; d < rank; d++) {
    assert(l[d] <= h[d]);
    lo[d] = l[d];
    hi[d] = h[d]; 
    dims[d] = h[d] - l[d] + 1;
    size *= dims[d];
  }
  matr = new T[size];
  if (matr==NULL) std::cerr << "Could not allocate memory!\n";
  int p = -lo[0];

  for (d = 1; d < rank ; d++) {
    p = p*dims[d] - lo[d];
  }
  matr_fast = matr + p;
}

// template<class T, int rank>
// ostream& operator<<(ostream& os, const Matrix<T, rank>& matr) 
// {
//   if (!matr) {
//     os << "(undef)" << endl;
//     return os;
//   }
//   os << "ostream& operator<<(ostream& os, const Matrix<T, rank>& matr)"
//      << endl << " not implemented" << endl;
//   return os;
// }
// 
template<class T>
ostream& operator<<(ostream& os, const Matrix<T, 2>& matr) 
{
  if (!matr) {
    os << "(undef)" << endl;
    return os;
  }
  
  for (int i = matr.getLow(0); i<= matr.getHigh(0); ++i)
  {
    for (int j = matr.getLow(1); j<= matr.getHigh(1); ++j)
    {
      os << i << " " << j << " " << matr(i,j)<< "\n";
    }
    os << std::endl;
  }
  return os;
}


// -----------------------------------------------------------------

template<class T, int rank>
NumMatrix<T, rank>::NumMatrix() : Matrix<T, rank>() {}

template<class T, int rank>
NumMatrix<T, rank>::NumMatrix(const int* d) : Matrix<T, rank>(d) {}

template<class T, int rank>
NumMatrix<T, rank>::NumMatrix(const int* l, const int* h)
  : Matrix<T, rank>(l, h) {}

template<class T, int rank>
NumMatrix<T, rank>::NumMatrix(const Matrix<T, rank>& _matr)
  : Matrix<T, rank>(_matr) {}

template<class T, int rank>
void NumMatrix<T, rank>::clear()
{ 
  for (int i = 0; i < this->size; i++)
    this->matr[i] = 0;
}

template<class T, int rank>
T NumMatrix<T, rank>::interpol(const int* pos, const int* diff, int r) const
{
  int dual, i;
  Index p(rank);

  T res = 0;

  for (dual = 0; dual < (1 << rank); dual++) {
    int fak = 1;
    for (i = 0; i < rank; i++) {
      if (dual & (1 << i)) {
	fak *= r - diff[i];
	p[i] = pos[i];
      } else {
	fak *= diff[i];
	p[i] = pos[i] + 1;
      }
    }
    res += T(double(fak) * (*this)[p]);
  }

  for (i = 0; i < rank; i++)
    res /= r;
  return res;
}

// template<class T, int rank>
// T NumMatrix<T, rank>::interpol3(const int* pos, const int* diff, int r) const
// {
//   cerr << "fatal: interpol3 not generally implemented!" << endl;
//   return (T)0;
// }

template<class T, int rank>
NumMatrix<T, rank>& NumMatrix<T, rank>::operator*=(T t)
{
  for (int i = 0; i < this->size; i++)
    this->matr[i] *= t;
  return *this;
}

template<class T, int rank>
NumMatrix<T, rank>& NumMatrix<T, rank>::operator/=(T t)
{
  T div = T(1./t);
  for (int i = 0; i < this->size; i++)
    this->matr[i] *= div;
  return *this;
}

template<class T, int rank>
NumMatrix<T, rank>& NumMatrix<T, rank>::operator+=(const NumMatrix<T, rank>& 
						   _matr)
{
  assert(this->size == _matr.size);
  for (int i = 0; i < this->size; i++)
    this->matr[i] += _matr.matr[i];
  return *this;
}

template<class T, int rank>
NumMatrix<T, rank>& NumMatrix<T, rank>::operator-=(const NumMatrix<T, rank>& 
						   _matr)
{
  assert(this->size == _matr.size);
  for (int i = 0; i < this->size; i++)
    this->matr[i] -= _matr.matr[i];
  return *this;
}

template<class T, int rank>
NumMatrix<T, rank> NumMatrix<T, rank>::
  operator+(const NumMatrix<T, rank>& matr1) const
{
  assert(this->size == matr1.size);
  NumMatrix<T, rank> tmp = (*this);
  tmp += matr1; 
  return tmp;
}

template<class T, int rank>
NumMatrix<T, rank> NumMatrix<T, rank>::
  operator-(const NumMatrix<T, rank>& matr1) const
{
  assert(this->size == matr1.size);
  NumMatrix<T, rank> tmp = (*this);
  tmp -= matr1; 
  return tmp;
}

template<class T, int rank>
NumMatrix<T, rank> NumMatrix<T, rank>::operator*(T t) const
{
  NumMatrix<T, rank> tmp = (*this);
  tmp *= t; 
  return tmp;
}

template<class T, int rank>
NumMatrix<T, rank> NumMatrix<T, rank>::operator/(T t) const
{
  NumMatrix<T, rank> tmp = (*this);
  tmp /= t; 
  return tmp;
}


// ------------------------------------------------------------
// specializations

inline double interpol3_w(int i, double x)
{
  switch (i) {
  case 0: return((x*x*x-6*x*x+11*x-6)/(-6));
  case 1: return((x*x*x-5*x*x+6*x)/2);
  case 2: return((x*x*x-4*x*x+3*x)/(-2));
  case 3: return((x*x*x-3*x*x+2*x)/6);
  }
  return -1;
}

