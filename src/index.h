// -*- C++ -*-
// $Id$

#ifndef INDEX_H
#define INDEX_H

#include <assert.h>
#include <iostream>

// using namespace std;

template<class T, int rank>
class FixedArray {
    private:
        T data[rank];
    public:
        T &operator[](int i) { return data[i]; }
        const T &operator[](int i) const { return data[i]; }
    
        FixedArray() {}
        FixedArray(const T& a) {
            for (int i=0; i<rank; ++i) data[i] =a;
        }
        
        FixedArray(const T& a1, const T& a2) {
            assert(rank==2);
            data[0]=a1;
            data[1]=a2;
        }
        
        FixedArray(const T& a1, const T& a2, 
                   const T& a3) {
            assert(rank==3);
            data[0]=a1;
            data[1]=a2;
            data[2]=a3;
        }
        FixedArray(const T& a1, const T& a2, 
                   const T& a3, const T& a4) {
            assert(rank==4);
            data[0]=a1;
            data[1]=a2;
            data[2]=a3;
            data[3]=a4;
        }
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
        
        FixedArray(const FixedArray& A) {
            for (int i=0; i<rank; ++i) data[i] = A[i];
        }
        
        FixedArray &operator=(const FixedArray& A) {
            for (int i=0; i<rank; ++i) data[i] = A[i];
            return *this;
        }

        template<class T2>
        FixedArray operator+(const FixedArray<T2,rank> & A) {
            FixedArray<T,rank> result;
            for (int i=0; i<rank; ++i) result[i] = data[i] + A[i];
            return result;
        }

        template<class T2>
        FixedArray operator-(const FixedArray<T2,rank>& A) {
            FixedArray<T,rank> result;
            for (int i=0; i<rank; ++i) result[i] = data[i] - A[i];
            return result;
        }

        template<class T2>
        FixedArray operator*(const FixedArray<T2,rank>& A) {
            FixedArray<T,rank> result;
            for (int i=0; i<rank; ++i) result[i] = data[i] * A[i];
            return result;
        }

        template<class T2>
        FixedArray operator/(const FixedArray<T2,rank>& A) {
            FixedArray<T,rank> result;
            for (int i=0; i<rank; ++i) result[i] = data[i] / A[i];
            return result;
        }

        FixedArray operator*(T A) {
            FixedArray<T,rank> result;
            for (int i=0; i<rank; ++i) result[i] = data[i]*A;
            return result;
        }

        FixedArray operator/(T A) {
            FixedArray<T,rank> result;
            for (int i=0; i<rank; ++i) result[i] = data[i]/A;
            return result;
        }

        FixedArray operator+(T A) {
            FixedArray<T,rank> result;
            for (int i=0; i<rank; ++i) result[i] = data[i]+A;
            return result;
        }

        FixedArray operator-(T A) {
            FixedArray<T,rank> result;
            for (int i=0; i<rank; ++i) result[i] = data[i]-A;
            return result;
        }

        T product() {
            T Result=1;
            for (int i=0; i<rank; ++i) Result *= data[i];
            return Result;
        }
        
        // This is not safe and should be replaced by something better
        FixedArray(const T *A) {
            for (int i=0; i<rank; ++i) data[i] = A[i];
        }
        FixedArray &operator=(T *A) {
            for (int i=0; i<rank; ++i) data[i] = A[i];
            return *this;
        }
        // This is not very nice and only for compatibility with the 
        // Matrix class
        T *Data() { return data; }
        
        operator T*()  { return data; }
        
        int length() { return rank; }
};


template<class T, int rank>
ostream &operator<<(ostream &O, FixedArray<T,rank> A) {
    O << A[0];
    for (int i=1; i<rank; ++i) O << ", " << A[i];
    return O;
}

template<int rank>
class MatrixIterator {
    private:
        FixedArray<int,rank> from, to;
        FixedArray<int,rank> pos;
    public:
        MatrixIterator( const FixedArray<int,rank> &from_,
                        const FixedArray<int,rank> &to_)
                        : from(from_), to(to_), pos(from) {}
                        
        MatrixIterator &operator++() {
            
            for (int i=rank-1; i > 0; --i) {
                if (++pos[i] < to[i] ) return *this;
                pos[i] = from[i];
            }
            ++pos[0];
            return *this;
        }
        operator bool() {
            return pos[0]<to[0];
        }
        int* operator*() {
            return pos.Data();
        }
        FixedArray<int,rank>& Pos() { return pos; }
};

#endif
