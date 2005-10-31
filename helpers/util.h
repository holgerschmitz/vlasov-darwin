#include <math.h>
#ifndef UTIL_H
#define UTIL_H

#if ! defined (__GNUG__) && ! defined (__HP_aCC)
template<class T>
inline T max(T a, T b)
{
  return a > b ? a:b;
}

template<class T>
inline T min(T a, T b)
{
  return a < b ? a:b;
}

template<class T>
inline T abs(T a)
{
  return (a > 0) ? a : -a;
}
#endif

template<class T>
inline T sqr(T a)
{
  return a*a;
}

inline double norm(double a)
{
  return a*a;
}

inline float norm(float a)
{
  return a*a;
}

template<class T>
inline T sign(T a, T b)
{
  return (b >= 0) ? abs(a) : -abs(a);
}

template<class T>
inline T log2(T a)
{ 
  const double inv_log2 = 1/log(2.);
  return T(inv_log2*log(a));
}

#endif
