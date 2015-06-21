/**
   Author: Connor Stone
   Description: The "Series" templates are for producing arrays of values with simplistic spacing. They may also be given a function as an argument which determines the spacing arbitrarily. 
 */

#ifndef CSUTIL_SERIES
#define CSUTIL_SERIES
#include <vector>
#include <stdexcept>

namespace CsUtil{
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> array<T,int>& Series(const T& end);
  template <class T> array<T,int>& Series(const T& start, const T& end);
  template <class T> array<T,int>& Series(const T& start, const T& end, const T& increment);
  template <class I, class T> array<T,I>& Series(const I& n, T (*a)(I));
  template <class I, class T> array<T,I>& Series(const T& start, const I& n, T (*a)(I));
  template <class I, class T> array<T,I>& Series(const T& start, const T& increment, const I& n, T (*a)(I));
  
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> array<T,int>& Series(const T& end){
    return Series((T) 0, end, (T) 1);
  }

  template <class T> array<T,int>& Series(const T& start, const T& end){
    return Series(start, end, (T) 1);
  }

  template <class T> array<T,int>& Series(const T& start, const T& end, const T& increment){
    if (end < start){return Series(end, start, increment);}
    else if (end == start){throw std::invalid_argument( "The inputs cannot be equivalent" );}

    array<T,int>* series = new array<T,int>(1,(end-start)/increment);
    T a = start;
    for (int count = 0; count < series->Size(); ++count){
      series->Assign(count, a);
      a += increment;
    }
    return *series;
  }

  template <class I, class T> array<T,I>& Series(const I& n, T (*a)(I)){
    return Series((T) 0, (T) 1, n, a );
  }

  template <class I, class T> array<T,I>& Series(const T& start, const I& n, T (*a)(I)){
    return Series(start, (T) 1, n, a );
  }
  
  template <class I, class T> array<T,I>& Series(const T& start, const T& increment, const I& n, T (*a)(I)){
    if (n < 1){throw std::invalid_argument( "The series must have at least one element" );}
    
    array<T,I>* series = new array<T,I>(1,n);
    I i = start;
    for (I count = 0; count < n; ++count){
      series->Assign(count, (*a)(i));
      i += increment;
    }
    return *series;
  }


}
#endif
