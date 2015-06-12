
//#include <iostream>
#include <vector>
#include <math.h>
#include "CsArray.hh"
//#include <cstdarg.h>

namespace CsUtil{
  const double PI = 3.14159265358979323;

  template <class T> void Normalize(T* begin, T* end);
  template <class T> void Normalize(std::vector< T > inputvector);
  template <class T, class I> void Normalize(array<T,I> inarray, I axis = -1);

  void SphericalVolumePDF(std::vector< double > binedges, std::vector< double > distribution, double radius = 850.0);

  void Histogram(std::vector< double > x, std::vector< double > binedges, std::vector< double > h, bool normed = true);
  template<class T> void Histogram(T* xbegin, T* xend, T* hbegin, T* hend, std::vector< T > binedges, bool normed = true);

  template <class T> double Rsquared(T* ybegin, T* yend, T* fbegin, T* fend);  

  template <class T> T Mean(T* xbegin, T* xend);

  template <class T> T DifferenceSquared(T* x1begin, T* x1end, T* x2begin, T* x2end);
  template <class T> T DifferenceSquared(T* x1begin, T* x1end, T C);

  template <class T> std::vector< T > Series(const T& end);
  template <class T> std::vector< T > Series(const T& start, const T& end);
  template <class T> std::vector< T > Series(const T& start, const T& end, const T& increment);

  template <class CT, class CI> array<CT,CI>* InnerProduct(const array<CT,CI>& leftarray, const array<CT,CI>& rightarray, CI ndimensions = 1);
  template <class CT, class CI> CT RecursiveInnerProduct(const array<CT,CI>* leftarray, const array<CT,CI>* rightarray, std::vector<CI>& leftindex, std::vector<CI>& rightindex, CI ndimensions);
  
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> void Normalize(T* begin, T* end){
    T sum = accumulate(begin, end, (T) 0 );
    for (unsigned int i = 0; i < end - begin + 1; ++i){
      *(begin+i) = *(begin+i)/sum;
    }
  }
  
  template <class T> void Normalize(std::vector< T > inputvector){
    Normalize(&(inputvector.front()), &(inputvector.back()));
  }

  template <class T, class I> void Normalize(array<T,I> inarray, I axis = -1){
    std::vector<T*> axispointers;
    if (axis == -1){axispointers = inarray.Axis((inarray.Shape()).size()-1);}
    else{axispointers = inarray.Axis(axis);}
    
    for (I i = 0; i < axispointers.size(); ++i){
      Normalize(axispointers[i][0], axispointers[i][1]);
    }
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  void SphericalVolumePDF(std::vector< double > binedges, std::vector< double > distribution, double radius){
    double spherevolume = (4.0/3.0)*PI*pow(radius,3);
  
    for (int edge = 0; edge < binedges.size()-1; ++edge){
      distribution[edge] = PI*(pow(radius,2)*(binedges[edge+1] - binedges[edge]) - (1.0/3.0)*(pow(binedges[edge+1],3) - pow(binedges[edge],3)))/spherevolume;
    }
  
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  void Histogram(std::vector< double > x, std::vector< double > binedges, std::vector< double > h, bool normed){
    int n_x = x.size(), n_b = binedges.size();
    /**
       y = mx +b
       y1 = mx1+b
       y2 = mx2+b

       y = (y2(x - x1) - y1(x - x2))/(x2 - x1)
     */
    for (int i = 0; i < n_x; ++i){
      int bin;
      if (x[i] <= binedges[0]) {bin = 0;}
      else if (x[i] >= binedges[n_b-1]) {bin = n_b-1;}
      else{
	int lowbin = 0, highbin = n_b-1;
	while (true){
	  bin = (highbin*(x[i] - x[highbin]) - lowbin*(x[i] - x[lowbin]))/(x[highbin] - x[lowbin]);
	  if (x[i] < binedges[bin]){highbin = bin;}
	  else if (x[i] > binedges[bin+1]){lowbin = bin+1;}
	  else {break;}
	}
      }
      h[bin]++;
    }

    if (normed){
      for (int bin = 0; bin < h.size(); ++bin){
	h[bin] = h[bin]/n_x;
      }
    }
    
  }

  template<class T> void Histogram(T* xbegin, T* xend, T* hbegin, T* hend, std::vector< T > binedges, bool normed = true){
    int n_x = xend - xbegin + 1, n_b = binedges.size();
    /**
       y = mx +b
       y1 = mx1+b
       y2 = mx2+b

       rearrange:
       y = (y2(x - x1) - y1(x - x2))/(x2 - x1)
     */
    for (int i = 0; i < n_x; ++i){
      int bin;
      if (*(xbegin + i) <= binedges[0]) {bin = 0;}
      else if (*(xbegin + i) >= binedges[n_b-1]) {bin = n_b-1;}
      else{
	int lowbin = 0, highbin = n_b-1;
	while (true){
	  bin = (highbin*(*(xbegin + i) - *(xbegin+highbin)) - lowbin*(*(xbegin + i) - *(xbegin+lowbin)))/(*(xbegin+highbin) - *(xbegin+lowbin));
	  if (*(xbegin + i) < binedges[bin]){highbin = bin;}
	  else if (*(xbegin + i) > binedges[bin+1]){lowbin = bin+1;}
	  else {break;}
	}
      }
      *(hbegin + bin)++;
    }

    if (normed){
      for (int bin = 0; bin < hend-hbegin+1; ++bin){
	*(hbegin + bin) = *(hbegin + bin)/n_x;
      }
    }
    
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> double Rsquared(T* ybegin, T* yend, T* fbegin, T* fend){   
    return 1.0 - DifferenceSquared(ybegin, yend, fbegin, fend)/DifferenceSquared(ybegin, yend, Mean(ybegin, yend));
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> T Mean(T* xbegin, T* xend){
    return accumulate(xbegin,xend,0.0)/(xend - xbegin + 1);
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> T DifferenceSquared(T* x1begin, T* x1end, T* x2begin, T* x2end){
    T SS = 0;
    for (int i = 0; i < x1end - x1begin + 1; ++i){
      SS += pow(*(x1begin + i) - *(x2begin + i), 2);
    }
    return SS;
  }

  template <class T> T DifferenceSquared(T* x1begin, T* x1end, T C){
    T SS = 0;
    for (int i = 0; i < x1end - x1begin + 1; ++i){
      SS += pow(*(x1begin + i) - C, 2);
    }
    return SS;
  }

  
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> std::vector< T > Series(const T& end){
    return Series((T) 0, end, (T) 1);
  }

  template <class T> std::vector< T > Series(const T& start, const T& end){
    return Series(start, end, (T) 1);
  }

  template <class T> std::vector< T > Series(const T& start, const T& end, const T& increment){
    if (end < start){return Series(end, start, increment);}
    else if (end == start){throw;}

    std::vector< T > series;
    T a = start;
    while (a < end){
      series.push_back(a);
      a += increment;
    }
    return series;
  }
  /**
  template <class I, class T> std::vector< T > Series(I end, T a(I)){
    return Series((I) 0, end, (I) 1, a );
  }

  template <class I, class T> std::vector< T > Series(I start, I end, T a(I)){
    return Series(start, end, (I) 1, a );
  }
  
  template <class I, class T> std::vector< T > Series(I start, I end, I increment, T a(I)){
    if (end < start){return Series(end, start, increment, a);}
    else if (end == start){throw;}
    
    std::vector< T > series;
    I i = start;
    while (i < end){
      series.push_back(a(i));
      i += increment;
    }
    return series;
  }
  */
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class CT, class CI> array<CT,CI>* InnerProduct(const array<CT,CI>& leftarray, const array<CT,CI>& rightarray, CI ndimensions){
    const std::vector<CI> leftshape = leftarray.Shape(), rightshape = rightarray.Shape();
    std::vector<CI> newshape(leftshape.size() + rightshape.size() - 2*ndimensions);
    for (CI i = 0; i < ndimensions; ++i){
      if (leftshape[leftshape.size()-1-i] != rightshape[i]){
      }
    }

    if (leftshape.size() + rightshape.size() - 2*ndimensions > 0){
      for (CI i = 0; i < leftshape.size()-ndimensions; ++i){newshape[i] = leftshape[i];}
      for (CI i = 0; i < rightshape.size()-ndimensions; ++i){newshape[newshape.size()-1-i] = rightshape[rightshape.size()-1-i];}
    }else{
      newshape.resize(1);
      newshape[0] = 1;
    }
    array<CT,CI>* newarray = new array<CT,CI>(newshape);

    std::vector<CI> newindex(newshape.size(),0), leftindex(leftshape.size(),0), rightindex(rightshape.size(),0);
  
    for (CI index = 0; index < newarray->Size(); ++index){
      newindex = newarray->GetIndex(index);
      for (CI i = 0; i < leftshape.size()-ndimensions; ++i){leftindex[i] = newindex[i];}
      for (CI i = 0; i < rightshape.size()-ndimensions; ++i){rightindex[ndimensions+i] = newindex[newindex.size()-rightshape.size() + ndimensions+i];}
      newarray->Assign(index, RecursiveInnerProduct(&leftarray, &rightarray, leftindex, rightindex, ndimensions));
    }

    return newarray;
  }

  template <class CT, class CI> CT RecursiveInnerProduct(const array<CT,CI>* leftarray, const array<CT,CI>* rightarray, std::vector<CI>& leftindex, std::vector<CI>& rightindex, CI ndimensions){
    CT newvalue = 0;
    CI leftactingindex = leftindex.size()-ndimensions;
    CI actingindexlength = leftarray->Shape(leftactingindex);
    if (ndimensions == 1){
      //leftindex[leftactingindex] = 0;
      //rightindex[0] = 0;
      CT* leftbegin = leftarray->GetValueP(leftindex), *rightbegin = rightarray->GetValueP(rightindex);
      CI columnstep = rightarray->Size(1);
      for (CI index = 0; index < actingindexlength; ++index){
	newvalue += (*(leftbegin + index))*(*(rightbegin + index*columnstep));
      }
    }else{
      for (CI index = 0; index < actingindexlength; ++index){
	leftindex[leftactingindex] = index;
	rightindex[ndimensions-1] = index;
	newvalue += RecursiveInnerProduct(leftarray, rightarray, leftindex, rightindex, ndimensions-1);
      }
    }
    return newvalue;
  }

}
