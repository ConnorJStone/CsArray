
#include <iostream>
#include <vector>
#include <math.h>
#include <stdexcept>
#include "CsArray.hh"


namespace CsUtil{
  namespace{}
  const double PI = 3.14159265358979323;

  template <class T> void Normalize(T* begin, T* end);
  template <class T> void Normalize(std::vector< T > inputvector);
  template <class T, class I> void Normalize(array<T,I>* inarray, I axis = -1);

  template <class T> T Sum(T* begin, T* end);
  template <class CT, class CI> array<CT,CI>& Sum(const array<CT,CI>& inarray, CI axis = -1);
  
  template <class T> T Mean(T* xbegin, T* xend);
  template <class CT, class CI> array<CT,CI>& Mean(const array<CT,CI>& inarray, CI axis = -1);
  
  template <class T> std::vector< T >& Series(const T& end);
  template <class T> std::vector< T >& Series(const T& start, const T& end);
  template <class T> std::vector< T >& Series(const T& start, const T& end, const T& increment);
  template <class I, class T> std::vector< T >& Series(const I& end, T (*a)(I));
  template <class I, class T> std::vector< T >& Series(const I& start, const I& end, T (*a)(I));
  template <class I, class T> std::vector< T >& Series(const I& start, const I& end, const I& increment, T (*a)(I));

  template <class CT, class CI> array<CT,CI>& InnerProduct(const array<CT,CI>& leftarray, const array<CT,CI>& rightarray, CI ndimensions = 1);
  namespace{template <class CT, class CI> CT RecursiveInnerProduct(const array<CT,CI>& leftarray, const array<CT,CI>& rightarray, std::vector<CI>& leftindex, std::vector<CI>& rightindex, const CI& ndimensions);}

  std::vector< double > SphericalVolumePDF(std::vector< double > binedges, double radius = 850.0);

  void Histogram(std::vector< double > x, std::vector< double > binedges, std::vector< double > h, bool normed = true);
  template<class T> void Histogram(T* xbegin, T* xend, T* hbegin, T* hend, std::vector< T > binedges, bool normed = true);
  
  template <class T> double Rsquared(int n, T* ybegin, T* fbegin);  

  template <class T> T DifferenceSquared(int n, T* x1begin, T* x2begin);
  template <class T> T DifferenceSquared(T* x1begin, T* x1end, T C);

  //template <class CT, class CI> array<CT,CI>* MatrixInverse(const array<CT,CI>& operand);
  //template <class CT, class CI> CT MatrixDeterminant(const array<CT,CI>& operand);
  //template <class CT, class CI> array<CT,CI>* MatrixCofactor(const array<CT,CI>& operand);
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> void Normalize(T* begin, T* end){
    T sum = Sum(begin, end);
    for (unsigned int i = 0; i < end - begin + 1; ++i){
      *(begin+i) = *(begin+i)/sum;
    }
  }
  
  template <class T> void Normalize(std::vector< T > inputvector){
    Normalize(&(inputvector.front()), &(inputvector.back()));
  }

  template <class T, class I> void Normalize(array<T,I>* inarray, I axis = -1){
    inarray->Apply(Normalize, axis);
  }
  /**
  template <class T, class I> void Normalize(const array<T,I>& inarray, I axis = -1){
    std::vector<T*> axispointers;
    if (axis == -1){axispointers = inarray.Axis((inarray.Shape()).size()-1);}
    else{axispointers = inarray.Axis(axis);}
    
    for (I i = 0; i < axispointers.size(); ++i){
      Normalize(axispointers[i][0], axispointers[i][1]);
    }
  }
  */
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> T Sum(T* begin, T* end){
    T total = 0;
    for (T* i = begin; i <= end; ++i){
      total += *i;
    }
    return total;
  }

  template <class CT, class CI> array<CT,CI>& Sum(const array<CT,CI>& inarray, CI axis){
    return inarray.RunFunc(Sum, axis);
  }
  /**
  template <class CT, class CI> array<CT,CI>* Sum(const array<CT,CI>& inarray, CI axis){
    std::vector<CT*> axispointers;
    std::vector<CI> newshape, inshape = inarray.Shape();
    if (axis == -1){axispointers = inarray.Axis(inshape.size()-1);}
    else{axispointers = inarray.Axis(axis);}

    if (inshape.size() == 1){
      newshape.resize(1);
      newshape[0] = 1;
    }else{
      newshape.resize(inshape.size()-1);
      newshape.assign(inshape.begin(), inshape.end()-1);
    }
    array<CT,CI>* newarray = new array<CT,CI>(newshape);
    
    for (CI i = 0; i < axispointers.size(); ++i){
      newarray->Assign(i, Sum(axispointers[i][0], axispointers[i][1]));
    }
    return newarray;
  }
  */
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  std::vector< double > SphericalVolumePDF(std::vector< double > binedges, double radius){
    double spherevolume = (4.0/3.0)*PI*pow(radius,3);
    std::vector< double > distribution(binedges.size());
    
    for (int edge = 0; edge < binedges.size()-1; ++edge){
      distribution[edge] = PI*(pow(radius,2)*(binedges[edge+1] - binedges[edge]) - (1.0/3.0)*(pow(binedges[edge+1],3) - pow(binedges[edge],3)))/spherevolume;
    }

    return distribution;
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
	  T x_i = *(xbegin + i);
	  bin = (highbin*(x_i - *(xbegin+highbin)) - lowbin*(x_i - *(xbegin+lowbin)))/(*(xbegin+highbin) - *(xbegin+lowbin));
	  if (x_i < binedges[bin]){highbin = bin;}
	  else if (x_i > binedges[bin+1]){lowbin = bin+1;}
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
  template <class T> double Rsquared(int n, T* ybegin, T* fbegin){   
    return 1.0 - DifferenceSquared(n, ybegin, fbegin)/DifferenceSquared(ybegin, ybegin+n-1, Mean(ybegin, ybegin+n-1));
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> T Mean(T* xbegin, T* xend){
    return Sum(xbegin,xend)/(xend - xbegin + 1);
  }

  template <class CT, class CI> array<CT,CI>& Mean(const array<CT,CI>& inarray, CI axis){
    return inarray.RunFunc(Mean, axis);
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> T DifferenceSquared(int n, T* x1begin, T* x2begin){
    T SS = 0;
    for (int i = 0; i < n; ++i){
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
  template <class T> std::vector< T >& Series(const T& end){
    return Series((T) 0, end, (T) 1);
  }

  template <class T> std::vector< T >& Series(const T& start, const T& end){
    return Series(start, end, (T) 1);
  }

  template <class T> std::vector< T >& Series(const T& start, const T& end, const T& increment){
    if (end < start){return Series(end, start, increment);}
    else if (end == start){throw;}

    std::vector< T >* series = new std::vector< T >;
    T a = start;
    while (a < end){
      series->push_back(a);
      a += increment;
    }
    return *series;
  }

  template <class I, class T> std::vector< T >& Series(const I& end, T (*a)(I)){
    return Series((I) 0, end, (I) 1, a );
  }

  template <class I, class T> std::vector< T >& Series(const I& start, const I& end, T (*a)(I)){
    return Series(start, end, (I) 1, a );
  }
  
  template <class I, class T> std::vector< T >& Series(const I& start, const I& end, const I& increment, T (*a)(I)){
    if (end < start){return Series(end, start, increment, a);}
    else if (end == start){throw;}
    
    std::vector< T >* series = new std::vector< T >;
    I i = start;
    while (i < end){
      series->push_back((*a)(i));
      i += increment;
    }
    return *series;
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class CT, class CI> array<CT,CI>& InnerProduct(const array<CT,CI>& leftarray, const array<CT,CI>& rightarray, CI ndimensions){
    const std::vector<CI> leftshape = leftarray.Shape(), rightshape = rightarray.Shape();
    std::vector<CI> newshape(leftshape.size() + rightshape.size() - 2*ndimensions);
    
    for (CI i = 0; i < ndimensions; ++i){
      if (leftshape[leftshape.size()-1-i] != rightshape[i]){
	throw std::invalid_argument( "These arrays cannot undergo inner product!" );
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
      newarray->Assign(index, RecursiveInnerProduct(leftarray, rightarray, leftindex, rightindex, ndimensions));
    }

    return *newarray;
  }

  namespace {
    template <class CT, class CI> CT RecursiveInnerProduct(const array<CT,CI>& leftarray, const array<CT,CI>& rightarray, std::vector<CI>& leftindex, std::vector<CI>& rightindex, const CI& ndimensions){
      CT newvalue = 0;
      CI leftactingindex = leftindex.size()-ndimensions;
      CI actingindexlength = leftarray.Shape(leftactingindex);
      if (ndimensions == 1){
	CT* leftbegin = leftarray.GetValueP(leftindex), *rightbegin = rightarray.GetValueP(rightindex);
	CI columnstep = rightarray.Size(1);
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

  /**
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class CT, class CI> array<CT,CI>* MatrixInverse(const array<CT,CI>& operand){
    //if (operand.Shape().size() != 2){throw std::invalid_argument( "This is not a matrix" );}
    //else if (operand.Shape()[0] != operand.Shape()[1]){throw std::invalid_argument( "This matrix is not square" );}
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class CT, class CI> CT MatrixDeterminant(const array<CT,CI>& operand){}

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class CT, class CI> array<CT,CI>* MatrixCofactor(const array<CT,CI>& operand){}
  */
}
