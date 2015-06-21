
#ifndef CSUTIL_BASE
#define CSUTIL_BASE
#include <iostream>
#include <vector>
#include <math.h>
#include <stdexcept>
#include "CsArray.hh"

namespace CsUtil{
  namespace{}
  //array<CT,CI>* TensorProduct(const array<CT,CI>& operand);//future
  
  template <class T> void Normalize(T* begin, T* end);
  template <class T> void Normalize(std::vector< T > inputvector);
  template <class T, class I> void Normalize(array<T,I>* inarray, I axis = -1);

  template <class T> T Sum(T* begin, T* end);
  template <class CT, class CI> array<CT,CI>& Sum(const array<CT,CI>& inarray, CI axis = -1);
  
  template <class T> T Mean(T* xbegin, T* xend);
  template <class CT, class CI> array<CT,CI>& Mean(const array<CT,CI>& inarray, CI axis = -1);
  
  template <class CT, class CI> void Apply(const array<CT,CI>& inarray, void (*f)(CT*,CT*), CI axis = -1);
  template <class CT, class CI> void Apply(const array<CT,CI>& inarray, void (*f)(CT*,CT*,CT*), CT* params, CI axis = -1);
  template <class CT, class CI> array<CT,CI>& RunFunc(const array<CT,CI>& inarray, CT (*f)(CT*,CT*), CI axis = -1);
  template <class CT, class CI> array<CT,CI>& RunFunc(const array<CT,CI>& inarray, CT (*f)(CT*,CT*,CT*), CT* params, CI axis = -1);

  template <class CT, class CI> array<CT,CI>& InnerProduct(const array<CT,CI>& leftarray, const array<CT,CI>& rightarray, CI ndimensions = 1);

  template <class CT, class CI> array<CT,CI>& Transpose(const array<CT,CI>& inarray, CI step = 1);
  
  template<class T> void Histogram(T* xbegin, T* xend, T* hbegin, T* hend, std::vector< T > binedges, bool normed = true);
  template<class T> void HistogramBinaryFill(T* xbegin, T* xend, T* hbegin, T* hend, std::vector< T > binedges, bool normed = true);
    
  template <class T> double Rsquared(int n, T* ybegin, T* fbegin);  

  template <class T> T DifferenceSquared(int n, T* x1begin, T* x2begin);
  template <class T> T DifferenceSquared(T* x1begin, T* x1end, T C);

  namespace{
    template <class CT, class CI> CT RecursiveInnerProduct(const array<CT,CI>& leftarray, const array<CT,CI>& rightarray, array<CI,CI>& leftindex, array<CI,CI>& rightindex, const CI& ndimensions);
    template <class T> std::vector<T> Permute(const std::vector<T>& index, const T& step);
    template <class T> std::vector<T> GeneratePermuteMap(const T& nindices, const T& step);
    template <class CI> array<CI,CI> Permute(const array<CI,CI>& index, const CI& step);//std::vector<CI>
  }

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
    Apply(*inarray, &Normalize, axis);
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template <class T> T Sum(T* begin, T* end){
    T total = 0;
    for (T* i = begin; i <= end; ++i){
      total += *i;
    }
    return total;
  }

  template <class CT, class CI> array<CT,CI>& Sum(const array<CT,CI>& inarray, CI axis){
    return RunFunc(inarray, &Sum, axis);
  }

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  template<class T> void Histogram(T* xbegin, T* xend, T* hbegin, T* hend, std::vector< T > binedges, bool normed = true){
    int n_x = xend - xbegin + 1, n_b = binedges.size();
    for (T* x_i = xbegin; x_i <= xend; ++x_i){
      int bin;
      if (*x_i <= binedges[0]) {bin = 0;}
      else if (*x_i >= binedges[n_b-1]) {bin = n_b-1;}
      else{
	int lowbin = 0, highbin = n_b-1;
	while (true){
	  bin = (highbin*(x_i - *(xbegin+highbin)) - lowbin*(x_i - *(xbegin+lowbin)))/(*(xbegin+highbin) - *(xbegin+lowbin));
	  if (x_i < binedges[bin]){highbin = bin;}
	  else if (x_i > binedges[bin+1]){lowbin = bin+1;}
	  else {break;}
	}
      }
      (*(hbegin + bin))++;
    }

    if (normed){
      for (T* hbin = hbegin; hbin <= hbegin; ++hbin){
	*hbin = *hbin/n_x;
      }
    }
    
  }

  template<class T> void HistogramBinaryFill(T* xbegin, T* xend, T* hbegin, T* hend, std::vector< T > binedges, bool normed = true){
    int n_b = binedges.size();

    for (T* x_i = xbegin; x_i <= xend; ++x_i){
      int bin;
      if (*x_i <= binedges[0]) {bin = 0;}
      else if (*x_i >= binedges[n_b-1]) {bin = n_b-1;}
      else{
	int lowbin = 0, highbin = n_b-1;
	while (true){
	  bin = (highbin+lowbin)/2;
	  if (x_i < binedges[bin]){highbin = bin;}
	  else if (x_i > binedges[bin+1]){lowbin = bin+1;}
	  else {break;}
	}
      }
      (*(hbegin + bin))++;
    }

    if (normed){
      int n_x = xend - xbegin + 1;
      for (T* hbin = hbegin; hbin <= hbegin; ++hbin){
	*hbin = *hbin/n_x;
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
    return RunFunc(inarray, &Mean, axis);
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
  template <class CT, class CI> array<CT,CI>& InnerProduct(const array<CT,CI>& leftarray, const array<CT,CI>& rightarray, CI ndimensions){
    std::vector<CI> newshape(leftarray.NDimensions() + rightarray.NDimensions() - 2*ndimensions);
    
    for (CI i = 0; i < ndimensions; ++i){
      if (leftarray.Shape(leftarray.NDimensions()-1-i) != rightarray.Shape(i)){
	throw std::invalid_argument( "These arrays cannot undergo inner product!" );
      }
    }

    if (newshape.size() > 0){
      for (CI i = 0; i < leftarray.NDimensions()-ndimensions; ++i){newshape[i] = leftarray.Shape(i);}
      for (CI i = 0; i < rightarray.NDimensions()-ndimensions; ++i){newshape[newshape.size()-1-i] = rightarray.Shape(rightarray.NDimensions()-1-i);}
    }else{
      newshape.resize(1);
      newshape[0] = 1;
    }
    array<CT,CI>* newarray = new array<CT,CI>(newshape);
    array<CI,CI> newindex(1,newshape.size()), leftindex(1,leftarray.NDimensions()), rightindex(1,rightarray.NDimensions());
    for (CI index = 0; index < newarray->Size(); ++index){
      newindex = newarray->GetIndex(index);
      for (CI i = 0; i < leftarray.NDimensions()-ndimensions; ++i){leftindex[i] = newindex[i];}
      for (CI i = 0; i < rightarray.NDimensions()-ndimensions; ++i){rightindex[ndimensions+i] = newindex[newarray->NDimensions()-rightarray.NDimensions() + ndimensions+i];}
      newarray->Assign(index, RecursiveInnerProduct(leftarray, rightarray, leftindex, rightindex, ndimensions));
    }

    return *newarray;
  }

  namespace {
    template <class CT, class CI> CT RecursiveInnerProduct(const array<CT,CI>& leftarray, const array<CT,CI>& rightarray, array<CI,CI>& leftindex, array<CI,CI>& rightindex, const CI& ndimensions){
      CT newvalue = 0;
      CI leftactingindex = leftindex.Size()-ndimensions;
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

    //-----------------------------------------------------------------------------
    template <class T>
    std::vector<T> Permute(const std::vector<T>& index, const T& step){
      std::vector<T> newindex(index.size());
      for (T i = 0; i < index.size(); ++i){
	newindex[(i+step)%index.size()] = index[i];
      }
      return newindex;
    }
    
    //-----------------------------------------------------------------------------
    template <class T>
    std::vector<T> GeneratePermuteMap(const T& nindices, const T& step){
      std::vector<T> map(nindices);
      for (T i = 0; i < nindices; ++i){
	map[i] = (i+step)%nindices;
      }
      return map;
    }

        //-----------------------------------------------------------------------------
    template <class CI>
    array<CI,CI> Permute(const array<CI,CI>& index, const CI& step){
      array<CI,CI> newindex(1,index.Size());
      for (CI i = 0; i < index.Size(); ++i){
	newindex[(i+step)%index.Size()] = index[i];
      }
      return newindex;
    }
  }
  
  //-----------------------------------------------------------------------------
  template <class CT, class CI>
  array<CT,CI>& Transpose(const array<CT,CI>& inarray, CI step){
    array<CT,CI>* newarray = new array<CT,CI>(Permute(inarray.Shape(), step));
    for (CI i = 0; i < inarray.Size(); ++i){
      newarray->Assign(newarray->GetIndex(Permute(inarray.GetIndex(i),step)), inarray.GetValue(i));
    }
    return *newarray;
  }

  //-----------------------------------------------------------------------------
  template <class CT, class CI>
  void Apply(const array<CT,CI>& inarray, void (*f)(CT*,CT*), CI axis){
    array<CT*,CI> axispointers = inarray.Axis(axis);
  
    for (CI i = 0; i < axispointers.Shape(0); ++i){
      (*f)(axispointers[2*i], axispointers[2*i + 1]);
    }
  }

  template <class CT, class CI>
  void Apply(const array<CT,CI>& inarray, void (*f)(CT*,CT*,CT*), CT* params, CI axis){
    array<CT*,CI> axispointers = inarray.Axis(axis);

    for (CI i = 0; i < axispointers.Shape(0); ++i){
      (*f)(axispointers[2*i], axispointers[2*i + 1], params);
    }
  }

  template <class CT, class CI>
  array<CT,CI>& RunFunc(const array<CT,CI>& inarray, CT (*f)(CT*,CT*), CI axis){
    if (axis < 0){axis = inarray.NDimensions()-1;}
    array<CT*,CI> axispointers = inarray.Axis(axis);
    std::vector<CI> newshape(axis == 0 ? 1:axis,1);

    for (CI i = 0; i < axis; ++i){newshape[i] = inarray.Shape(i);}
  
    array<CT,CI>* newarray = new array<CT,CI>(newshape);
  
    for (CI i = 0; i < axispointers.Shape(0); ++i){
      newarray->Assign(i, (*f)(axispointers[2*i], axispointers[2*i + 1]));
    }
    return *newarray;
  }

  template <class CT, class CI>
  array<CT,CI>& RunFunc(const array<CT,CI>& inarray, CT (*f)(CT*,CT*,CT*), CT* params, CI axis){
    if (axis == -1){axis = inarray.NDimensions()-1;}
    array<CT*,CI> axispointers = inarray.Axis(axis);
    std::vector<CI> newshape(axis == 0 ? 1:axis,1);

    for (CI i = 0; i < axis; ++i){newshape[i] = inarray.Shape(i);}
  
    array<CT,CI>* newarray = new array<CT,CI>(newshape);
  
    for (CI i = 0; i < axispointers.Shape(0); ++i){
      newarray->Assign(i, (*f)(axispointers[2*i], axispointers[2*i + 1], params));
    }
    return *newarray;
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

#endif
