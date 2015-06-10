
#include <vector>
#include <iostream>

//-----------------------------------------------------------------------------
template <class CT, class CI>
class array {
  //--------------------------------------------------
  std::vector<CI> shape;
  CT* origin;
  std::vector<CI> X;
  
  void FillX();
  
  std::vector<CI> GetLocation(CI i);
  
  CI GetLocation(std::vector<CI> index);
  
  bool InRange(CI location, std::vector< std::vector<CI> > indices);
  
  bool InRange(std::vector<CI> location, std::vector< std::vector<CI> > indices);
  
  std::vector<CI> Permute(std::vector<CI> index, CI step);
  
  CI Permute(CI index, CI step);

public:

  array();

  array(std::vector<CI> D);

  ~array();

  CT GetElement(std::vector<CI> index);

  CT GetElement(CI index);

  //std::vector< std::vector<CT*> > GetSubArray(std::vector< std::vector<CI> > indices);

  //array<CT,CI> SubArray(std::vector< std::vector<CI> > indices);

  CI Size(CI dimension = 0);

  std::vector< std::vector<CT*> > Axis(CI dimension);

  CT* Begin();

  std::vector<CT*> Begin(CI dimension);

  CT* End();

  std::vector<CT*> End(CI dimension);

  void Assign(CT* selfstartat, CT* databegin, CT* dataend);

  void Assign(CI selfstartatindex, CT* databegin, CT* dataend);

  void Assign(CT* databegin, CT* dataend);

  void Assign(std::vector<CI> index, CT value);

  void Assign(CI index, CT value);

  array<CT,CI> Transpose(CI step = 1);

  //string Print();

  //array<CT,CI> TensorContract(array<CT,CI> operand);

  //array<CT,CI> TensorProduct(array<CT,CI> operand);//future

  /**
  array<CT,CI> operator*(array<CT,CI> operand);

  array<CT,CI> operator*(CT operand);

  void operator*=(array<CT,CI> operand);

  void operator*=(CT operand);

  array<CT,CI> operator+(array<CT,CI> operand);

  array<CT,CI> operator+(CT operand);

  void operator+=(array<CT,CI> operand);

  void operator+=(CT operand);

  array<CT,CI> operator-(array<CT,CI> operand);

  array<CT,CI> operator-(CT operand);

  void operator-=(array<CT,CI> operand);

  void operator-=(CT operand);

  array<CT,CI> operator/(array<CT,CI> operand);

  array<CT,CI> operator/(CT operand);

  void operator/=(array<CT,CI> operand);

  void operator/=(CT operand);

  array<CT,CI> operator[](CI index);
  */
};


//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::FillX(){
  X.resize(shape.size()+1);
  for (CI a = 0; a < X.size(); ++a){
    CI distance = 1;
    for (CI d = a; d < shape.size(); ++d){distance *= shape[d];}
    X[a] = distance;
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector<CI> array<CT, CI>::GetLocation(CI i){
  std::vector<CI> location(shape.size());
  for (CI a = 0; a < location.size(); ++a){
    location[a] = (i%X[a])/X[a+1];
  }
  return location;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CI array<CT, CI>::GetLocation(std::vector<CI> index){
  CI location = 0;
  for (CI i = 0; i < index.size(); ++i){
    location += index[i]*X[i+1];
  }
  return location;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
bool array<CT, CI>::InRange(CI location, std::vector< std::vector<CI> > indices){
  return InRange(GetLocation(location), indices);
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
bool array<CT, CI>::InRange(std::vector<CI> location, std::vector< std::vector<CI> > indices){
  for (CI i = 0; i < location.size(); ++i){
    if (location[i] < indices[i][0] || location[i] > indices[i][1]){return false;}
  }
  return true;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector<CI> array<CT, CI>::Permute(std::vector<CI> index, CI step){
  std::vector<CI> newindex(index.size());
  for (CI i = 0; i < index.size(); ++i){
    newindex[(i+step)%index.size()] = index[i];
  }
  return newindex;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CI array<CT, CI>::Permute(CI index, CI step){
  return GetLocation(Permute(GetLocation(index), step));
}



//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT, CI>::array(){
  std::cout << "array constructed" << std::endl;
}


//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT, CI>::array(std::vector<CI> D){
  shape.resize(D.size());
  shape.assign(D.begin(), D.end());
  FillX();
  origin = new CT[X[0]];
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT, CI>::~array(){delete[] origin;}


//-----------------------------------------------------------------------------
template <class CT, class CI>
CT array<CT, CI>::GetElement(std::vector<CI> index){
  return *(origin + GetLocation(index));
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CT array<CT, CI>::GetElement(CI index){
  return *(origin + index);
}
    
//-----------------------------------------------------------------------------
/**
template <class CT, class CI>
std::vector< std::vector<CT*> > array<CT, CI>::GetSubArray(std::vector< std::vector<CI> > indices){
  std::vector< std::vector<CT*> > subarray;
  std::vector<CI> startat(indices.size()),endat(indices.size());
  for (CI i = 0; i < startat.size(); ++i){
    startat[i] = indices[i][0];
    endat[i] = indices[i][1];
  }
  CI start = GetLocation(startat), end = GetLocation(endat), basestep = indices[indices.size()-1][1] - indices[indices.size()-1][0];
      
  for (CI i = start; i <= end; ++i){
    if (InRange(i,indices)){
      subarray.push_back(std::vector<CT*> (2,origin+i));
      i += basestep;
      subarray[subarray.size()-1][1] = origin + i;
    }
  }
  return subarray;
}
*/
//-----------------------------------------------------------------------------
/**
template <class CT, class CI>
array<CT,CI> array<CT, CI>::SubArray(std::vector< std::vector<CI> > indices){
  std::vector<CI> subshape(indices.size());
  for (CI i = 0; i < subshape.size(); ++i){
    subshape[i] = indices[i][1] - indices[i][0] + 1;
  }
  array<CT,CI> subarray(subshape);
      
  return subarray;
}
*/
//-----------------------------------------------------------------------------
template <class CT, class CI>
CI array<CT, CI>::Size(CI dimension){
  return X[dimension];
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector< std::vector<CT*> > array<CT, CI>::Axis(CI dimension){
  std::vector< std::vector<CT*> > axis(X[0]/X[dimension]);
  for (CI i = 0; i < axis.size(); ++i){
    axis[i].resize(2);
    axis[i][0] = origin + i*X[dimension];
    axis[i][1] = origin + (i+1)*X[dimension] - 1;
  }
  return axis;
}
    

//-----------------------------------------------------------------------------
template <class CT, class CI>
CT* array<CT, CI>::Begin(){return origin;}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector<CT*> array<CT, CI>::Begin(CI dimension){
  std::vector<CT*> beginnings(X[0]/X[dimension]);
  for (CI i = 0; i < beginnings.size(); ++i){
    beginnings[i] = origin+i*X[dimension];
  }
  return beginnings;
}
    

//-----------------------------------------------------------------------------
template <class CT, class CI>
CT* array<CT, CI>::End(){return origin+X[0]-1;}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector<CT*> array<CT, CI>::End(CI dimension){
  std::vector<CT*> endings(X[0]/X[dimension]);
  for (CI i = 0; i < endings.size(); ++i){
    endings[i] = origin+(i+1)*X[dimension] - 1;
  }
  return endings;
}
    

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(CT* selfstartat, CT* databegin, CT* dataend){
  for (CI i = 0; i <= dataend-databegin; ++i){
    *(selfstartat+i) = *(databegin + i);
  }
}
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(CI selfstartatindex, CT* databegin, CT* dataend){
  for (CI i = 0; i <= dataend-databegin; ++i){
    *(origin+selfstartatindex+i) = *(databegin + i);
  }
}
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(CT* databegin, CT* dataend){
  for (CI i = 0; i <= dataend-databegin; ++i){
    *(origin+i) = *(databegin + i);
  }
}
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(std::vector<CI> index, CT value){
  *(origin + GetLocation(index)) = value;
}
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(CI index, CT value){
  *(origin + index) = value;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI> array<CT, CI>::Transpose(CI step){
  array<CT,CI> newarray(Permute(shape, step));
  for (CI i = 0; i < X[0]; ++i){
    newarray.Assign(Permute(i,step), GetElement(i));
  }
  return newarray;
}

//-----------------------------------------------------------------------------
/**
template <class CT, class CI>
std::string array<CT, CI>::Print(){
  std::string out = "";
  return out;
}
*/
//-----------------------------------------------------------------------------
/**
template <class CT, class CI>
array<CT,CI> array<CT, CI>::TensorContract(array<CT,CI> operand){
  return *this;
}
*/

/**
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI> array<CT, CI>::operator*(array<CT,CI> operand){
  array<CT,CI> product(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand.Begin(), *productbegin = product.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(productbegin + i) = (*(selfbegin + i))*(*(operandbegin+i));
  }
  return product;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI> array<CT, CI>::operator*(CT operand){
  array<CT,CI> product(shape);
      
  CT* selfbegin = Begin(), *productbegin = product.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(productbegin + i) = (*(selfbegin + i))*operand;
  }
  return product;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI> array<CT, CI>::operator+(array<CT,CI> operand){
  array<CT,CI> addition(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand.Begin(), *additionbegin = addition.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(additionbegin + i) = (*(selfbegin + i))+(*(operandbegin+i));
  }
  return addition;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI> array<CT, CI>::operator+(CT operand){
  array<CT,CI> addition(shape);
      
  CT* selfbegin = Begin(), *additionbegin = addition.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(additionbegin + i) = (*(selfbegin + i))+operand;
  }
  return addition;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI> array<CT, CI>::operator-(array<CT,CI> operand){
  array<CT,CI> difference(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand.Begin(), *differencebegin = difference.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(differencebegin + i) = (*(selfbegin + i)) - (*(operandbegin+i));
  }
  return difference;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI> array<CT, CI>::operator-(CT operand){
  array<CT,CI> difference(shape);
      
  CT* selfbegin = Begin(), *differencebegin = difference.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(differencebegin + i) = (*(selfbegin + i)) - operand;
  }
  return difference;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI> array<CT, CI>::operator/(array<CT,CI> operand){
  array<CT,CI> quotient(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand.Begin(), *quotientbegin = quotient.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(quotientbegin + i) = (*(selfbegin + i))/(*(operandbegin+i));
  }
  return quotient;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI> array<CT, CI>::operator/(CT operand){
  array<CT,CI> quotient(shape);
      
  CT* selfbegin = Begin(), *quotientbegin = quotient.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(quotientbegin + i) = (*(selfbegin + i))/operand;
  }
  return quotient;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator*=(array<CT,CI> operand){
  CT* selfbegin = Begin(), *operandbegin = operand.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))*(*(operandbegin+i));
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator*=(CT operand){
  CT* selfbegin = Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))*operand;
  }
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator+=(array<CT,CI> operand){
  CT* selfbegin = Begin(), *operandbegin = operand.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))+(*(operandbegin+i));
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator+=(CT operand){
  CT* selfbegin = Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))+operand;
  }
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator-=(array<CT,CI> operand){
  CT* selfbegin = Begin(), *operandbegin = operand.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))-(*(operandbegin+i));
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator-=(CT operand){
  CT* selfbegin = Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))-operand;
  }
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator/=(array<CT,CI> operand){
  CT* selfbegin = Begin(), *operandbegin = operand.Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))/(*(operandbegin+i));
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator/=(CT operand){
  CT* selfbegin = Begin();
  for (CI i = 0; i < X[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))/operand;
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI> array<CT, CI>::operator[](CI index){
  std::vector<CI> newshape(shape.size()-1);
  newshape.assign(shape.begin()+1,shape.end());
  array<CT,CI> newarray(newshape);
  newarray.Assign(origin+index*X[1], origin+(index+1)*X[1]-1);
  return newarray;
}
*/
