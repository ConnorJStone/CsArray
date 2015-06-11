#include <stdarg.h>
#include <vector>
#include <iostream>

template <class CT, class CI>
class array {
  //Contains the lengths of each dimension
  std::vector<CI> shape;
  //Pointer to the beginning of the array data
  CT* origin;
  //The lengths in memory of each dimension. offsets[0] is the length of the whole array, offsets[1] is the length along the first dimension,...,offsets[n+1] is always 1 (to represent the length of one data object) 
  std::vector<CI> offsets;

  //Given the requested dimensions, it fills the shape and offsets objects, also creates the array in memory and assigns it to the origin
  void InitializeSelf(std::vector<CI> D);

  //Given an index along the memory, it returns a vector which has the distances along each dimension
  std::vector<CI> GetIndex(CI index);

  //Given the distances along each dimension, this will return the true index of the desired data point
  CI GetIndex(std::vector<CI> index);

  //bool InRange(CI location, std::vector< std::vector<CI> > indices);
  
  //bool InRange(std::vector<CI> location, std::vector< std::vector<CI> > indices);

  //Given a vector of indices, this will return a new vector with the indices permuted by an amount equal to step
  std::vector<CI> Permute(std::vector<CI> index, CI step);

  //Given a index in the memory, this will return the permuted index in memory
  CI Permute(CI index, CI step);

  //Creates a vector where accessing an element returns the permuted index
  std::vector<CI> GeneratePermuteMap(CI nindices, CI step);

  //Recursively applies the calculations needed for n dimensional inner product
  CT RecursiveInnerProduct(array<CT,CI>* operand, std::vector<CI> THISINDEX, std::vector<CI> OPERANDINDEX, CI ndimensions);
    
public:

  //Empty initializer
  array();

  //Initializes the array with no guarintee about the values in its indices
  array(std::vector<CI> D);

  //Initializes the array and assigns a value to every position
  //array(std::vector<CI> D, CT value);

  //Destructor, will call delete on origin
  ~array();

  //Returns the value at an index location
  CT* GetValue(std::vector<CI> index);

  //Returns the value at the index
  CT* GetValue(CI index);

  //Returns the size of the array along the given dimension. note1: size(0) gives the size of the whole array in memory, note2: size(n) is always 1
  CI Size(CI dimension = 0);

  //Returns the vector of values describing how long each dimension is
  std::vector<CI> Shape();

  //Returns the length of one of the dimensions, as specified
  CI Shape(CI dimension);

  //Returns the beginning and end pointer for each segment of the array according to the dimension requested 
  std::vector< std::vector<CT*> > Axis(CI dimension);

  //Returns a pointer to the beginning of the array
  CT* Begin();

  //Returns a vector of pointers to the beginning of each segment on a specified dimension
  std::vector<CT*> Begin(CI dimension);

  //Returns a pointer to the end of the array
  CT* End();

  //Returns a vector of pointers to the end of each segment on a specified dimension
  std::vector<CT*> End(CI dimension);

  //Assign the values between the given begin/end pointers starting at memory index: origin + selfstartatindex
  void Assign(CI selfstartatindex, CT* databegin, CT* dataend);

  //Assign the values at the pointed to locations to the array
  void Assign(CT* databegin, CT* dataend);

  //Assign the given value to an index
  void Assign(std::vector<CI> index, CT value);

  //Assign the given value to an index
  void Assign(CI index, CT value);

  //Assign the given value to every element in the array
  void Assign(CT value);

  //Assign a list of values to the array, the first argument is the number of values to assign
  //  void Assign(CI n, ...);

  //Returns a new array which is the transpose of the old one by a given number of steps (if you are using matrices, the default is what you want)
  array<CT,CI>* Transpose(CI step = 1);

  //string Print();

  //Computes the inner product of two arrays, for matrices this is known as matrix multiplication, it scales to n dimensions 
  array<CT,CI>* InnerProduct(array<CT,CI>* operand, CI ncontractionindices = 1);

  //array<CT,CI> TensorProduct(array<CT,CI> operand);//future
  
  array<CT,CI>* operator*(array<CT,CI>* operand);

  array<CT,CI>* operator*(CT operand);
  
  void operator*=(array<CT,CI>* operand);
  
  void operator*=(CT operand);

  array<CT,CI>* operator+(array<CT,CI>* operand);

  array<CT,CI>* operator+(CT operand);

  void operator+=(array<CT,CI>* operand);

  void operator+=(CT operand);

  array<CT,CI>* operator-(array<CT,CI>* operand);

  array<CT,CI>* operator-(CT operand);

  void operator-=(array<CT,CI>* operand);

  void operator-=(CT operand);

  array<CT,CI>* operator/(array<CT,CI>* operand);

  array<CT,CI>* operator/(CT operand);

  void operator/=(array<CT,CI>* operand);

  void operator/=(CT operand);

  array<CT,CI>* operator[](CI index);

};


//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::InitializeSelf(std::vector<CI> D){
  shape.resize(D.size());
  shape.assign(D.begin(), D.end());
  offsets.resize(shape.size()+1);
  for (CI a = 0; a < offsets.size(); ++a){
    CI axislength = 1;
    for (CI d = a; d < shape.size(); ++d){axislength *= shape[d];}
    offsets[a] = axislength;
  }
  origin = new CT[offsets[0]];
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector<CI> array<CT, CI>::GetIndex(CI index){
  std::vector<CI> newindex(shape.size());
  for (CI i = 0; i < newindex.size(); ++i){
    newindex[i] = (index%offsets[i])/offsets[i+1];
  }
  return newindex;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CI array<CT, CI>::GetIndex(std::vector<CI> index){
  CI newindex = 0;
  for (CI i = 0; i < index.size(); ++i){
    newindex += index[i]*offsets[i+1];
  }
  return newindex;
}

/**
//-----------------------------------------------------------------------------
template <class CT, class CI>
bool array<CT, CI>::InRange(CI location, std::vector< std::vector<CI> > indices){
  return InRange(GetIndex(location), indices);
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
bool array<CT, CI>::InRange(std::vector<CI> location, std::vector< std::vector<CI> > indices){
  for (CI i = 0; i < location.size(); ++i){
    if (location[i] < indices[i][0] || location[i] > indices[i][1]){return false;}
  }
  return true;
}
*/
  
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
  return GetIndex(Permute(GetIndex(index), step));
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector<CI> array<CT, CI>::GeneratePermuteMap(CI nindices, CI step){
  std::vector<CI> map(nindices);
  for (CI i = 0; i < nindices; ++i){
    map[i] = (i+step)%nindices;
  }
  return map;
}


//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT, CI>::array(){
  //std::cout << "array constructed" << std::endl;
}


//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT, CI>::array(std::vector<CI> D){
  InitializeSelf(D);
}

/**
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT, CI>::array(std::vector<CI> D, CT value){
  InitializeSelf(D);
  Assign(value);
}
*/
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT, CI>::~array(){delete[] origin;}


//-----------------------------------------------------------------------------
template <class CT, class CI>
CT* array<CT, CI>::GetValue(std::vector<CI> index){
  return (origin + GetIndex(index));
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CT* array<CT, CI>::GetValue(CI index){
  return (origin + index);
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CI array<CT, CI>::Size(CI dimension){
  return offsets[dimension];
}
  
//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector<CI> array<CT, CI>::Shape(){
  return shape;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CI array<CT, CI>::Shape(CI dimension){
  return shape[dimension];
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector< std::vector<CT*> > array<CT, CI>::Axis(CI dimension){
  std::vector< std::vector<CT*> > axis(offsets[0]/offsets[dimension]);
  for (CI i = 0; i < axis.size(); ++i){
    axis[i].resize(2);
    axis[i][0] = origin + i*offsets[dimension];
    axis[i][1] = origin + (i+1)*offsets[dimension] - 1;
  }
  return axis;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
CT* array<CT, CI>::Begin(){return origin;}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector<CT*> array<CT, CI>::Begin(CI dimension){
  std::vector<CT*> beginnings(offsets[0]/offsets[dimension]);
  for (CI i = 0; i < beginnings.size(); ++i){
    beginnings[i] = origin+i*offsets[dimension];
  }
  return beginnings;
}
    

//-----------------------------------------------------------------------------
template <class CT, class CI>
CT* array<CT, CI>::End(){return origin+offsets[0]-1;}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector<CT*> array<CT, CI>::End(CI dimension){
  std::vector<CT*> endings(offsets[0]/offsets[dimension]);
  for (CI i = 0; i < endings.size(); ++i){
    endings[i] = origin+(i+1)*offsets[dimension] - 1;
  }
  return endings;
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
  std::cout << "got it wrong, using pointers for both" << std::endl;
  for (CI i = 0; i <= dataend-databegin; ++i){
    *(origin+i) = *(databegin + i);
  }
}
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(std::vector<CI> index, CT value){
  std::cout << "um wtf?" << std::endl;
  *(origin + GetIndex(index)) = value;
}
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(CI index, CT value){
  std::cout << " assigning value: " << value << std::endl;
  *(origin + index) = value;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(CT value){
  for (CI i = 0; i < offsets[0]; ++i){
    *(origin + i) = value;
  }
}

/**
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(CI n, ...){
  va_list values;
  va_start(values, n);
  for (CI i = 0; i < n; ++i){
    *(origin+i) = va_arg(values, CT);
  }
  va_end(values);
}
*/
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::Transpose(CI step){
  array<CT,CI>* newarray = new array<CT,CI>(Permute(shape, step));
  for (CI i = 0; i < offsets[0]; ++i){
    newarray->Assign(Permute(i,step), *GetValue(i));
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

template <class CT, class CI>
array<CT,CI>* array<CT, CI>::InnerProduct(array<CT,CI>* operand, CI ndimensions){
  std::vector<CI> newshape(shape.size() + operand->Shape().size() - 2*ndimensions), operandshape;
  operandshape = operand->Shape();

  for (CI i = 0; i < ndimensions; ++i){
    if (shape[shape.size()-1-i] != operandshape[i]){
      std::cout << "The required dimensions did not match!" << std::endl;
    }
  }

  if (shape.size()-ndimensions > 0){
    for (CI i = 0; i < shape.size()-ndimensions; ++i){newshape[i] = shape[i];}
    for (CI i = 0; i < operandshape.size() - ndimensions; ++i){newshape[newshape.size()-1-i] = operandshape[operandshape.size()-1-i];}
  }else{
    newshape.resize(1);
    newshape[0] = 1;
  }
  array<CT,CI>* newarray = new array<CT,CI>(newshape);

  std::vector<CI> NEWINDEX(newshape.size()), THISINDEX(shape.size()), OPERANDINDEX(operandshape.size());
  
  for (CI index = 0; index < newarray->Size(); ++index){
    NEWINDEX = newarray->GetIndex(index);
    for (CI i = 0; i < shape.size()-ndimensions; ++i){THISINDEX[i] = NEWINDEX[i];}
    for (CI i = 0; i < operandshape.size() - ndimensions; ++i){OPERANDINDEX[ndimensions+i] = NEWINDEX[NEWINDEX.size()-operandshape.size() + ndimensions+i];}
    newarray->Assign(index, RecursiveInnerProduct(operand, THISINDEX, OPERANDINDEX, ndimensions));
  }

  return newarray;
}

template <class CT, class CI>
CT array<CT, CI>::RecursiveInnerProduct(array<CT,CI>* operand, std::vector<CI> THISINDEX, std::vector<CI> OPERANDINDEX, CI ndimensions){
  CT newvalue = 0;
  if (ndimensions == 1){
    CI lastnumber = THISINDEX.size()-1;
    THISINDEX[lastnumber] = 0;
    OPERANDINDEX[0] = 0;
    CT* thisbegin = GetValue(THISINDEX), *operandbegin = GetValue(OPERANDINDEX);
    CI operandstep = operand->Size(1);
    for (CI index = 0; index < shape[lastnumber]; ++index){
      newvalue += (*(thisbegin + index))*(*(operandbegin + index*operandstep));
    }
  }else{
    CI intermediatedimension = shape[shape.size()-ndimensions];
    for (CI index = 0; index < intermediatedimension; ++index){
      THISINDEX[THISINDEX.size()-ndimensions] = index;
      OPERANDINDEX[ndimensions-1] = index;
      newvalue += RecursiveInnerProduct(operand, THISINDEX, OPERANDINDEX, ndimensions-1);
    }
  }
  return newvalue;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator*(array<CT,CI>* operand){
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand->Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))*(*(operandbegin+i));
  }
  return result;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator*(CT operand){
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))*operand;
  }
  return result;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator+(array<CT,CI>* operand){
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand->Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))+(*(operandbegin+i));
  }
  return result;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator+(CT operand){
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))+operand;
  }
  return result;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator-(array<CT,CI>* operand){
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand->Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i)) - (*(operandbegin+i));
  }
  return result;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator-(CT operand){
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i)) - operand;
  }
  return result;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator/(array<CT,CI>* operand){
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand->Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))/(*(operandbegin+i));
  }
  return result;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator/(CT operand){
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))/operand;
  }
  return result;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator*=(array<CT,CI>* operand){
  CT* selfbegin = Begin(), *operandbegin = operand->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))*(*(operandbegin+i));
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator*=(CT operand){
  CT* selfbegin = Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))*operand;
  }
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator+=(array<CT,CI>* operand){
  CT* selfbegin = Begin(), *operandbegin = operand->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))+(*(operandbegin+i));
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator+=(CT operand){
  CT* selfbegin = Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))+operand;
  }
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator-=(array<CT,CI>* operand){
  CT* selfbegin = Begin(), *operandbegin = operand->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))-(*(operandbegin+i));
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator-=(CT operand){
  CT* selfbegin = Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))-operand;
  }
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator/=(array<CT,CI>* operand){
  CT* selfbegin = Begin(), *operandbegin = operand->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))/(*(operandbegin+i));
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator/=(CT operand){
  CT* selfbegin = Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(selfbegin + i) = (*(selfbegin + i))/operand;
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator[](CI index){
  std::vector<CI> newshape;
  if (shape.size() > 1){
    newshape.resize(shape.size()-1);
    newshape.assign(shape.begin()+1,shape.end());
  }else{
    newshape.push_back(1);
  }  
  array<CT,CI>* newarray = new array(newshape);
  newarray->Assign(origin+index*offsets[1], origin+(index+1)*offsets[1]-1);
  return newarray;
}

