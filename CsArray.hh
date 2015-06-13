#ifndef __CS_ARRAY__
#define __CS_ARRAY__

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

  //Will not print an array if it is longer than this value
  static const CI printlimit = 150;
  
  //Given the requested dimensions, it fills the shape and offsets objects, also creates the array in memory and assigns it to the origin
  void InitializeSelf(std::vector<CI> D);

  //Given a vector of indices, this will return a new vector with the indices permuted by an amount equal to step
  std::vector<CI> Permute(const std::vector<CI>& index, const CI& step) const;

  //Given a index in the memory, this will return the permuted index in memory
  CI Permute(const CI& index, const CI& step) const;

  //Creates a vector where accessing an element returns the permuted index
  std::vector<CI> GeneratePermuteMap(const CI& nindices, const CI& step) const;
    
public:

  //Empty initializer
  array();

  //Initializes the array with no guarintee about the values in its indices
  array(const std::vector<CI>& D);

  //Initializes the array and assigns a value to every position
  array(std::vector<CI> D, CT value);

  //Destructor, will call delete on origin
  ~array();

  //Given an index along the memory, it returns a vector which has the distances along each dimension
  std::vector<CI> GetIndex(const CI& index) const;

  //Given the distances along each dimension, this will return the true index of the desired data point
  CI GetIndex(const std::vector<CI>& index) const;

  //Returns the value at an index location
  CT* GetValueP(const std::vector<CI>& index) const;

  //Returns the value at the index
  CT* GetValueP(const CI& index) const;

  //Returns the value at an index location
  CT GetValue(const std::vector<CI>& index) const;

  //Returns the value at the index
  CT GetValue(const CI& index) const;

  //Returns the size of the array along the given dimension. note1: size(0) gives the size of the whole array in memory, note2: size(n) is always 1
  const CI Size(const CI& dimension = 0) const;

  //Returns the vector of values describing how long each dimension is
  const std::vector<CI> Shape() const;

  //Returns the length of one of the dimensions, as specified
  const CI Shape(const CI& dimension) const;

  //Returns the beginning and end pointer for each segment of the array according to the dimension requested 
  const std::vector< std::vector<CT*> > Axis(const CI& dimension) const;

  //Returns a pointer to the beginning of the array
  const CT* Begin() const;

  //Returns a vector of pointers to the beginning of each segment on a specified dimension
  const std::vector<CT*> Begin(const CI& dimension) const;

  //Returns a pointer to the end of the array
  const CT* End() const;

  //Returns a vector of pointers to the end of each segment on a specified dimension
  const std::vector<CT*> End(const CI& dimension) const;

  //Assign all the values from a vector
  void Assign(const std::vector<CT>& values);

  //Assign the values at the pointed to locations to the array
  void Assign(CT* databegin, CT* dataend);

  //Assign the given value to an index
  void Assign(const std::vector<CI>& index, const CT& value);

  //Assign the given value to an index
  void Assign(const CI& index, const CT& value);

  //Assign the given value to every element in the array
  void Assign(const CT& value);

  //Assign a list of values to the array, the first argument is the number of values to assign
  void ManualAssign(const CI& n, ...);

  //Returns a new array which is the transpose of the old one by a given number of steps (if you are using matrices, the default is what you want)
  array<CT,CI>* Transpose(CI step = 1) const;

  void Print() const;

  //array<CT,CI>* TensorProduct(const array<CT,CI>& operand);//future
  
  array<CT,CI>* operator*(const array<CT,CI>* operand) const;

  array<CT,CI>* operator*(CT operand) const;
  
  void operator*=(const array<CT,CI>* operand);
  
  void operator*=(CT operand);

  array<CT,CI>* operator+(const array<CT,CI>* operand) const;

  array<CT,CI>* operator+(CT operand) const;

  void operator+=(const array<CT,CI>* operand);

  void operator+=(CT operand);

  array<CT,CI>* operator-(const array<CT,CI>* operand) const;

  array<CT,CI>* operator-(CT operand) const;

  void operator-=(const array<CT,CI>* operand);

  void operator-=(CT operand);

  array<CT,CI>* operator/(const array<CT,CI>* operand) const;

  array<CT,CI>* operator/(CT operand) const;

  void operator/=(const array<CT,CI>* operand);

  void operator/=(CT operand);

  array<CT,CI>* operator[](CI index) const;

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
std::vector<CI> array<CT, CI>::GetIndex(const CI& index) const{
  std::vector<CI> newindex(shape.size());
  for (CI i = 0; i < newindex.size(); ++i){
    newindex[i] = (index%offsets[i])/offsets[i+1];
  }
  return newindex;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CI array<CT, CI>::GetIndex(const std::vector<CI>& index) const{
  CI newindex = 0;
  for (CI i = 0; i < index.size(); ++i){
    newindex += index[i]*offsets[i+1];
  }
  return newindex;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector<CI> array<CT, CI>::Permute(const std::vector<CI>& index, const CI& step) const{
  std::vector<CI> newindex(index.size());
  for (CI i = 0; i < index.size(); ++i){
    newindex[(i+step)%index.size()] = index[i];
  }
  return newindex;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CI array<CT, CI>::Permute(const CI& index, const CI& step) const{
  return GetIndex(Permute(GetIndex(index), step));
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
std::vector<CI> array<CT, CI>::GeneratePermuteMap(const CI& nindices, const CI& step) const{
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
array<CT, CI>::array(const std::vector<CI>& D){
  InitializeSelf(D);
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT, CI>::array(std::vector<CI> D, CT value){
  InitializeSelf(D);
  Assign(value);
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT, CI>::~array(){delete[] origin;}


//-----------------------------------------------------------------------------
template <class CT, class CI>
CT* array<CT, CI>::GetValueP(const std::vector<CI>& index) const{
  return (origin + GetIndex(index));
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CT* array<CT, CI>::GetValueP(const CI& index) const{
  return (origin + index);
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CT array<CT, CI>::GetValue(const std::vector<CI>& index) const{
  return *(origin + GetIndex(index));
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
CT array<CT, CI>::GetValue(const CI& index) const{
  return *(origin + index);
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
const CI array<CT, CI>::Size(const CI& dimension) const{
  return offsets[dimension];
}
  
//-----------------------------------------------------------------------------
template <class CT, class CI>
const std::vector<CI> array<CT, CI>::Shape() const{
  return shape;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
const CI array<CT, CI>::Shape(const CI& dimension) const{
  return shape[dimension];
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
const std::vector< std::vector<CT*> > array<CT, CI>::Axis(const CI& dimension) const{
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
const CT* array<CT, CI>::Begin() const{return origin;}

//-----------------------------------------------------------------------------
template <class CT, class CI>
const std::vector<CT*> array<CT, CI>::Begin(const CI& dimension) const{
  std::vector<CT*> beginnings(offsets[0]/offsets[dimension]);
  for (CI i = 0; i < beginnings.size(); ++i){
    beginnings[i] = origin+i*offsets[dimension];
  }
  return beginnings;
}
    

//-----------------------------------------------------------------------------
template <class CT, class CI>
const CT* array<CT, CI>::End() const{return origin+offsets[0]-1;}

//-----------------------------------------------------------------------------
template <class CT, class CI>
const std::vector<CT*> array<CT, CI>::End(const CI& dimension) const{
  std::vector<CT*> endings(offsets[0]/offsets[dimension]);
  for (CI i = 0; i < endings.size(); ++i){
    endings[i] = origin+(i+1)*offsets[dimension] - 1;
  }
  return endings;
}
/**
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(CI selfstartatindex, CT* databegin, CT* dataend){
  for (CI i = 0; i <= dataend-databegin; ++i){
    *(origin+selfstartatindex+i) = *(databegin + i);
  }
}
*/

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(const std::vector<CT>& values){
  for (CI i = 0; i < values.size(); ++i){
    *(origin+i) = values[i];
  }
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(CT* databegin, CT* dataend){
  for (CI i = 0; i < dataend-databegin+1; ++i){
    *(origin+i) = *(databegin + i);
  }
}
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(const std::vector<CI>& index, const CT& value){
  *(origin + GetIndex(index)) = value;
}
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(const CI& index, const CT& value){
  *(origin + index) = value;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::Assign(const CT& value){
  for (CI i = 0; i < offsets[0]; ++i){
    *(origin + i) = value;
  }
}


//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::ManualAssign(const CI& n, ...){
  va_list values;
  va_start(values, n);
  for (CI i = 0; i < n; ++i){
    *(origin+i) = va_arg(values, CT);
  }
  va_end(values);
}

    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::Transpose(CI step) const{
  array<CT,CI>* newarray = new array<CT,CI>(Permute(shape, step));
  for (CI i = 0; i < offsets[0]; ++i){
    newarray->Assign(Permute(i,step), GetValue(i));
  }
  return newarray;
}

//-----------------------------------------------------------------------------

template <class CT, class CI>
void array<CT, CI>::Print() const{
  if (offsets[0] > printlimit){std::cout<<"array too long: " << offsets[0] << " elements" << std::endl; return;}
  
  bool comma = true;
  for (CI i = 0; i < offsets[0]; ++i){
    comma = true;
    for (CI dim = 0; dim < offsets.size()-1; ++dim){
      if (i%offsets[dim] == 0){std::cout << "|";comma = false;}
    }
    if (i != 0 && comma){std::cout << ", ";}
    std::cout << GetValue(i);
  }
  for (CI dim = 0; dim < offsets.size()-1; ++dim){std::cout << "|";}
  std::cout<<std::endl;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator*(const array<CT,CI>* operand) const{
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand->Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))*(*(operandbegin+i));
  }
  return result;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator*(CT operand) const{
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))*operand;
  }
  return result;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator+(const array<CT,CI>* operand) const{
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand->Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))+(*(operandbegin+i));
  }
  return result;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator+(CT operand) const{
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))+operand;
  }
  return result;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator-(const array<CT,CI>* operand) const{
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand->Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i)) - (*(operandbegin+i));
  }
  return result;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator-(CT operand) const{
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i)) - operand;
  }
  return result;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator/(const array<CT,CI>* operand) const{
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *operandbegin = operand->Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))/(*(operandbegin+i));
  }
  return result;
}

//-----------------------------------------------------------------------------
template <class CT, class CI>
array<CT,CI>* array<CT, CI>::operator/(CT operand) const{
  array<CT,CI>* result = new array(shape);
      
  CT* selfbegin = Begin(), *resultbegin = result->Begin();
  for (CI i = 0; i < offsets[0]; ++i){
    *(resultbegin + i) = (*(selfbegin + i))/operand;
  }
  return result;
}
    
//-----------------------------------------------------------------------------
template <class CT, class CI>
void array<CT, CI>::operator*=(const array<CT,CI>* operand){
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
void array<CT, CI>::operator+=(const array<CT,CI>* operand){
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
void array<CT, CI>::operator-=(const array<CT,CI>* operand){
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
void array<CT, CI>::operator/=(const array<CT,CI>* operand){
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
array<CT,CI>* array<CT, CI>::operator[](CI index) const{
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

#endif
