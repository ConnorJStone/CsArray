#ifndef __CS_ARRAY__
#define __CS_ARRAY__

#include <stdarg.h>
#include <vector>
#include <iostream>
#include <stdexcept>
template<class T1, class T2> T1 min(const T1& a, const T2& b){return a>b ? a:b;}//fixme deleteme

template <class DATATYPE, class INDEXTYPE>
class array {
  //Contains the lengths of each dimension
  std::vector<INDEXTYPE> shape;
  //Pointer to the beginning of the array data
  DATATYPE* origin;
  //The lengths in memory of each dimension. offsets[0] is the length of the whole array, offsets[1] is the length along the first dimension,...,offsets[n+1] is always 1 (to represent the length of one data object) 
  std::vector<INDEXTYPE> offsets;

  //Will not print an array if it is longer than this value
  static const INDEXTYPE printlimit = 150;
  
  //Given the requested dimensions, it fills the shape and offsets objects, also creates the array in memory and assigns it to the origin
  void InitializeSelf(std::vector<INDEXTYPE> D);

  //Instead of creating new memory, it is possible to initialize the object as a subset of another array
  void InitializeAsReference(std::vector<INDEXTYPE> D, DATATYPE* begin);
  
  //Given a vector of indices, this will return a new vector with the indices permuted by an amount equal to step
  std::vector<INDEXTYPE> Permute(const std::vector<INDEXTYPE>& index, const INDEXTYPE& step) const;

  //Given a index in the memory, this will return the permuted index in memory
  INDEXTYPE Permute(const INDEXTYPE& index, const INDEXTYPE& step) const;

  //Creates a vector where accessing an element returns the permuted index
  std::vector<INDEXTYPE> GeneratePermuteMap(const INDEXTYPE& nindices, const INDEXTYPE& step) const;

  //used by operators to guarintee consistancy of implimentation
  inline array<DATATYPE,INDEXTYPE>& operate(const array<DATATYPE,INDEXTYPE>& operand, DATATYPE (array<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&)const) const;
  inline array<DATATYPE,INDEXTYPE>& operate(const DATATYPE& operand, DATATYPE (array<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&)const) const;
  inline void operateeq(const array<DATATYPE,INDEXTYPE>& operand, DATATYPE (array<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&)const) const;
  inline void operateeq(const DATATYPE& operand, DATATYPE (array<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&)const) const;
  inline DATATYPE mul(const DATATYPE& x1, const DATATYPE& x2) const{return x1*x2;}
  inline DATATYPE div(const DATATYPE& x1, const DATATYPE& x2) const{return x1/x2;}
  inline DATATYPE add(const DATATYPE& x1, const DATATYPE& x2) const{return x1+x2;}
  inline DATATYPE sub(const DATATYPE& x1, const DATATYPE& x2) const{return x1-x2;}  
  
public:

  //Empty initializer
  array();

  //Initializes the array with no guarintee about the values in its indices
  array(const std::vector<INDEXTYPE>& D);

  //Initializes the array and assigns a value to every position
  array(std::vector<INDEXTYPE> D, DATATYPE value);

  //Constructor for manual assignment of shape
  array(const INDEXTYPE& n, ...);
  
  //Destructor, will call delete on origin
  ~array();

  //Given an index along the memory, it returns a vector which has the distances along each dimension
  std::vector<INDEXTYPE> GetIndex(const INDEXTYPE& index) const;

  //Given the distances along each dimension, this will return the true index of the desired data point
  INDEXTYPE GetIndex(const std::vector<INDEXTYPE>& index) const;

  //Returns the value at an index location
  DATATYPE* GetValueP(const std::vector<INDEXTYPE>& index) const;

  //Returns the value at the index
  DATATYPE* GetValueP(const INDEXTYPE& index) const;

  //Returns the value at an index location
  DATATYPE GetValue(const std::vector<INDEXTYPE>& index) const;

  //Returns the value at the index
  DATATYPE GetValue(const INDEXTYPE& index) const;

  //Returns the size of the array along the given dimension. note1: size(0) gives the size of the whole array in memory, note2: size(n) is always 1
  INDEXTYPE Size(const INDEXTYPE& dimension = 0) const;

  //Returns the vector of values describing how long each dimension is
  std::vector<INDEXTYPE> Shape() const;

  //Returns the length of one of the dimensions, as specified
  INDEXTYPE Shape(const INDEXTYPE& dimension) const;

  //Returns the beginning and end pointer for each segment of the array according to the dimension requested 
  std::vector< std::vector<DATATYPE*> > Axis(const INDEXTYPE& dimension) const;

  //Applies a function to the elements, segmented by the axis. The function is given the begin/end pointers for the data in that segment and allowed to make changes
  void Apply(void (*f)(DATATYPE*,DATATYPE*), INDEXTYPE axis = -1);

  //Simmilar to "Apply" except the output of the function must be a DATATYPE type, the results from each function run are written to a new array. the user is not allowed to change the original data
  array<DATATYPE,INDEXTYPE>& RunFunc(DATATYPE (*f)(DATATYPE*,DATATYPE*), INDEXTYPE axis = -1) const;
  
  //Returns a pointer to the beginning of the array
  DATATYPE* Begin() const;

  //Returns a vector of pointers to the beginning of each segment on a specified dimension
  std::vector<DATATYPE*> Begin(const INDEXTYPE& dimension) const;

  //Returns a pointer to the end of the array
  DATATYPE* End() const;

  //Returns a vector of pointers to the end of each segment on a specified dimension
  std::vector<DATATYPE*> End(const INDEXTYPE& dimension) const;

  //Assign all the values from a vector
  void Assign(const std::vector<DATATYPE>& values);

  //Assign the values at the pointed to locations to the array
  void Assign(const DATATYPE* databegin, const DATATYPE* dataend);

  //Assign the given value to an index
  void Assign(const std::vector<INDEXTYPE>& index, const DATATYPE& value);

  //Assign the given value to an index
  void Assign(const INDEXTYPE& index, const DATATYPE& value);

  //Assign the given value to every element in the array
  void Assign(const DATATYPE& value);

  //Assign a list of values to the array, the first argument is the number of values to assign
  void ManualAssign(const INDEXTYPE& n, ...);

  void Resize(const std::vector<INDEXTYPE>& D);

  void Resize(const INDEXTYPE& n, ...);

  //Returns a new array which is the transpose of the old one by a given number of steps (if you are using matrices, the default is what you want)
  array<DATATYPE,INDEXTYPE>* Transpose(INDEXTYPE step = 1) const;

  void Print() const;

  //array<DATATYPE,INDEXTYPE>* TensorProduct(const array<DATATYPE,INDEXTYPE>& operand);//future
  
  array<DATATYPE,INDEXTYPE>& operator*(const array<DATATYPE,INDEXTYPE>& operand) const;

  array<DATATYPE,INDEXTYPE>& operator*(const DATATYPE& operand) const;
  
  void operator*=(const array<DATATYPE,INDEXTYPE>& operand);
  
  void operator*=(const DATATYPE& operand);

  array<DATATYPE,INDEXTYPE>& operator+(const array<DATATYPE,INDEXTYPE>& operand) const;

  array<DATATYPE,INDEXTYPE>& operator+(const DATATYPE& operand) const;

  void operator+=(const array<DATATYPE,INDEXTYPE>& operand);

  void operator+=(const DATATYPE& operand);

  array<DATATYPE,INDEXTYPE>& operator-(const array<DATATYPE,INDEXTYPE>& operand) const;

  array<DATATYPE,INDEXTYPE>& operator-(const DATATYPE& operand) const;

  void operator-=(const array<DATATYPE,INDEXTYPE>& operand);

  void operator-=(const DATATYPE& operand);

  array<DATATYPE,INDEXTYPE>& operator/(const array<DATATYPE,INDEXTYPE>& operand) const;

  array<DATATYPE,INDEXTYPE>& operator/(const DATATYPE& operand) const;

  void operator/=(const array<DATATYPE,INDEXTYPE>& operand);

  void operator/=(const DATATYPE& operand);

  array<DATATYPE,INDEXTYPE>& operator()(const INDEXTYPE& index) const;

  DATATYPE operator[](const INDEXTYPE& index) const;
  
  array<DATATYPE,INDEXTYPE>& operator=(const array<DATATYPE,INDEXTYPE>& operand);
};


//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::InitializeSelf(std::vector<INDEXTYPE> D){
  for (INDEXTYPE i = 0; i < D.size(); ++i){
    if (D[i] == (DATATYPE) 0){throw std::invalid_argument( "Array cannot have dimension length of 0" );}
  }
  
  if (D.size() == 0){
    this->shape.resize(1);
    this->shape[0] = (INDEXTYPE) 1;
    this->offsets.resize(1);
    this->offsets[0] = (INDEXTYPE) 1;
  }else{
    this->shape.resize(D.size());
    this->shape.assign(D.begin(), D.end());
    this->offsets.resize(this->shape.size()+1);
    for (INDEXTYPE a = 0; a < this->offsets.size(); ++a){
      INDEXTYPE axislength = 1;
      for (INDEXTYPE d = a; d < this->shape.size(); ++d){axislength *= this->shape[d];}
      this->offsets[a] = axislength;
    }
  }
  this->origin = new DATATYPE[this->offsets[0]];
}

template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::InitializeAsReference(std::vector<INDEXTYPE> D, DATATYPE* begin){
  if (D.size() == 0){
    this->shape.resize(1);
    this->shape[0] = (INDEXTYPE) 1;
    this->offsets.resize(1);
    this->offsets[0] = (INDEXTYPE) 1;
  }else{
    this->shape.resize(D.size());
    this->shape.assign(D.begin(), D.end());
    this->offsets.resize(this->shape.size()+1);
    for (INDEXTYPE a = 0; a < this->offsets.size(); ++a){
      INDEXTYPE axislength = 1;
      for (INDEXTYPE d = a; d < this->shape.size(); ++d){axislength *= this->shape[d];}
      this->offsets[a] = axislength;
    }
  }
  this->origin = begin;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
std::vector<INDEXTYPE> array<DATATYPE, INDEXTYPE>::GetIndex(const INDEXTYPE& index) const{
  std::vector<INDEXTYPE> newindex(this->shape.size());
  for (INDEXTYPE i = 0; i < newindex.size(); ++i){
    newindex[i] = (index%this->offsets[i])/this->offsets[i+1];
  }
  return newindex;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
INDEXTYPE array<DATATYPE, INDEXTYPE>::GetIndex(const std::vector<INDEXTYPE>& index) const{
  INDEXTYPE newindex = 0;
  for (INDEXTYPE i = 0; i < index.size(); ++i){
    newindex += index[i]*this->offsets[i+1];
  }
  return newindex;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
std::vector<INDEXTYPE> array<DATATYPE, INDEXTYPE>::Permute(const std::vector<INDEXTYPE>& index, const INDEXTYPE& step) const{
  std::vector<INDEXTYPE> newindex(index.size());
  for (INDEXTYPE i = 0; i < index.size(); ++i){
    newindex[(i+step)%index.size()] = index[i];
  }
  return newindex;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
INDEXTYPE array<DATATYPE, INDEXTYPE>::Permute(const INDEXTYPE& index, const INDEXTYPE& step) const{
  return this->GetIndex(this->Permute(this->GetIndex(index), step));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
std::vector<INDEXTYPE> array<DATATYPE, INDEXTYPE>::GeneratePermuteMap(const INDEXTYPE& nindices, const INDEXTYPE& step) const{
  std::vector<INDEXTYPE> map(nindices);
  for (INDEXTYPE i = 0; i < nindices; ++i){
    map[i] = (i+step)%nindices;
  }
  return map;
}


//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE, INDEXTYPE>::array(){
  std::vector<INDEXTYPE> D(1,1);
  this->InitializeSelf(D);
}


//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE, INDEXTYPE>::array(const std::vector<INDEXTYPE>& D){
  this->InitializeSelf(D);
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE, INDEXTYPE>::array(std::vector<INDEXTYPE> D, DATATYPE value){
  this->InitializeSelf(D);
  this->Assign(value);
}

template <class DATATYPE, class INDEXTYPE>
array<DATATYPE, INDEXTYPE>::array(const INDEXTYPE& n, ...){
  std::vector<INDEXTYPE> D(n);
  va_list values;
  va_start(values, n);
  for (INDEXTYPE i = 0; i < n; ++i){
    D[i] = va_arg(values, INDEXTYPE);
  }
  InitializeSelf(D);
  va_end(values);
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE, INDEXTYPE>::~array(){delete[] this->origin;}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE* array<DATATYPE, INDEXTYPE>::GetValueP(const std::vector<INDEXTYPE>& index) const{
  return (this->origin + this->GetIndex(index));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE* array<DATATYPE, INDEXTYPE>::GetValueP(const INDEXTYPE& index) const{
  return (this->origin + index);
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE array<DATATYPE, INDEXTYPE>::GetValue(const std::vector<INDEXTYPE>& index) const{
  return *(this->origin + this->GetIndex(index));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE array<DATATYPE, INDEXTYPE>::GetValue(const INDEXTYPE& index) const{
  return *(this->origin + index);
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
INDEXTYPE array<DATATYPE, INDEXTYPE>::Size(const INDEXTYPE& dimension) const{
  return this->offsets[dimension];
}
  
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
std::vector<INDEXTYPE> array<DATATYPE, INDEXTYPE>::Shape() const{
  return this->shape;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
INDEXTYPE array<DATATYPE, INDEXTYPE>::Shape(const INDEXTYPE& dimension) const{
  return this->shape[dimension];
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
std::vector< std::vector<DATATYPE*> > array<DATATYPE, INDEXTYPE>::Axis(const INDEXTYPE& dimension) const{
  std::vector< std::vector<DATATYPE*> > axis(this->offsets[0]/this->offsets[dimension]);
  for (INDEXTYPE i = 0; i < axis.size(); ++i){
    axis[i].resize(2);
    axis[i][0] = this->origin + i*this->offsets[dimension];
    axis[i][1] = this->origin + (i+1)*this->offsets[dimension] - 1;
  }
  return axis;
}

///**
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE,INDEXTYPE>::Apply(void (*f)(DATATYPE*,DATATYPE*), INDEXTYPE axis){
  std::vector< std::vector<DATATYPE*> > axispointers;
  if (axis == -1){axispointers = this->Axis(this->shape.size()-1);}
  else{axispointers = this->Axis(axis);}
  
  for (INDEXTYPE i = 0; i < axispointers.size(); ++i){
    (*f)(axispointers[i][0], axispointers[i][1]);
  }
}

template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>& array<DATATYPE,INDEXTYPE>::RunFunc(DATATYPE (*f)(DATATYPE*,DATATYPE*), INDEXTYPE axis) const{
  std::vector< std::vector<DATATYPE*> > axispointers;
  std::vector<INDEXTYPE> newshape;
  if (axis == -1){axispointers = this->Axis(this->shape.size()-1);}
  else{axispointers = this->Axis(axis);}

  //fixme, it should work for all axes values
  if (this->shape.size() == 1){
    newshape.resize(1);
    newshape[0] = 1;
  }else{
    newshape.resize(this->shape.size()-1);
    newshape.assign(this->shape.begin(), this->shape.end()-1);
  }
  array<DATATYPE,INDEXTYPE>* newarray = new array<DATATYPE,INDEXTYPE>(newshape);
  
  for (INDEXTYPE i = 0; i < axispointers.size(); ++i){
    newarray->Assign(i, (*f)(axispointers[i][0], axispointers[i][1]));
  }
  return *newarray;
}
//*/
/**
template <class DATATYPE, class INDEXTYPE>
const array<DATATYPE*,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::Axis(const INDEXTYPE& dimension) const{
  std::vector<INDEXTYPE> axisshape(2);
  axisshape[0] = this->offsets[0]/this->offsets[dimension];
  axisshape[1] = (INDEXTYPE) 2;
  array<DATATYPE*,INDEXTYPE>* axis = new array<DATATYPE,INDEXTYPE>(axisshape);
  for (INDEXTYPE i = 0; i < axis.size(); ++i){
    axis->Assign(2*i, this->origin + i*this->offsets[dimension]);
    axis->Assign(2*i + 1, this->origin + (i+1)*this->offsets[dimension] - 1);
  }
  return *axis;
}
*/

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE* array<DATATYPE, INDEXTYPE>::Begin() const{return this->origin;}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
std::vector<DATATYPE*> array<DATATYPE, INDEXTYPE>::Begin(const INDEXTYPE& dimension) const{
  std::vector<DATATYPE*> beginnings(this->offsets[0]/this->offsets[dimension]);
  for (INDEXTYPE i = 0; i < beginnings.size(); ++i){
    beginnings[i] = this->origin+i*this->offsets[dimension];
  }
  return beginnings;
}
    

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE* array<DATATYPE, INDEXTYPE>::End() const{return this->origin+this->offsets[0]-1;}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
std::vector<DATATYPE*> array<DATATYPE, INDEXTYPE>::End(const INDEXTYPE& dimension) const{
  std::vector<DATATYPE*> endings(this->offsets[0]/this->offsets[dimension]);
  for (INDEXTYPE i = 0; i < endings.size(); ++i){
    endings[i] = this->origin+(i+1)*this->offsets[dimension] - 1;
  }
  return endings;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::Assign(const std::vector<DATATYPE>& values){
  for (INDEXTYPE i = 0; i < values.size(); ++i){
    *(this->origin+i) = values[i];
  }
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::Assign(const DATATYPE* databegin, const DATATYPE* dataend){
  for (INDEXTYPE i = 0; i < min(dataend-databegin+1, this->offsets[0]); ++i){
    *(this->origin+i) = *(databegin + i);
  }
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::Assign(const std::vector<INDEXTYPE>& index, const DATATYPE& value){
  *(this->origin + this->GetIndex(index)) = value;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::Assign(const INDEXTYPE& index, const DATATYPE& value){
  *(this->origin + index) = value;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::Assign(const DATATYPE& value){
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    *(this->origin + i) = value;
  }
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::ManualAssign(const INDEXTYPE& n, ...){
  va_list values;
  va_start(values, n);
  for (INDEXTYPE i = 0; i < n; ++i){
    *(this->origin+i) = va_arg(values, DATATYPE);
  }
  va_end(values);
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::Resize(const std::vector<INDEXTYPE>& D){
  bool allequal = true;
  if (D.size() != this->shape.size()){allequal = false;}
  else{for (INDEXTYPE i = 0; i < D.size(); ++i){if (D[i] != this->shape[i]){allequal = false; break;}}}
  
  if (!allequal){
    std::vector<INDEXTYPE> oldshape = this->shape, oldoffsets = this->offsets;
    DATATYPE* oldorigin = this->origin;
    this->InitializeSelf(D);
    this->Assign(oldorigin, oldorigin + min(oldoffsets[0], this->offsets[0]) - 1);
    delete[] oldorigin;
  }
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::Resize(const INDEXTYPE& n, ...){
  std::vector<INDEXTYPE> D(n);
  va_list values;
  va_start(values, n);
  for (INDEXTYPE i = 0; i < n; ++i){
    D[i] = va_arg(values, INDEXTYPE);
  }
  InitializeSelf(D);
  va_end(values);
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>* array<DATATYPE, INDEXTYPE>::Transpose(INDEXTYPE step) const{
  array<DATATYPE,INDEXTYPE>* newarray = new array<DATATYPE,INDEXTYPE>(this->Permute(this->shape, step));
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    newarray->Assign(this->Permute(i,step), this->GetValue(i));
  }
  return newarray;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::Print() const{
  if (this->offsets[0] > printlimit && printlimit >= 0){std::cout<<"array too long: " << this->offsets[0] << " elements" << std::endl; return;}
  
  bool comma = true;
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    comma = true;
    for (INDEXTYPE dim = 0; dim < this->offsets.size()-1; ++dim){
      if (i%this->offsets[dim] == 0){std::cout << "|";comma = false;}
    }
    if (i != 0 && comma){std::cout << ", ";}
    std::cout << this->GetValue(i);
  }
  for (INDEXTYPE dim = 0; dim < this->offsets.size()-1; ++dim){std::cout << "|";}
  std::cout<<std::endl;
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operator*(const array<DATATYPE,INDEXTYPE>& operand) const{
  return operate(operand, (&array<DATATYPE, INDEXTYPE>::mul));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operator*(const DATATYPE& operand) const{
  return operate(operand, (&array<DATATYPE, INDEXTYPE>::mul));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operator+(const array<DATATYPE,INDEXTYPE>& operand) const{
  return operate(operand, (&array<DATATYPE, INDEXTYPE>::add));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operator+(const DATATYPE& operand) const{
  return operate(operand, (&array<DATATYPE, INDEXTYPE>::add));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operator-(const array<DATATYPE,INDEXTYPE>& operand) const{
  return operate(operand, (&array<DATATYPE, INDEXTYPE>::sub));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operator-(const DATATYPE& operand) const{
  return operate(operand, (&array<DATATYPE, INDEXTYPE>::sub));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operator/(const array<DATATYPE,INDEXTYPE>& operand) const{
  return operate(operand, (&array<DATATYPE, INDEXTYPE>::div));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operator/(const DATATYPE& operand) const{
  return operate(operand, (&array<DATATYPE, INDEXTYPE>::div));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::operator*=(const array<DATATYPE,INDEXTYPE>& operand){
  return operateeq(operand, (&array<DATATYPE, INDEXTYPE>::mul));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::operator*=(const DATATYPE& operand){
  return operateeq(operand, (&array<DATATYPE, INDEXTYPE>::mul));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::operator+=(const array<DATATYPE,INDEXTYPE>& operand){
  return operateeq(operand, (&array<DATATYPE, INDEXTYPE>::add));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::operator+=(const DATATYPE& operand){
  return operateeq(operand, (&array<DATATYPE, INDEXTYPE>::add));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::operator-=(const array<DATATYPE,INDEXTYPE>& operand){
  return operateeq(operand, (&array<DATATYPE, INDEXTYPE>::sub));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::operator-=(const DATATYPE& operand){
  return operateeq(operand, (&array<DATATYPE, INDEXTYPE>::sub));
}
    
//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::operator/=(const array<DATATYPE,INDEXTYPE>& operand){
  return operateeq(operand, (&array<DATATYPE, INDEXTYPE>::div));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
void array<DATATYPE, INDEXTYPE>::operator/=(const DATATYPE& operand){
  return operateeq(operand, (&array<DATATYPE, INDEXTYPE>::div));
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operator()(const INDEXTYPE& index) const{
  std::vector<INDEXTYPE> newshape;
  if (this->shape.size() > 1){
    newshape.assign(this->shape.begin()+1,this->shape.end());
  }else{
    throw std::invalid_argument( "Cannot construct sub-array from 1D array, (Perhaps you want to use element access [] operator?)" );
  }
  
  array<DATATYPE,INDEXTYPE>* newarray = new array();
  newarray->InitializeAsReference(newshape, origin + index*offsets[1]);
  return *newarray;
}


//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
DATATYPE array<DATATYPE, INDEXTYPE>::operator[](const INDEXTYPE& index) const{
  return *(this->origin + index);
}

//-----------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operator=(const array<DATATYPE,INDEXTYPE>& operand){
  if (this != &operand){
    delete[] this->origin;
    this->InitializeSelf(operand.Shape());
    this->Assign(operand.Begin(), operand.End());
  }
  return *this;
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------
template <class DATATYPE, class INDEXTYPE>
inline array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operate(const array<DATATYPE,INDEXTYPE>& operand, DATATYPE (array<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&) const) const{
  array<DATATYPE,INDEXTYPE>* result = new array(this->shape);
      
  DATATYPE *operandbegin = operand.Begin(), *resultbegin = result->Begin();
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    *(resultbegin + i) = (this->*op)(*(this->origin + i),*(operandbegin+i));
  }
  return *result;
}

template <class DATATYPE, class INDEXTYPE>
inline array<DATATYPE,INDEXTYPE>& array<DATATYPE, INDEXTYPE>::operate(const DATATYPE& operand, DATATYPE (array<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&) const) const{
  array<DATATYPE,INDEXTYPE>* result = new array(this->shape);
      
  DATATYPE *resultbegin = result->Begin();
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    *(resultbegin + i) = (this->*op)(*(this->origin + i), operand);
  }
  return *result;
}

template <class DATATYPE, class INDEXTYPE>
inline void array<DATATYPE, INDEXTYPE>::operateeq(const array<DATATYPE,INDEXTYPE>& operand, DATATYPE (array<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&) const) const{
  DATATYPE *operandbegin = operand.Begin();
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    *(this->origin + i) = (this->*op)(*(this->origin + i),*(operandbegin+i));
  }
}

template <class DATATYPE, class INDEXTYPE>
inline void array<DATATYPE, INDEXTYPE>::operateeq(const DATATYPE& operand, DATATYPE (array<DATATYPE, INDEXTYPE>::*op)(const DATATYPE&, const DATATYPE&) const) const{
  for (INDEXTYPE i = 0; i < this->offsets[0]; ++i){
    *(this->origin + i) = (this->*op)(*(this->origin + i), operand);
  }
}

#endif
