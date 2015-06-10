
#include <vector>

namespace CsUtil{
  template <class CT, class CI>
  class array {
    //--------------------------------------------------
    std::vector<CI> shape;
    CT* origin;
    std::vector<CI> X;
    
    //--------------------------------------------------
    void FillX(){
      X.resize(shape.size()+1);
      for (CI a = 0; a < X.size(); ++a){
	CI distance = 1;
	for (CI d = a; d < shape.size(); ++d){distance *= shape[d];}
	X[a] = distance;
      }
    }
    
    //--------------------------------------------------
    //template <class I>
    std::vector<int> GetLocation(int i){
      std::vector<int> location(shape.size());
      for (CI a = 0; a < location.size(); ++a){
	location[a] = (i%X[a])/X[a+1];
      }
      return location;
    }
    //template <class I>
    int GetLocation(std::vector<int> index){
      int location = 0;
      for (CI i = 0; i < index.size(); ++i){
	location += index[i]*X[i+1];
      }
      return location;
    }
    
    //--------------------------------------------------
    //template <class I>
    bool InRange(int location, std::vector< std::vector<int> > indices){
      return InRange(GetLocation(location), indices);
    }
    //template <class I>
    bool InRange(std::vector<int> location, std::vector< std::vector<int> > indices){
      for (CI i = 0; i < location.size(); ++i){
	if (location[i] < indices[i][0] || location[i] > indices[i][1]){return false;}
      }
      return true;
    }
    
    //--------------------------------------------------
    //template <class I, class S>
    std::vector<int> Permute(std::vector<int> index, int step){
      std::vector<int> newindex(index.size());
      for (CI i = 0; i < index.size(); ++i){
	newindex[(i+step)%index.size()] = index[i];
      }
      return newindex;
    }
    //template <class I, class S>
    int Permute(int index, int step){
      return GetLocation(Permute(GetLocation(index), step));
    }
    
  public:
    //--------------------------------------------------
    array(){}//shape.resize(1);shape[0] = 1;FillX();origin = new CT[X[0]];
    
    ////template <class I>
    array(std::vector<int> D){
      shape.resize(D.size());
      shape.assign(D.begin(), D.end());
      FillX();
      origin = new CT[X[0]];
    }
    
    ~array(){delete[] origin;}
    
    //--------------------------------------------------
    //template <class I>
    CT GetElement(std::vector<int> index){
      return *(origin + GetLocation(index));
    }
    //template <class I>
    CT GetElement(int index){
      return *(origin + index);
    }
    
    //--------------------------------------------------
    //this can be done better, fixme
    std::vector< std::vector<CT*> > GetSubArray(std::vector< std::vector<int> > indices){
      std::vector< std::vector<CT*> > subarray;
      std::vector<int> startat(indices.size()),endat(indices.size());
      for (int i = 0; i < startat.size(); ++i){
	startat[i] = indices[i][0];
	endat[i] = indices[i][1];
      }
      int start = GetLocation(startat), end = GetLocation(endat), basestep = indices[indices.size()-1][1] - indices[indices.size()-1][0];
      
      for (int i = start; i <= end; ++i){
	if (InRange(i,indices)){
	  subarray.push_back(std::vector<CT*> (2,origin+i));
	  i += basestep;
	  subarray[subarray.size()-1][1] = origin + i;
	}
      }
      return subarray;
    }
    
    //--------------------------------------------------
    //not working yet, fixme
    array<CT,CI> SubArray(std::vector< std::vector<int> > indices){
      std::vector<int> subshape(indices.size());
      for (int i = 0; i < subshape.size(); ++i){
	subshape[i] = indices[i][1] - indices[i][0] + 1;
      }
      array<CT,CI> subarray(subshape);
      
      return subarray;
    }
    
    //--------------------------------------------------
    //template <class I>
    CI Size(int dimension = 0){
      return X[dimension];
    }

    //--------------------------------------------------
    //template <class I>
    std::vector< std::vector<CT*> > Axis(int dimension){
      std::vector< std::vector<CT*> > axis(X[0]/X[dimension]);
      for (CI i = 0; i < axis.size(); ++i){
	axis[i].resize(2);
	axis[i][0] = origin + i*X[dimension];
	axis[i][1] = origin + (i+1)*X[dimension] - 1;
      }
      return axis;
    }
    
    //--------------------------------------------------
    CT* Begin(){return origin;}
    //template <class I>
    std::vector<CT*> Begin(int dimension){
      std::vector<CT*> beginnings(X[0]/X[dimension]);
      for (CI i = 0; i < beginnings.size(); ++i){
	beginnings[i] = origin+i*X[dimension];
      }
      return beginnings;
    }
    
    //--------------------------------------------------
    CT* End(){return origin+X[0]-1;}
    //template <class I>
    std::vector<CT*> End(int dimension){
      std::vector<CT*> endings(X[0]/X[dimension]);
      for (CI i = 0; i < endings.size(); ++i){
	endings[i] = origin+(i+1)*X[dimension] - 1;
      }
      return endings;
    }
    
    //--------------------------------------------------
    void Assign(CT* selfstartat, CT* databegin, CT* dataend){
      for (CI i = 0; i <= dataend-databegin; ++i){
	*(selfstartat+i) = *(databegin + i);
      }
    }
    void Assign(int selfstartatindex, CT* databegin, CT* dataend){
      for (CI i = 0; i <= dataend-databegin; ++i){
	*(origin+selfstartatindex+i) = *(databegin + i);
      }
    }
    void Assign(CT* databegin, CT* dataend){
      for (CI i = 0; i <= dataend-databegin; ++i){
	*(origin+i) = *(databegin + i);
      }
    }
    //template <class I>
    void Assign(std::vector<int> index, CT value){
      *(origin + GetLocation(index)) = value;
    }
    //template <class I>
    void Assign(int index, CT value){
      *(origin + index) = value;
    }
    
    //--------------------------------------------------
    //template <class I>
    array<CT,CI> Transpose(int step = 1){
      array<CT,CI> newarray(Permute(shape, step));
      for (CI i = 0; i < X[0]; ++i){
	newarray.Assign(Permute(i,step), GetElement(i));
      }
      return newarray;
    }

    //--------------------------------------------------
    string Print(){
      string out = "";
      return out;
    }

    //fixme, look at tensor contraction
    //template <class OCT, class OCI>
    array<CT,CI> TensorContract(array<CT,CI> operand){
      return *this;
    }
    //--------------------------------------------------
    //--------------------------------------------------
    //template <class OCT, class OCI>
    array<CT,CI> operator*(array<CT,CI> operand){
      array<CT,CI> product(shape);
      
      CT* selfbegin = Begin(), *operandbegin = operand.Begin(), *productbegin = product.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(productbegin + i) = (*(selfbegin + i))*(*(operandbegin+i));
      }
      return product;
    }
    //template <class T>
    array<CT,CI> operator*(CT operand){
      array<CT,CI> product(shape);
      
      CT* selfbegin = Begin(), *productbegin = product.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(productbegin + i) = (*(selfbegin + i))*operand;
      }
      return product;
    }
    
    //--------------------------------------------------
    //template <class OCT, class OCI>
    array<CT,CI> operator+(array<CT,CI> operand){
      array<CT,CI> addition(shape);
      
      CT* selfbegin = Begin(), *operandbegin = operand.Begin(), *additionbegin = addition.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(additionbegin + i) = (*(selfbegin + i))+(*(operandbegin+i));
      }
      return addition;
    }
    //template <class T>
    array<CT,CI> operator+(CT operand){
      array<CT,CI> addition(shape);
      
      CT* selfbegin = Begin(), *additionbegin = addition.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(additionbegin + i) = (*(selfbegin + i))+operand;
      }
      return addition;
    }
    
    //--------------------------------------------------
    //template <class OCT, class OCI>
    array<CT,CI> operator-(array<CT,CI> operand){
      array<CT,CI> difference(shape);
      
      CT* selfbegin = Begin(), *operandbegin = operand.Begin(), *differencebegin = difference.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(differencebegin + i) = (*(selfbegin + i)) - (*(operandbegin+i));
      }
      return difference;
    }
    //template <class T>
    array<CT,CI> operator-(CT operand){
      array<CT,CI> difference(shape);
      
      CT* selfbegin = Begin(), *differencebegin = difference.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(differencebegin + i) = (*(selfbegin + i)) - operand;
      }
      return difference;
    }
    
    //--------------------------------------------------
    //template <class OCT, class OCI>
    array<CT,CI> operator/(array<CT,CI> operand){
      array<CT,CI> quotient(shape);
      
      CT* selfbegin = Begin(), *operandbegin = operand.Begin(), *quotientbegin = quotient.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(quotientbegin + i) = (*(selfbegin + i))/(*(operandbegin+i));
      }
      return quotient;
    }
    //template <class T>
    array<CT,CI> operator/(CT operand){
      array<CT,CI> quotient(shape);
      
      CT* selfbegin = Begin(), *quotientbegin = quotient.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(quotientbegin + i) = (*(selfbegin + i))/operand;
      }
      return quotient;
    }
    
    //--------------------------------------------------
    //--------------------------------------------------
    //template <class OCT, class OCI>
    void operator*=(array<CT,CI> operand){
      CT* selfbegin = Begin(), *operandbegin = operand.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(selfbegin + i) = (*(selfbegin + i))*(*(operandbegin+i));
      }
    }
    template <class T>
    void operator*=(T operand){
      CT* selfbegin = Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(selfbegin + i) = (*(selfbegin + i))*operand;
      }
    }
    
    //--------------------------------------------------
    //template <class OCT, class OCI>
    void operator+=(array<CT,CI> operand){
      CT* selfbegin = Begin(), *operandbegin = operand.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(selfbegin + i) = (*(selfbegin + i))+(*(operandbegin+i));
      }
    }
    //template <class T>
    void operator+=(CT operand){
      CT* selfbegin = Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(selfbegin + i) = (*(selfbegin + i))+operand;
      }
    }
    
    //--------------------------------------------------
    //template <class OCT, class OCI>
    void operator-=(array<CT,CI> operand){
      CT* selfbegin = Begin(), *operandbegin = operand.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(selfbegin + i) = (*(selfbegin + i))-(*(operandbegin+i));
      }
    }
    //template <class T>
    void operator-=(CT operand){
      CT* selfbegin = Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(selfbegin + i) = (*(selfbegin + i))-operand;
      }
    }
    
    //--------------------------------------------------
    //template <class OCT, class OCI>
    void operator/=(array<CT,CI> operand){
      CT* selfbegin = Begin(), *operandbegin = operand.Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(selfbegin + i) = (*(selfbegin + i))/(*(operandbegin+i));
      }
    }
    //template <class T>
    void operator/=(CT operand){
      CT* selfbegin = Begin();
      for (CI i = 0; i < X[0]; ++i){
	*(selfbegin + i) = (*(selfbegin + i))/operand;
      }
    }

    //--------------------------------------------------
    //--------------------------------------------------
    //template <class I>
    array<CT,CI> operator[](int index){
      std::vector<CI> newshape(shape.size()-1);
      newshape.assign(shape.begin()+1,shape.end());
      array<CT,CI> newarray(newshape);
      newarray.Assign(origin+index*X[1], origin+(index+1)*X[1]-1);
      return newarray;
    }
  };
}
