#ifndef VVECTOR_H
#define VVECTOR_H

#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

template<class Type>
class vvector : public std::vector<Type>
{
    private:
        typedef std::vector<Type> Base;
    public:
        vvector() { }

        vvector(int size, Type val = 0) : Base(size, val) { }

		template<class C>
		vvector& operator=( vvector<C>& vv){
			if((*this).size() != vv.size()){
				std::cerr<<"Different dim"<<std::endl;
				return *this;
			}
			for(size_t i=0;i!=(*this).size();++i){
				(*this)[i]=vv[i];
			}
			return *this;
		}


		vvector& operator+=(const vvector& vv) { 
			for(size_t i=0;i!=(*this).size();++i){
				(*this)[i]+=vv[i];
			}
			return *this; 
		}

		vvector& operator-=(const vvector& vv) { 
			for(size_t i=0;i!=(*this).size();++i){
				(*this)[i]-=vv[i];
			}
			return *this; 
		}

		template<class C> vvector& operator*=(const C& CC) { 
			for(size_t i=0;i!=(*this).size();++i){
				(*this)[i]*=CC; 
			}
			return *this; 
		}

		template<class C> vvector& operator/=(const C& CC) { 
			for(size_t i=0;i!=(*this).size();++i){
				(*this)[i]/=CC; 
			}
			return *this; 
		}
        
		Type Norm1(){
			Type a(0);
			for(size_t i=0;i!=(*this).size();++i){
				a+=std::abs((*this)[i]);
			}
			return a;
		}

		Type Norm2(){
			Type a(0);
			for(size_t i=0;i!=(*this).size();++i){
				a+=pow((*this)[i],2);
			}
			return sqrt(a);
		}

		Type NormInfinite(){
			Type a(0);
			for(size_t i=0;i!=(*this).size();++i){
				if(std::abs((*this)[i])>a){
					a=std::abs((*this)[i]);
				}
			}
			return a;
		}


//		Type arrmultiply(const vvector<C>& c){
//			Type m(0);
//			for(size_t i=0;i!=(*this).size();++i){
//				m+=(*this)[i]*c[i];
//			}
//			return m;
//		}

		template<class C> 
		vvector<int> sign(){
			vvector<int> s((*this).size());
			for(size_t i=0;i!=(*this).size();++i){
				if((*this)[i]>0) {
					s[i]=1;
				}
				else s[i]=-1;
			}
			return s;
		}
		
		Type Householder(vvector<Type>& v){//x的Householder变换,返回beta
			int n=(*this).size(); 
			Type eta=NormInfinite();
			if(eta==0){
				return 0;
			}
			(*this)/=eta;
			Type sigma(0);
			for(int i=1;i!=n;++i){
				sigma+=pow((*this)[i],2);
			}
			v=(*this);
			if(std::abs(sigma)>pow(10,-14)){
				Type alpha(sqrt(sigma+pow((*this)[0],2)));
				if((*this)[0] <= 0) { v[0]=(*this)[0]-alpha; }
				else v[0]=-sigma/((*this)[0]+alpha);
				Type beta=2*pow(v[0],2)/(sigma+pow(v[0],2));
				Type temp=v[0];
				v/=temp;
				return beta;
			}
			else{
				return 0;
			}

		}

        template<class T> friend std::ostream& operator<<(std::ostream&, const vvector&);

        template<class T> friend std::istream& operator>>(std::istream&, vvector&);
};

template<class T> 
std::istream& operator>>(std::istream& in, vvector<T>& vv){
	std::cout<<"Please input numbers"<<std::endl;
	T value;
	while(in>>value){
		vv.push_back(value);
	}
	return in;
}

template<class T> 
std::ostream& operator<<(std::ostream& os, const vvector<T>& vv){
	os.precision(4);
	os<<std::showpos;
	os.setf(std::ios::scientific);
	for(size_t i=0;i!=vv.size();++i){
		os<<std::setw(15)<<vv[i];
	}

	return os;
}

template<class T>
vvector<T> operator+(const vvector<T>& v1, const vvector<T>& v2)
{
	if(v1.size() != v2.size()){
		std::cerr<<"diferent dim"<<std::endl;
		return v1;
	}
    vvector<T> v3(v1);
    v3 += v2;
    return v3;
}

template<class T>
vvector<T> operator-(const vvector<T>& v1, const vvector<T>& v2)
{
	if(v1.size() != v2.size()){
		std::cerr<<"diferent dim"<<std::endl;
		return v1;
	}
    vvector<T> v3(v1);
    v3 -= v2;
    return v3;
}

template<class T>
vvector<T> operator*(const vvector<T>& v1, const vvector<T>& v2)
{
	if(v1.size() != 3 || v2.size() != 3){
		std::cerr<<"Size don't match"<<std::endl;
		return v1;
	}
	vvector<T> v3(3);
	v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
	v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
	v3[2] = v1[0]*v2[1] - v1[1]*v2[0];
	return v3;
}


template<class T,class T1>
vvector<T> operator*(const vvector<T>& v1, const T1& a)
{
    vvector<T> v3(v1);
    v3 *= a;
    return v3;
}

template<class T,class T1>
vvector<T> operator*(const T1& a, const vvector<T>& v1)
{
    vvector<T> v3(v1);
    v3 *= a;
    return v3;
}

template<class T,class T1>
vvector<T> operator/(const vvector<T>& v1, const T1& a)
{
	if(a==0){
		std::cerr<<"Cannot be divided by 0."<<std::endl;
		return v1;
	}
    vvector<T> v3(v1);
    v3 /= a;
    return v3;
}

#endif //VVECTOR_H 
