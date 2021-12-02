#include<iostream>
#include<vector>
#include<cmath>
#include<cstdlib>
#include"vvector.h"
#include <iomanip>
#include <fstream>


extern "C" void dgetrf_(int*, int*, double*, int*, int*, int*);
extern "C" void dgetri_(int* , double* , int* , int* , double* , int* , int* );
extern "C" void dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);
extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);

using namespace std;
const double pi = 3.1415926535897932384626;

double weight(double r,double epsilon){
    if(r > epsilon) return 0;
    else return pow((1-r/epsilon),4);
}

void inverse(double* A, int N)
{
    int *IPIV = new int[N+1];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;
    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
    delete[] IPIV;
    delete[] WORK;
}

//basis funtions
//
/////////////////////////////////////////////////////////////////////////
//case 1: quadratic in space, static in time
double u_exact(double x, double y){
  //return x*x+y*y+t*t;//cos(myB*pi*x)*sin(myB*pi*y);//t*cos(myB*pi*x)*sin(myB*pi*y);
  return cos(pi*x)*cos(pi*y);  
  //return x*x+y*y;
}



double Ffun(double x,double y){
  //return 2*t-4;//4*Bulk*myB*pi*myB*pi*cos(myB*pi*x)*sin(myB*pi*y)/5;//4*Bulk*myB*pi*myB*pi*t*cos(myB*pi*x)*sin(myB*pi*y)/5;
  return 2*pi*pi*cos(pi*x)*cos(pi*y);  
  //return -4;
}



////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////








vector<vector<double>>  Preprocess(double* x_0, double* y_0, int N,int k, double delta,vector<vector<int>> nei)
{
  int cenx,cenx1;
  vvector<double> xi(2),zeta(2),gamma(2), eta(2);
  vector<vector<double>> omega((N+1)*(N+1));  
    double h = 1./N;
    for(int i=0;i<N+1;i++){
        for(int j=0;j<N+1;j++){
            cenx = (i+k)*(N+1+2*k)+j+k;
            omega[i*(N+1)+j].resize(nei[i*(N+1)+j].size());
            // generate quadrature weights (double* omega)
            double* B = new double[6*(nei[i*(N+1)+j].size())];
            double* g = new double[6]; 
            for (int s = 0; s< nei[i*(N+1)+j].size(); s++){
                B[0*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*1; 
                B[1*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*(x_0[nei[i*(N+1)+j][s]]-x_0[cenx]); 
                B[2*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*(y_0[nei[i*(N+1)+j][s]]-y_0[cenx]); 
                B[3*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*pow(x_0[nei[i*(N+1)+j][s]]-x_0[cenx],2); 
                B[4*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*pow(y_0[nei[i*(N+1)+j][s]]-y_0[cenx],2); 
                B[5*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*(x_0[nei[i*(N+1)+j][s]]-x_0[cenx])*(y_0[nei[i*(N+1)+j][s]]-y_0[cenx]);
            }
            g[0] = 8/pi/pow(delta*h,4)*pi*(delta*h)*(delta*h);
            g[1] = 0; 
            g[2] = 0; 
            g[3] = 8/pi/pow(delta*h,4)*1/4*pi*pow(delta*h, 4);
            g[4] = 8/pi/pow(delta*h,4)*1/4*pi*pow(delta*h, 4);
            g[5] = 0;   
            
            double* S = new double[6*6]; 
            
            for (int a = 0; a<6; a++){
              for (int b = 0; b<6; b++){
                S[a*6+b]=0;
                for (int ss = 0; ss<nei[i*(N+1)+j].size(); ss++){
                  S[a*6+b] += B[a*nei[i*(N+1)+j].size()+ss]*B[b*nei[i*(N+1)+j].size()+ss];
                }
              }
            }

            //pseudo_inverse(S, 6); 
            inverse(S,6);

            double* K = new double[nei[i*(N+1)+j].size()*6]; 
    
            for (int a = 0; a<nei[i*(N+1)+j].size(); a++){
              for (int b = 0; b<6; b++){
                K[a*6+b]=0;
                for (int ss = 0; ss<6; ss++){
                  K[a*6+b] += B[ss*nei[i*(N+1)+j].size()+a]*S[ss*6+b];
                }
              }
            }

            for (int a = 0; a<nei[i*(N+1)+j].size(); a++){
              omega[i*(N+1)+j][a]=0;
              for (int ss = 0; ss<6; ss++){
                  omega[i*(N+1)+j][a] += K[a*6+ss]*g[ss];  
                }
            }
            

        }
    }
    return omega;

}





int main(int argc, char* argv[]){
    if(argc < 3){
        cerr<<"Error: Not enough input variables"<<" ";
        cerr<<" Usage: <N>  <delta>"<<endl;
     }


    int N = atoi(argv[1]);                        
    double delta = atof(argv[2]);
      //      double *g = new double [2*N+2];
    double h = 1./N;            // lattice constant, all the nodes are evenly distributed 
    int k = floor(delta);
    int dim = (N+1+2*k)*(N+1+2*k);
    //int basedim = 10;      //linear basis funtions
    double* x_0 = new double [dim];    //initial configuration 
    double* y_0 = new double [dim];
    double* u = new double [dim];   //displacement in x direction and y direction 
    double* u0 = new double [dim];
    double* rhs = new double [dim];
    double* stiffKT = new double[dim*dim];
    int* IPIV = new int[dim];
    int info;
    int nrhs = 1;
    char C='T';
    char D = 'A';

    vector< vector<int> > nei((N+1)*(N+1)); //neighborlist
  vvector<double> xi(2),zeta(2),gamma(2), eta(2);
    
    //Build the neighborhood list 
  //set initial configuration 
  for(int i=0;i<N+1+2*k;i++){
    for(int j=0;j<N+1+2*k;j++){
      x_0[i*(N+1+2*k)+j] = (j-k)*h;
      y_0[i*(N+1+2*k)+j] = (i-k)*h;
    }
  }


  for(int i=0;i<N+1+2*k;i++){
    for(int j=0;j<N+1+2*k;j++){
      u0[i*(N+1+2*k)+j]=u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j]);
    }
  }


  for(int i=0;i<dim*dim;i++){
    stiffKT[i] = 0;
  }

  //set neighborhood list for every particle 
  for(int i=0;i<N+1;i++){
    for(int j=0;j<N+1;j++){
      for(int s1=-k;s1<k+1;s1++){
        for(int s2=-k;s2<k+1;s2++){
          if((i+k+s1)*(N+1+2*k)+j+k+s2 < (N+1+2*k)*(N+1+2*k) && (i+k+s1)*(N+1+2*k)+j+k+s2 > -1){
            xi[0] = x_0[(i+k+s1)*(N+1+2*k)+j+k+s2] - x_0[(i+k)*(N+1+2*k)+j+k];
            xi[1] = y_0[(i+s1+k)*(N+1+2*k)+j+k+s2] - y_0[(i+k)*(N+1+2*k)+j+k];
            if(xi.Norm2() < delta*h){
              nei[i*(N+1)+j].push_back((i+k+s1)*(N+1+2*k)+j+k+s2);
            }
          }
        }
      }       
    }
  }

    
    
    //generate quadrature weights 
    vector<vector<double>> omega = Preprocess(x_0,y_0,N,k,delta,nei);

    //assemble matrix
    for(int i=0;i<N+1;i++){
        for(int j=0;j<N+1;j++){
            int cenx = (i+k)*(N+1+2*k)+j+k;
            for(int s=0;s<nei[i*(N+1)+j].size();s++){
                stiffKT[cenx*dim+nei[i*(N+1)+j][s]] -= 8/pi/pow(delta*h,4)*omega[i*(N+1)+j][s];
                stiffKT[cenx*dim+cenx] += 8/pi/pow(delta*h,4)*omega[i*(N+1)+j][s];
            }
            rhs[cenx] = Ffun(x_0[cenx],y_0[cenx]);//+u0[cenx]/Dt;

        }
    }  
  
    //apply boundary condition 
    for(int i=0;i<N+1+2*k;i++){
    if(i < k || i > N+k){
      for(int j=0;j<N+1+2*k;j++){
        int cenx = i*(N+1+2*k)+j;
        stiffKT[cenx*dim+cenx] = 1;
        rhs[cenx] = u_exact(x_0[cenx],y_0[cenx]);
      }
    } 
    else{
      for(int j=0;j<k;j++){
        int cenx = i*(N+1+2*k)+j;
        int cenx1 = i*(N+1+2*k)+N+2*k-j;
        stiffKT[cenx*dim+cenx] = 1;
        stiffKT[cenx1*dim+cenx1] = 1;
        rhs[cenx] = u_exact(x_0[cenx],y_0[cenx]);
        rhs[cenx1] = u_exact(x_0[cenx1],y_0[cenx1]);
      }
    }
  }
  dgetrf_(&dim,&dim,stiffKT,&dim,IPIV,&info);
  dgetrs_(&C,&dim,&nrhs,stiffKT,&dim, IPIV,rhs,&dim,&info);
  
  for(int i=0;i<dim;i++){
      u0[i] = rhs[i];
  }
    
    //post process, calculate error
    double error = 0.;
      for(int i=0;i<N+1+2*k;i++){
        for(int j=0;j<N+1+2*k;j++){
          error += pow(u0[i*(N+1+2*k)+j]-u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j]),2);
        }
      }
      double errorinf = 0;
      for(int i=0;i<N+1+2*k;i++){
        for(int j=0;j<N+1+2*k;j++){
          if(fabs(u0[i*(N+1+2*k)+j]-u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j])) > errorinf){
            errorinf = fabs(u0[i*(N+1+2*k)+j]-u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j])); 
          }
        }
      }


      error /= (N+1)*(N+1);
    error = sqrt(error);
      cout<<" error L2: "<<error<<"  errorinf: "<<errorinf<<endl;
      //ofstream ufile,ufile1;
      //ufile1.open("vel_x");
      //for(int i=0;i<N+1;i++){
      //  int j=N;
      //  int cenx = (i+k)*(N+1+2*k)+j+k;
      //  ufile1<<du0[cenx]<<" ";
     // }

      //ufile1.close();
      //ufile.close();

      return 0;  
}
