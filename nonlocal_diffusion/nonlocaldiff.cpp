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
double u_exact(double x, double y,double t){
  return x*x+y*y+t*t;//cos(myB*pi*x)*sin(myB*pi*y);//t*cos(myB*pi*x)*sin(myB*pi*y);
}



double Ffun(double x,double y,double t){
  return 2*t-4;//4*Bulk*myB*pi*myB*pi*cos(myB*pi*x)*sin(myB*pi*y)/5;//4*Bulk*myB*pi*myB*pi*t*cos(myB*pi*x)*sin(myB*pi*y)/5;
}



////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////





int N = 80; 
const double delta=3.8; //ratio delta/h
double h = 1./N;            // lattice constant, all the nodes are evenly distributed 
int k = floor(delta);
int dim = (N+1+2*k)*(N+1+2*k);
//int basedim = 10;      //linear basis funtions
double* x_0 = new double [dim];    //initial configuration 
double* y_0 = new double [dim];
double* u = new double [dim];   //displacement in x direction and y direction 
double* u0 = new double [dim];
double* du0 = new double [dim];
double* rhs = new double [dim];
double* stiffK = new double [dim*dim];
double* stiffKT = new double[dim*dim];

vector< vector<int> > nei((N+1)*(N+1)); //neighborlist

double diff = 100;
const double totalT=0.1;
const int nstep = N*N/40;//time-step 
const double Dt = totalT/nstep;//0.0001; // time-stepsize 

int* IPIV = new int[dim];
int info;
int nrhs = 1;
char C='T';
char D = 'A';

void Preprocess()
{
  int cenx,cenx1;
  vvector<double> xi(2),zeta(2),gamma(2), eta(2);

  //set initial configuration 
  for(int i=0;i<N+1+2*k;i++){
    for(int j=0;j<N+1+2*k;j++){
      x_0[i*(N+1+2*k)+j] = (j-k)*h;
      y_0[i*(N+1+2*k)+j] = (i-k)*h;
    }
  }


  for(int i=0;i<N+1+2*k;i++){
    for(int j=0;j<N+1+2*k;j++){
      u0[i*(N+1+2*k)+j]=u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j],0.);
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

    for(int i=0;i<N+1;i++){
        for(int j=0;j<N+1;j++){
            cenx = (i+k)*(N+1+2*k)+j+k;
            
            // generate quadrature weights (double* omega)
            double* B = new double[6*(nei[i*(N+1)+j].size())];
            double* g = new double[6]; 
            for (int s = 0; s< nei[i*(N+1)+j].size(); s++){
                B[0*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*0.5*1; 
                B[1*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*0.5*(x_0[nei[i*(N+1)+j][s]]-x_0[cenx]); 
                B[2*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*0.5*(y_0[nei[i*(N+1)+j][s]]-y_0[cenx]); 
                B[3*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*0.5*pow(x_0[nei[i*(N+1)+j][s]]-x_0[cenx],2); 
                B[4*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*0.5*pow(y_0[nei[i*(N+1)+j][s]]-y_0[cenx],2); 
                B[5*nei[i*(N+1)+j].size()+s] = 8/pi/pow(delta*h,4)*0.5*(x_0[nei[i*(N+1)+j][s]]-x_0[cenx])*(y_0[nei[i*(N+1)+j][s]]-y_0[cenx]);
            }
            g[0] = 8/pi/pow(delta*h,4)*0.5*pi*(delta*h)*(delta*h);
            g[1] = 0; 
            g[2] = 0; 
            g[3] = 8/pi/pow(delta*h,4)*0.5*1/4*pi*pow(delta*h, 4);
            g[4] = 8/pi/pow(delta*h,4)*0.5*1/4*pi*pow(delta*h, 4);
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

            double* omega = new double[nei[i*(N+1)+j].size()];
            for (int a = 0; a<nei[i*(N+1)+j].size(); a++){
              omega[a]=0;
              for (int ss = 0; ss<6; ss++){
                  omega[a] += K[a*6+ss]*g[ss];
                  
                }
                //cout << omega[a] << endl; 
            }

            for(int s=0;s<nei[i*(N+1)+j].size();s++){
                stiffKT[cenx*dim+nei[i*(N+1)+j][s]] -= 8/pi/pow(delta*h,4)*omega[s];
                stiffKT[cenx*dim+cenx] += 8/pi/pow(delta*h,4)*omega[s];
            }
            stiffKT[cenx*dim+cenx] += 1./Dt;
            

        }
    }



  for(int i=0;i<N+1+2*k;i++){
    if(i < k || i > N+k){
      for(int j=0;j<N+1+2*k;j++){
        cenx = i*(N+1+2*k)+j;
        stiffKT[cenx*dim+cenx] = 1;
      }
    } 
    else{
      for(int j=0;j<k;j++){
        cenx = i*(N+1+2*k)+j;
        cenx1 = i*(N+1+2*k)+N+2*k-j;
        stiffKT[cenx*dim+cenx] = 1;
        stiffKT[cenx1*dim+cenx1] = 1;
      }
    }
  }

  //for(int i=0;i<dim;i++){
  //	for(int j=0;j<dim;j++){
	//	cout<<stiffKT[i*dim+j]<<" ";	
  //  	}
//	cout<<endl;
//  }

}

//void BE(int s,npArray g, npArray v, int UPDATE)
void BE(int s, int UPDATE)
{
  int cenx,cenx1;

/*  auto g_buffer = g.request();
  double *g_ptr = (double *)g_buffer.ptr;
  auto v_buffer = v.request();
  double *v_ptr = (double *)v_buffer.ptr;*/

  //apply boundary condition
  for(int i=0;i<N+1;i++){
    for(int j=0;j<N+1;j++){
      cenx = (i+k)*(N+1+2*k)+j+k;
      rhs[cenx] = Ffun(x_0[cenx],y_0[cenx],s*Dt)+u0[cenx]/Dt;
    }
  }


  //Dirichlet boundary condition 
  for(int i=0;i<N+1+2*k;i++){
    if(i < k || i > N+k){
      for(int j=0;j<N+1+2*k;j++){
        rhs[i*(N+1+2*k)+j] = u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j],s*Dt);
      }
    }
    else{
      for(int j=0;j<k;j++){
        rhs[i*(N+1+2*k)+j] = u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j],s*Dt);
        rhs[i*(N+1+2*k)+N+2*k-j] = u_exact(x_0[i*(N+1+2*k)+N+2*k-j],y_0[i*(N+1+2*k)+N+2*k-j],s*Dt);
      }
    }
  }




  for(int i=0;i<dim*dim;i++){
    stiffK[i] = stiffKT[i];
  }



  dgetrf_(&dim,&dim,stiffK,&dim,IPIV,&info);
  dgetrs_(&C,&dim,&nrhs,stiffK,&dim, IPIV,rhs,&dim,&info);

  for(int i=0;i<N+1;i++){
    for(int j=0;j<N+1;j++){
      cenx = (i+k)*(N+1+2*k)+j+k;
      du0[cenx] = (rhs[cenx]-u0[cenx])/Dt;
    }
  }

  if(UPDATE)
  {
    for(int i=0;i<dim;i++){
      u0[i] = rhs[i];
    }
    }

  double error=0.;
  double errorinf=0.;
  if(UPDATE)
  {
      for(int i=0;i<N+1+2*k;i++){
        for(int j=0;j<N+1+2*k;j++){
          error += pow(u0[i*(N+1+2*k)+j]-u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j],s*Dt),2);
        
        //cout<< x_0[i*(N+1+2*k)+j] <<" "<<  y_0[i*(N+1+2*k)+j] <<" "<<pow(u0[i*(N+1+2*k)+j]-u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j],s*Dt),2)<<endl;
        }
      }
      for(int i=0;i<N+1+2*k;i++){
        for(int j=0;j<N+1+2*k;j++){
          if(fabs(u0[i*(N+1+2*k)+j]-u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j],s*Dt)) > errorinf){
            errorinf = fabs(u0[i*(N+1+2*k)+j]-u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j],s*Dt)); 
          }
        }
      }

      error /= (N+1)*(N+1);
      error = sqrt(error);
  }
  cout<<"timestep: "<<s<<" L2 error: "<<error<<"   errorinf: "<<errorinf<<endl;

}


//void Newark(int s,npArray g, npArray v, int UPDATE)



int main(int argc, char* argv[]){
      /*  if(argc < 3){
          cerr<<"Error: Not enough input variables"<<" ";
          cerr<<" Usage: <N>  <delta>"<<endl;
         }


          int N = atoi(argv[1]);                        
          double delta = atof(argv[2]);*/
      //      double *g = new double [2*N+2];
      int totaliter=1;
      Preprocess();

      for(int s=0;s<nstep;s++){
        for(int subiter=0; subiter<totaliter; subiter++)
        {
          if (subiter<(totaliter-1))
            BE(s+1,0);
          else
            BE(s+1,1);
        }
      }
      double error = 0.;
      for(int i=0;i<N+1+2*k;i++){
        for(int j=0;j<N+1+2*k;j++){
          error += pow(u0[i*(N+1+2*k)+j]-u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j],nstep*Dt),2);
          //  error += pow(u[i*(N+1+2*k)+j]-pow(x[i*(N+1+2*k)+j],2),2) + pow(u[i*(N+1+2*k)+j+(N+1+2*k)*(N+1+2*k)]-pow(x[i*(N+1+2*k)+j+(N+1+2*k)*(N+1+2*k)],2),2);
          //error += pow(u[i*(N+1+2*k)+j] - pow(x[i*(N+1+2*k)+j],2),2);
        }
      }
      double errorinf = 0;
      for(int i=0;i<N+1+2*k;i++){
        for(int j=0;j<N+1+2*k;j++){
          if(fabs(u0[i*(N+1+2*k)+j]-u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j],nstep*Dt)) > errorinf){
            errorinf = fabs(u0[i*(N+1+2*k)+j]-u_exact(x_0[i*(N+1+2*k)+j],y_0[i*(N+1+2*k)+j],nstep*Dt)); 
          }
        }
      }


      error /= (N+1)*(N+1);
      error = sqrt(error);
      cout<<"AT final step T="<<totalT<<" error L2: "<<error<<"  errorinf: "<<errorinf<<endl;
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
