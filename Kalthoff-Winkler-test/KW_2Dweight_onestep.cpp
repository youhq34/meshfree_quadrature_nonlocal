#include<iostream>
#include<vector>
#include<cmath>
#include<cstdlib>
#include"vvector.h"

extern "C" void dgetrf_(int*, int*, double*, int*, int*, int*);
extern "C" void dgetri_(int* , double* , int* , int* , double* , int* , int* );
extern "C" void dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);
extern "C" void dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);


using namespace std;
const double pi = 3.141592653;


//basis funtions, use second-order basis here
vvector<double> phi(int m,vvector<double> x){
    vvector<double> y(2);
    if(m == 0){
        y[0] = x[0];
        y[1] = 0;
        return y;
    }
    else if(m == 1){
        y[0] = 0;
        y[1] = x[0];
        return y;
    }
    else if(m == 2){
        y[0] = x[1];
        y[1] = 0;
        return y;
    }
    else if(m == 3){
        y[0] = 0;
        y[1] = x[1];
        return y;
    }
    else if(m == 4){
        y[0] = x[0]*x[0];
        y[1] = 0;
        return y;
    }
    else if(m == 5){
        y[0] = 0;
        y[1] = x[0]*x[0];
        return y;
    }
    else if(m == 6){
        y[0] = x[1]*x[1];
        y[1] = 0;
        return y;
    }
    else if(m == 7){
        y[0] = 0;
        y[1] = x[1]*x[1];
        return y;
    }
    else if(m == 8){
        y[0] = x[0]*x[1];
        y[1] = 0;
        return y;
    }
    else{
        y[0] = 0;
        y[1] = x[0]*x[1];
        return y;
    }
}

double u_exact(double x, double y){
    return sin(x)*cos(y);
}
double v_exact(double x, double y){
    return cos(x)*sin(y);
}

double fx_exact(double x, double y, double Bulk){
    return -18*Bulk/5*sin(x)*cos(y);//-3*Bulk/5*cos(x)*sin(y);
}
double fy_exact(double x, double y, double Bulk){
    return -18*Bulk/5*cos(x)*sin(y);//-3*Bulk*cos(x)*sin(y);
}


int main(int argc, char* argv[]){
    if(argc < 4){
        cerr<<"Error: Not enough input variables"<<" ";
        cerr<<" Usage: <N>  <delta> <random_perturbe>"<<endl;
    }

    
    int N = atoi(argv[1]);                        
    double delta = atof(argv[2]);
    
    double Mx = 20.0;
    double My = 10.0;               
    double h = My/(N-1);            // lattice constant, all the nodes are evenly distributed 
    int k = floor(delta);
    int dim = 2*N*(2*N-1);
    cout<<"total number of points: "<<dim/2<<endl;
    int basedim = 10;      //linear basis funtions
    double* x = new double [dim/2];    //initial configuration 
    double* y = new double [dim/2];    //initial configuration 
    double* u = new double [dim];   //displacement in x direction and y direction 
    double* rhs = new double [dim];
    double* stiffK = new double [dim*dim];
    vector< vector<int> > nei(N*(2*N-1)); //neighborlist
    vector< vector<double> > weight(N*(2*N-1)); //quadrature weights
    vvector<double> xi(2),zeta(2),gamma(2), eta(2);

    double diff = 100;
    double error = 0;
    double load = 100;
    double Bulk = 4e4; 
    // constant in 2d PMB model 
    double c = 72*Bulk/5/pi/pow(delta*h,3);
    double slope = 5*load/3/Bulk;
    int cenx,ceny,cenx1,ceny1;
    double dpi = 0;
    int nstep = 0;
    int s = 0;

    int* IPIV = new int[dim];
    int info;
    int nrhs = 1;
    char C='T';
    char D = 'A';

    //set initial configuration
    cout<<"set initial configuration"<<endl;
    double movex = -10.0;
    double movey = -10.0;
    double centerx = 0.0;
    double centery = -5.0;
    vector<int> ID(N*(2*N-1));
    for(int i=0;i<2*N-1;i++){
        for(int j=0;j<N;j++){
            x[i*N+j] = h*i-10.0;
            y[i*N+j] = h*j-10.0;
            ID[i*N+j] = 0;
            if(abs(x[i*N+j]-centerx) > 10-delta*h) ID[i*N+j] = 1;
            if(abs(y[i*N+j]-centery) > 5-delta*h) ID[i*N+j] = 1;
        }
    }

    
    //random perturbation 
    /*
    double random_coe = atof(argv[3]);
    for(int i=0;i<N+1+2*k;i++){
        for(int j=0;j<N+1+2*k;j++){
            x[i*(N+1+2*k)+j] += (2*rand()/double(RAND_MAX)-1)*random_coe*h;
            x[i*(N+1+2*k)+j+(N+1+2*k)*(N+1+2*k)] += (2*rand()/double(RAND_MAX)-1)*random_coe*h ;
        }
    }
    */
    //set neighborhood list for every particle 
    cout<<"build neighborhood list"<<endl;
    for(int i=0;i<2*N-1;i++){
        for(int j=0;j<N;j++){
            if(ID[i*N+j] == 0){
                for(int s1=-k;s1<k+1;s1++){
                    for(int s2=-k;s2<k+1;s2++){
                        xi[0] = x[(i+s1)*N+j+s2] - x[i*N+j];
                        xi[1] =  y[(i+s1)*N+j+s2] - y[i*N+j];
                        if(xi.Norm2() < delta*h && !( s1==0 && s2==0)){
                            nei[i*N+j].push_back((i+s1)*(N)+j+s2);
                        }
                    }
                }
            }       
        }
    }
    cout<<"generate quadrature weights"<<endl;
    //solve for quadrature weights
    double* WORK = new double[1000];
    for(int i=0;i<2*N-1;i++){
        for(int j=0;j<N;j++){
            if(ID[i*N+j] == 0){
            s = nei[i*N+j].size();
            weight[i*N+j].resize(s,0);
            double* stiffweight = new double[(s+2*basedim)*(s+2*basedim)];
            double* rhsweight = new double[s+2*basedim];
            double* U = new double[(s+2*basedim)*(s+2*basedim)];
            double* VT = new double[(s+2*basedim)*(s+2*basedim)];
            double* S = new double[s+2*basedim];
            int INFO;
            int LWORK = 1000;
            int dim1 = s+2*basedim;
            for(int t=0;t<(s+2*basedim)*(s+2*basedim);t++){
                stiffweight[t] = 0;
            }
            for(int t=0;t<s+2*basedim;t++){
                rhsweight[t] = 0;
            }
            for(int t=0;t<s;t++){
                zeta[0] = x[i*N+j];
                zeta[1] = y[i*N+j];
                gamma[0] = x[nei[i*N+j][t]];
                gamma[1] = y[nei[i*N+j][t]];
                xi[0] = gamma[0] - zeta[0]; 
                xi[1] = gamma[1] - zeta[1];
                stiffweight[t*(s+2*basedim)+t] = 2;
                for(int m=0;m<basedim;m++){
                    eta = phi(m,gamma) - phi(m,zeta);
                    stiffweight[t*(s+2*basedim)+s+m] = c*(xi[0]*xi[0]*eta[0] + xi[0]*xi[1]*eta[1])/pow(xi.Norm2(),3);
                    stiffweight[t*(s+2*basedim)+s+basedim+m] = c*(xi[0]*xi[1]*eta[0] + xi[1]*xi[1]*eta[1])/pow(xi.Norm2(),3);
                    stiffweight[(s+m)*(s+2*basedim)+t] = stiffweight[t*(s+2*basedim)+s+m];
                    stiffweight[(s+basedim+m)*(s+2*basedim)+t] = stiffweight[t*(s+2*basedim)+s+basedim+m];
                }
                
            }
            
            rhsweight[s+4] = c*pi*pow(delta*h,3)/4;
            rhsweight[s+6] = c*pi*pow(delta*h,3)/12;
            rhsweight[s+9] = c*pi*pow(delta*h,3)/12;
            rhsweight[s+10+5] = c*pi*pow(delta*h,3)/12;
            rhsweight[s+10+7] = c*pi*pow(delta*h,3)/4;
            rhsweight[s+10+8] = c*pi*pow(delta*h,3)/12;
            
            dgesvd_(&D,&D,&dim1,&dim1,stiffweight,&dim1,S,U,&dim1,VT,&dim1,WORK,&LWORK,&INFO);
            for(int m=0;m<dim1;m++){
                if(S[m] < 1e-6) S[m] = 0;
            }
            for(int m=0;m<dim1;m++){
                if(S[m] != 0) S[m] = 1./S[m];
            }
            
            
            for(int t=0;t<s;t++){
                for(int l=0;l<dim1;l++){
                    for(int m=0;m<dim1;m++){
                            weight[i*N+j][t] += VT[t*dim1+l]*S[l]*U[l*dim1+m]*rhsweight[m];
                    }
                }
            }
            delete[] stiffweight;
            delete[] rhsweight;
            delete[] U;
            delete[] S;
            delete[] VT;

            }
        }
    }
    
            delete[] WORK;
        
        for(int i=0;i<dim;i++){
            for(int j=0;j<dim;j++){
                stiffK[i*dim+j] = 0;
            }
            rhs[i] = 0; 
        }
        
        cout<<"assemble matrix"<<endl; 
        for(int i=0;i<2*N-1;i++){
            for(int j=0;j<N;j++){
                cenx = i*N+j;
                ceny = i*N+j+(2*N-1)*N;
                //cout<<x[cenx]<<" "<<y[cenx]<<endl;
                if(ID[cenx] == 0){
                    //cout<<nei[cenx].size()<<endl;
                    for(int s=0;s<nei[cenx].size();s++){
                        //cout<<weight[cenx][s]<<" ";
                        xi[0] = x[nei[cenx][s]] - x[cenx];
                        xi[1] = y[nei[cenx][s]] - y[cenx];
                        stiffK[cenx*dim+cenx] += -c*xi[0]*xi[0]*weight[cenx][s]/pow(xi.Norm2(),3);
                        stiffK[cenx*dim+nei[cenx][s]] += c*xi[0]*xi[0]*weight[cenx][s]/pow(xi.Norm2(),3);
                        stiffK[cenx*dim+ceny] += -c*xi[0]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3);
                        stiffK[cenx*dim+nei[cenx][s]+dim/2] += c*xi[0]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3);

                        stiffK[ceny*dim+ceny] += -c*xi[1]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3);
                        stiffK[ceny*dim+nei[cenx][s]+dim/2] += c*xi[1]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3);
                        stiffK[ceny*dim+cenx] += -c*xi[0]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3);
                        stiffK[ceny*dim+nei[cenx][s]] += c*xi[0]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3);
                    }
                    //cout<<endl;
                    rhs[cenx] = fx_exact(x[cenx],y[cenx],Bulk);
                    rhs[ceny] = fy_exact(x[cenx],y[cenx],Bulk);//+9*Bulk*pow(delta*h,2)/25;
                }
                else{
                    rhs[cenx] = u_exact(x[cenx],y[cenx]);
                    rhs[ceny] = v_exact(x[cenx],y[cenx]);
                    stiffK[cenx*dim+cenx] = 1;
                    stiffK[ceny*dim+ceny] = 1; 
                }
                   
            }
        }
    /*
    ofstream mfile;
    mfile.open("matrix.txt");
    for(int i=0;i<dim;i++){
            for(int j=0;j<dim;j++){
                //u[i] += stiffK[i*dim+j]*rhs[j];
            mfile<<stiffK[i*dim+j]<<" ";
            }
            mfile<<endl;
        }
mfile.close();
*/
  //      for(int i=0;i<N+1+2*k;i++){
  //          for(int j=0;j<N+1+2*k;j++){
  //              cout<<u[i*(N+1+2*k)+j]<<" ";
  //          }
  //          cout<<endl;
  //      }
        cout<<"solve for displacement"<<endl;
        dgetrf_(&dim,&dim,stiffK,&dim,IPIV,&info);
        dgetrs_(&C,&dim,&nrhs,stiffK,&dim, IPIV,rhs,&dim,&info);
        
        diff = 0;
       



        for(int i=0;i<dim;i++){
            u[i] = rhs[i];
        }    
        

    for(int i=0;i<2*N-1;i++){
        for(int j=0;j<N;j++){
            cenx = i*N+j;
            ceny = i*N+j+(2*N-1)*N;
            error += pow( u[cenx] - u_exact(x[cenx],y[cenx]),2) + pow( u[ceny] - v_exact(x[cenx],y[cenx]),2); 
        }
    }

    
    error /= N*(2*N-1);
    error = sqrt(error);

    cout<<" error: "<<error<<endl;
    ofstream ufile,ufile1;
    ufile.open("displacement_y.csv");
    ufile1.open("displacement_x.csv");
    ufile<<"x y v"<<endl;
    ufile1<<"x y u"<<endl;
    for(int i=0;i<2*N-1;i++){
        for(int j=0;j<N;j++){
            cenx = i*N+j;
            ceny = i*N+j+(2*N-1)*N;
            ufile1<<x[cenx]<<" "<<y[cenx]<<" "<<u[cenx]<<" "<<endl;
            ufile<<x[cenx]<<" "<<y[cenx]<<" "<<u[ceny]<<" "<<endl;
        }
        ufile1<<endl;
        ufile<<endl;
    }
    
    delete[] x;
    delete[] u;
    delete[] rhs;
    delete[] stiffK;
    delete[] IPIV;
    return 0;  
}
