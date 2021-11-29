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
    double* damage = new double [dim/2];     
    double* u = new double [dim];   //displacement in x direction and y direction 
    double* u_old = new double [dim];   //displacement in x direction and y direction 
    double* u_oldold = new double [dim];   //displacement in x direction and y direction 
    double* rhs = new double [dim];
    double* stiffK = new double [dim*dim];
    vector< vector<int> > nei(N*(2*N-1)); //neighborlist
    vector< vector<int> > Bondbroke(N*(2*N-1)); //bondbroke
    vector< vector<double> > weight(N*(2*N-1)); //quadrature weights
    vvector<double> xi(2),zeta(2),gamma(2), eta(2);

    double diff = 100;
    double error = 0;
    double load = 100;
    double Bulk = 191e4;
    double smax = 0.0099/sqrt(delta*h);
    int Nsteps = 400;
    double dt = 0.001;
    double den = 8e-3;
    // constant in 2d PMB model 
    double c = 72*Bulk/5/pi/pow(delta*h,3);
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
            x[i*N+j] = h*i+movex;
            y[i*N+j] = h*j+movey;
            ID[i*N+j] = 0;
            if(abs(x[i*N+j]-centerx) > 10-delta*h) ID[i*N+j] = 1;
            if(abs(y[i*N+j]-centery) > 5-delta*h) ID[i*N+j] = 1;
            if(ID[i*N+j] == 1 && abs(x[i*N+j]) <= 2.5+0.725 && y[i*N+j] > -delta*h-1e-10) ID[i*N+j] = 2;
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
                            Bondbroke[i*N+j].push_back(false);
                        }
                    }
                }
            }       
        }
    }

    cout<<"break bonds that pass through boundary or crack"<<endl;
    for(int i=0;i<2*N-1;i++){
        for(int j=0;j<N;j++){
            if(ID[i*N+j] == 0){
                for(int s=0;s<nei[i*N+j].size();s++){
                    //break the left crack
                    double xc1 = -2.5-0.725;
                    if((x[i*N+j]-xc1)*(x[nei[i*N+j][s]]-xc1) <= 0 ){
                        double slope = (xc1-x[i*N+j])/(x[nei[i*N+j][s]]-x[i*N+j]);
                        double yc = y[i*N+j] + slope*(y[nei[i*N+j][s]]-y[i*N+j]);
                        if(yc >= -5.0) Bondbroke[i*N+j][s] = true;
                    }
                    //right crack    
                    double xc2 = 2.5+0.725;
                    if((x[i*N+j]-xc2)*(x[nei[i*N+j][s]]-xc2) <= 0 ){
                        double slope = (xc2-x[i*N+j])/(x[nei[i*N+j][s]]-x[i*N+j]);
                        double yc = y[i*N+j] + slope*(y[nei[i*N+j][s]]-y[i*N+j]);
                        if(yc >= -5.0) Bondbroke[i*N+j][s] = true;
                    }
                    //break other boundarys as free surface
                    if((y[i*N+j]<-2*delta*h) && (ID[nei[i*N+j][s]] == 1)){
                        Bondbroke[i*N+j][s] = true;
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
        
        
             
    ofstream ufile,ufile1,dfile;
            
        for(int step = 1;step <=Nsteps;step++ ){    
            if(step == 1){
                //print out initial damage file 
                //postprocess damage 
                for(int i=0;i<2*N-1;i++){
                    for(int j=0;j<N;j++){
                        int nneigh = 0;
                        double sum = 0;
                        if(ID[i*N+j] == 0){
                            for(int s=0;s<nei[i*N+j].size();s++){
                                if(Bondbroke[i*N+j][s]) sum++;
                                nneigh++;
                            }
                            damage[i*N+j] = sum/double(nneigh);
                        }
                    }
                }
                dfile.open("initial_damage.csv");
                dfile<<"x y damage"<<endl;
                for(int i=0;i<2*N-1;i++){
                    for(int j=0;j<N;j++){
                        if(ID[i*N+j] == 0){
                            dfile<<x[i*N+j]<<" "<<y[i*N+j]<<" "<<damage[i*N+j]<<endl;
                        }
                    }
                }
                dfile.close();
            }    
            double t = dt*step;
           
            
            
            
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
                        double broke = double(!Bondbroke[i*N+j][s]);
                        stiffK[cenx*dim+cenx] += c*xi[0]*xi[0]*weight[cenx][s]/pow(xi.Norm2(),3)*broke;
                        stiffK[cenx*dim+nei[cenx][s]] += -c*xi[0]*xi[0]*weight[cenx][s]/pow(xi.Norm2(),3)*broke;
                        stiffK[cenx*dim+ceny] += c*xi[0]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3)*broke;
                        stiffK[cenx*dim+nei[cenx][s]+dim/2] += -c*xi[0]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3)*broke;

                        stiffK[ceny*dim+ceny] += c*xi[1]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3)*broke;
                        stiffK[ceny*dim+nei[cenx][s]+dim/2] += -c*xi[1]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3)*broke;
                        stiffK[ceny*dim+cenx] += c*xi[0]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3)*broke;
                        stiffK[ceny*dim+nei[cenx][s]] += -c*xi[0]*xi[1]*weight[cenx][s]/pow(xi.Norm2(),3)*broke;
                    }
                    //cout<<endl;
                    rhs[cenx] = (2*u_old[cenx]-u_oldold[cenx])/dt/dt*den;
                    rhs[ceny] = (2*u_old[ceny]-u_oldold[ceny])/dt/dt*den;
                    stiffK[cenx*dim+cenx] += den/dt/dt;
                    stiffK[ceny*dim+ceny] += den/dt/dt; 
                }
                else{
                    rhs[cenx] = 0.0;
                    rhs[ceny] = -double(ID[cenx] == 2)*3.2*t;
                    stiffK[cenx*dim+cenx] = 1;
                    stiffK[ceny*dim+ceny] = 1; 
                }
                   
            }
        }
        cout<<"solve for displacement"<<endl;
        dgetrf_(&dim,&dim,stiffK,&dim,IPIV,&info);
        dgetrs_(&C,&dim,&nrhs,stiffK,&dim, IPIV,rhs,&dim,&info);
        
        diff = 0;
       



        for(int i=0;i<dim;i++){
            u[i] = rhs[i];
        }    
        //update broken bonds
        int Nbroke = 0;
        for(int i=0;i<2*N-1;i++){
            for(int j=0;j<N;j++){
                if(ID[i*N+j] == 0){
                    for(int s=0;s<nei[i*N+j].size();s++){
                        xi[0] = x[nei[i*N+j][s]] + u[nei[i*N+j][s]]-x[i*N+j] -u[i*N+j];
                        xi[1] = y[nei[i*N+j][s]] + u[nei[i*N+j][s]+dim/2]-y[i*N+j] - u[i*N+j+dim/2];
                        gamma[0] = x[nei[i*N+j][s]] - x[i*N+j];
                        gamma[1] = y[nei[i*N+j][s]] - y[i*N+j];
                        double stretch = xi.Norm2()/gamma.Norm2()-1.0;
                            if(!Bondbroke[i*N+j][s] && stretch > smax && ID[nei[i*N+j][s]] == 0){
                                Bondbroke[i*N+j][s] = true;
                                Nbroke += 1;
                            }
                    }
                }
            }
        }
    cout<<"bonds broken at step "<<step<<": "<<Nbroke<<endl;
    //update uold and uoldold 
    for(int i=0;i<dim;i++){
        u_oldold[i] = u_old[i];
        u_old[i] = u[i];
    }

    /*
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
    */
    if(step % 10 == 0){
    ufile.open("displacement_y"+std::to_string(step)+".csv");
    ufile1.open("displacement_x"+std::to_string(step)+".csv");
    ufile<<"x y v"<<endl;
    ufile1<<"x y u"<<endl;
        for(int i=0;i<2*N-1;i++){
            for(int j=0;j<N;j++){
                if(ID[i*N+j] == 0){
                cenx = i*N+j;
                ceny = i*N+j+(2*N-1)*N;
                ufile1<<x[cenx]<<" "<<y[cenx]<<" "<<u[cenx]<<" "<<endl;
                ufile<<x[cenx]<<" "<<y[cenx]<<" "<<u[ceny]<<" "<<endl;
                }
            }
        //ufile1<<endl;
        //ufile<<endl;
        }
        ufile.close();
        ufile1.close();
        for(int i=0;i<2*N-1;i++){
            for(int j=0;j<N;j++){
                int nneigh = 0;
                double sum = 0;
                if(ID[i*N+j] == 0){
                for(int s=0;s<nei[i*N+j].size();s++){
                    if(Bondbroke[i*N+j][s]) sum++;
                        nneigh++;
                        }
                        damage[i*N+j] = sum/double(nneigh);
                    }
                }
            }
            ofstream dfile;
            dfile.open("damage_"+std::to_string(step)+".csv");
            dfile<<"x y damage"<<endl;
            for(int i=0;i<2*N-1;i++){
                for(int j=0;j<N;j++){
                    if(ID[i*N+j] == 0){
                        dfile<<x[i*N+j]<<" "<<y[i*N+j]<<" "<<damage[i*N+j]<<endl;
                    }
                }
            }
            dfile.close();
    
    }
    }
    delete[] x;
    delete[] u;
    delete[] rhs;
    delete[] stiffK;
    delete[] IPIV;
    return 0;
}
