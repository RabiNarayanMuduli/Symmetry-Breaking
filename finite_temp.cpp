/* Libraries */
#include <execution>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <cstring>
using namespace std;
using execution::par; // Parallel execution policy

/* Constants */
const double pi = M_PI;
const double T = 0.100;
const double D = 1;
const double w2 = 0.25;
const double mass = 0.0;
const double eps = 1e-15;
const double epslog = -w2*log(1e-28);

const int Ny = 101;
const int Nm = 101;

typedef vector<double> VECTOR;
typedef vector<VECTOR> MATRIX;

// returns the weight and abscissa of legendre polynomial
void gauleg(const double x1, const double x2, VECTOR &x,VECTOR &w){
    const int N=x.size();  // number of points
    const double EPS=1e-14; // tolerance
    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for(i=0;i<N/2;i++){
        
        z=cos(pi*(i+0.75)/(N+0.5));
        while(abs(z-z1)>EPS){  // Newton's method of root finding
            p1 =1;
            p2 = 0;
            for(j=0;j<N;j++){
                p3 = p2;
                p2 = p1;
                p1 = ((2*j+1)*z*p2-j*p3)/(j+1);
            }
            pp = N*(z*p1-p2)/(z*z-1);
            z1 = z;
            z = z1-p1/pp;
        }
        x[i] = xm-xl*z;
        x[N-1-i] = xm+xl*z;
        w[i] = 2*xl/((1-z*z)*pp*pp);
        w[N-1-i] = w[i];
    }
}
 // Returns the Matsubura Frequencies
double omega(int n){
	n-=Nm/2;
	return pi*T*(2*n+1);
}
// Returns the frequency difference
double Omega(int n,int m){
	return omega(n)-omega(m);
}

// Returns the square of frequency difference
double Omega2(int n,int m){
	return pow(Omega(n,m),2);
}

// returns the distance between the momentum 
double distance(double t,double k,double l,int n,int m){
	return k*k+l*l-2*k*l*t+Omega2(n,m);
}

double exp_dist(double dist){
	return 4*pi*pi*D * exp(-dist/w2)*pow(w2,-3);
}

double no_t(double k,double l, int n, int m){
    double c = 2*k*l/w2;
    double dist1 = distance(1,k,l,n,m);
    double dist2 = distance(-1,k,l,n,m);
    return (exp_dist(dist1)-exp_dist(dist2))/c;
}

double one_t(double k,double l,int n,int m){
    double c = 2*k*l/w2;
    double dist1 = distance(1,k,l,n,m);
    double dist2 = distance(-1,k,l,n,m);
    return (exp_dist(dist1)+exp_dist(dist2)-no_t(k,l,n,m))/c;
}

double two_t(double k,double l,int n,int m){
    double c = 2*k*l/w2;
    double dist1 = distance(1,k,l,n,m);
    double dist2 = distance(-1,k,l,n,m);
    return (exp_dist(dist1)-exp_dist(dist2)-2*one_t(k,l,n,m))/c;
}

// The functions below returns the integration with respect to t in A,B,C
//
//
double funcA1(double k,double l,int n,int m){
	double a1 = distance(0,k,l,n,m);
	double b1 = a1*l/k* one_t(k,l,n,m)-2*l*l*two_t(k,l,n,m);
	double b2 = 2*k*l*one_t(k,l,n,m)-2*l*l*no_t(k,l,n,m)-2*l*l*two_t(k,l,n,m)+2*pow(l,3)/k*one_t(k,l,n,m);
	return b1+b2;
}

double funcA2( double k,double l,int n,int m){
	return 2*Omega(n,m)*omega(m)*(no_t(k,l,n,m)-l/k*one_t(k,l,n,m));  
}

double funcB( double k,double l,int n,int m){
	double a1 = distance(0,k,l,n,m);
	double b1 = a1*no_t(k,l,n,m);
	double b2 = -2*k*l*one_t(k,l,n,m);
	return b1+b2;
}

double funcC1( double k,double l,int n,int m){
	return Omega(n,m)*2*(k*l*one_t(k,l,n,m)-l*l*no_t(k,l,n,m));

}

double funcC2( double k,double l,int n,int m){
	double a1 = distance(0,k,l,n,m);
	double b1 = a1*no_t(k,l,n,m);
	double b2 = -2*k*l*one_t(k,l,n,m);
	double b3 = 2*Omega2(n,m)*no_t(k,l,n,m);
	return omega(m)*(b1+b2+b3);
}

void print(int no,VECTOR& y,MATRIX&A,MATRIX &B,MATRIX &C){
    string file_name = "output_temp_"+to_string(no)+".txt";
    ofstream file(file_name);
	file <<"# log(p) n A(p,n) B(p,n) C(p,n)\n";
	for(int ind=0 ;ind<y.size();ind++){
		for(int jj=0;jj<Nm;jj++)
			file <<2*log10(y[ind])<<" "<<jj-Nm/2<<" "<<A[ind][jj] <<" "<<B[ind][jj]<<" "<< C[ind][jj]<<endl; 
	}
}

int main(){
	cout <<"Starting main\n";
	VECTOR y(Ny),wy(Ny);
	gauleg(-2,2,y,wy);
	u_int16_t index[Ny];
	iota(index,index+Ny,0);
    for_each(par,index,index+Ny,[&](auto ind){
        y[ind] = exp10(y[ind]);
        wy[ind] *= y[ind]*log(10.0);
    });
    MATRIX A(Ny,VECTOR(Nm,1));
	MATRIX B(Ny,VECTOR(Nm,1));
	MATRIX C(Ny,VECTOR(Nm,1));
	MATRIX A1(Ny,VECTOR(Nm,1));
	MATRIX B1(Ny,VECTOR(Nm,1));
	MATRIX C1(Ny,VECTOR(Nm,1));
    double error = 0;
    auto sigma_A = [&](int j,int m){
		if (m<0 or m>=Nm) return (double)1/(y[j]*y[j]+pow(omega(m),2)+mass*mass);
		return A[j][m]/(pow(y[j]*A[j][m],2)+pow(omega(m)*C[j][m],2)+pow(B[j][m],2));
	};	
	auto sigma_B = [&](int j,int m){
		double ll;	
		if (m<0 or m>=Nm) {
			ll =  mass/(y[j]*y[j]+pow(omega(m),2)+mass*mass);
			// cout <<y[j]<<endl;	
		}
		else {
			ll =  B[j][m]/(pow(y[j]*A[j][m],2)+pow(omega(m)*C[j][m],2)+pow(B[j][m],2));
			// cout <<"Executed Usual\n";	
		}
		// if (ll !=0)cout << "j = "<<j<<", m = "<<m<<", B = "<<setprecision(15)<<ll<<endl;
		return ll; 
	};
	// sigma_C
	auto sigma_C = [&](int j,int m){
		if (m<0 or m>= Nm) return (double)1/(y[j]*y[j]+pow(omega(m),2)+mass*mass);
		double ll = C[j][m]/(pow(y[j]*A[j][m],2)+pow(omega(m)*C[j][m],2)+pow(B[j][m],2));
		// cout << "Inside = "<< m <<endl;
		return ll;
	};
    bool converged = false;
	int no_iter = 1;
	cout <<"Loop About to start\n";
	while (!converged ){
		converged=true;
		cout <<"Iteration no = "<<no_iter ;
		no_iter++;
        error = 0;
		for_each(par,index,index+Ny,[&](auto ind){
				double sumAn,sumBn,sumCn;
				/* double sumA,sumB,sumC; */
				for(int n = 0;n<Nm;n++){
					sumAn = sumBn = sumCn = 0;					
					int m = n;
					while (  Omega2(n, m)<epslog){
					// cout << m << endl;

						int j = ind;
						while(j>-1 and pow(y[j]-y[ind],2)+Omega2(n,m)<epslog){
							sumAn += wy[j]*y[j]*y[j]*(funcA1( y[ind], y[j], n, m)*sigma_A(  j, m)+funcA2( y[ind], y[j], n, m)*sigma_C(  j, m));
							
							double cbc = wy[j]*y[j]*y[j]*funcB(y[ind], y[j], n, m)*sigma_B(  j,m);
							sumBn += cbc;
							// cout << "sum Bn = "<<cbc <<endl ;
							sumCn += wy[j]*y[j]*y[j]*(funcC1( y[ind], y[j], n, m)*sigma_A( j, m)+funcC2( y[ind], y[j], n, m)*sigma_C( j, m));
							j--;
						}
						j = ind+1;
						while(j<Ny and pow(y[j]-y[ind],2)+Omega2(n, m)<epslog){
							sumAn += wy[j]*y[j]*y[j]*(funcA1( y[ind], y[j], n, m)*sigma_A(  j, m)+funcA2( y[ind], y[j], n, m)*sigma_C( j, m));
							sumBn += wy[j]*y[j]*y[j]*funcB( y[ind], y[j], n, m)*sigma_B( j, m);
							sumCn += wy[j]*y[j]*y[j]*(funcC1( y[ind], y[j], n, m)*sigma_A(j,m)+ funcC2( y[ind], y[j], n, m)*sigma_C( j, m));
							j++;
						};
						m--;

					}
					// exit(0);
					m = n+1;
					while (  Omega2(n, m)<epslog){
						int j = ind;

						while(j>-1 and pow(y[j]-y[ind],2)+Omega2(n, m)<epslog){

							sumAn += wy[j]*y[j]*y[j]*  (funcA1( y[ind], y[j], n, m)*sigma_A( j, m)+funcA2( y[ind], y[j], n, m)*sigma_C(  j, m));

							sumBn += wy[j]*y[j]*y[j]*  funcB( y[ind], y[j], n,m)*sigma_B( j, m);
							sumCn += wy[j]*y[j]*y[j]*  (funcC1( y[ind], y[j], n, m)*sigma_A( j, m)+funcC2( y[ind], y[j], n, m)*sigma_C(  j, m));
							j--;
						}

						j = ind+1;
						while(j<Ny and pow(y[j]-y[ind],2)+Omega2(n, m)<epslog){
							sumAn += wy[j]*y[j]*y[j]*  (funcA1( y[ind], y[j], n, m)*sigma_A( j, m)+funcA2( y[ind], y[j], n, m)*sigma_C( j, m));
							sumBn += wy[j]*y[j]*y[j]*  funcB( y[ind], y[j], n, m)*sigma_B( j, m);
							sumCn += wy[j]*y[j]*y[j]*  (funcC1( y[ind], y[j], n, m)*sigma_A( j, m)+funcC2( y[ind], y[j], n, m)*sigma_C( j, m));
							j++;
						};
						m++;
					}

					sumAn *= T/3/pi/pi;
					sumBn *= T/pi/pi;
					sumCn *= T/3/pi/pi/omega(n);
					error += abs(A[ind][n]-1-sumAn);
					error += abs(B[ind][n]-mass-sumBn);
					error += abs(C[ind][n]-1-sumCn);
					A1[ind][n] = 1+sumAn;
					B1[ind][n] = mass+sumBn;
					C1[ind][n] = 1+sumCn;
				}	
		});
		error = error /Ny/Nm;
        cout <<", Error = "<<setprecision(10)<<error << endl;
		std :: swap(A,A1);
		std :: swap(B,B1);
		std :: swap(C,C1);
		converged = (error<eps);
	}	
	cout << "converged = "<< converged << endl;
	print(0,y,A,B,C);
	A.clear();
	B.clear();
	C.clear();
	A1.clear();
	B1.clear();
	C1.clear();
	y.clear();
	wy.clear();
	return 0;
}