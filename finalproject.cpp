#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
const int Nx=129;
const int Ny=129;

int main(int argc, char** argv) {
	
	char filenameuc[128];
	FILE *tempfileuc;	
	sprintf(filenameuc, "Uc.dat");
	tempfileuc = fopen(filenameuc, "w+");

	char filenamemin[128];
	FILE *tempfilemin;	
	sprintf(filenamemin, "Min.dat");
	tempfilemin = fopen(filenamemin, "w+");
	
		
	int imin,jmin;
	double slmin,vormin;
	
	int i,j,k,l,m,n;
	int itime=0;
	double Lx,Ly,Lz;
	double e;
	double dt;
	double dx,dy,dz;
	double rhs;
	double convnx,convny,convnz;
	double convox,convoy,convoz;
	double diffx,diffy,diffz;
	double Hn,Ho,diffall;
	double ro;
	double Re;
	double eps1;
	double eps2;
	double ai,bi,ci,at,bt,ct;
	double invRe,invdx,invdx2,invdt;
	double invdy,invdy2,invdz,invdz2;
	double buf;
	double w;
	
	// U velocity
	double **Uold;
	double **Un;
	double **U1;
	double **Ud;
	double **Uh;
	double **Unew;
	
	// V velocity
	double **Vold;
	double **Vn;
	double **V1;
	double **Vd;
	double **Vh;
	double **Vnew;
	
	// Pseudo-pressure
	double **fi;
	double **fibuf;
	double **finew;
	
	// Pressure
	double **P;
	
	//Postproc
	double **SL;
	double **SL1;
	double **Vor;
		
	double ***Pcof;
	
	// Allocation
	
	//U velocity
	Uold=new double*[Nx+1];
	for(i=0;i<Nx+1;i++){
		Uold[i]=new double[Ny+2];
	}
	
	Un=new double*[Nx+1];
	for(i=0;i<Nx+1;i++){
		Un[i]=new double[Ny+2];
	}
	
	U1=new double*[Nx+1];
	for(i=0;i<Nx+1;i++){
		U1[i]=new double[Ny+2];
	}
	
	Ud=new double*[Nx+1];
	for(i=0;i<Nx+1;i++){
		Ud[i]=new double[Ny+2];
	}

	Uh=new double*[Nx+1];
	for(i=0;i<Nx+1;i++){
		Uh[i]=new double[Ny+2];
	}
	
	Unew=new double*[Nx+1];
	for(i=0;i<Nx+1;i++){
		Unew[i]=new double[Ny+2];
	}
	
	//V velocity
	Vold=new double*[Nx+2];
	for(i=0; i<Nx+2; i++){
		Vold[i]=new double[Ny+1];
	}

	Vn=new double*[Nx+2];
	for(i=0; i<Nx+2; i++){
		Vn[i]=new double[Ny+1];
	}

	V1=new double*[Nx+2];
	for(i=0; i<Nx+2; i++){
		V1[i]=new double[Ny+1];
	}
	
	Vd=new double*[Nx+2];
	for(i=0; i<Nx+2; i++){
		Vd[i]=new double[Ny+1];
	}

	Vh=new double*[Nx+2];
	for(i=0; i<Nx+2; i++){
		Vh[i]=new double[Ny+1];
	}
	
	Vnew=new double*[Nx+2];
	for(i=0; i<Nx+2; i++){
		Vnew[i]=new double[Ny+1];
	}
	
	SL=new double*[Nx+2];
	for(i=0;i<Nx+2;i++){
		SL[i]=new double[Ny+2];
	}	
	
	SL1=new double*[Nx+2];
	for(i=0;i<Nx+2;i++){
		SL1[i]=new double[Ny+2];
	}
	
	Vor=new double*[Nx+2];
	for(i=0;i<Nx+2;i++){
		Vor[i]=new double[Ny+2];
	}	
	
	// Pseudo-pressure
	fi=new double*[Nx+2];
	for(i=0;i<Nx+2;i++){
		fi[i]=new double[Ny+2];
	}
	
	finew=new double*[Nx+2];
	for(i=0;i<Nx+2;i++){
		finew[i]=new double[Ny+2];
	}		

	fibuf=new double*[Nx+2];
	for(i=0;i<Nx+2;i++){
		fibuf[i]=new double[Ny+2];
	}
	
	// Pressure
	P=new double*[Nx+2];
	for(i=0;i<Nx+2;i++){
		P[i]=new double[Ny+2];
	}	
	
	Pcof=new double**[Nx+2];
	for(i=0;i<Nx+2;i++){
		Pcof[i]=new double*[Ny+2];
	}
	for(i=0;i<Nx+2;i++){
		for(j=0;j<Ny+2;j++){
			Pcof[i][j]=new double[5]; 
		}
	}

	
	//Initialization
	Lx=1.0;
	Ly=1.0;
	
	dx=Lx/Nx;
	dy=Ly/Ny;	
	dt=0.5*dx;
	
	ro=1.0;
	Re=100.0;
	eps1=1.0;
	eps2=1.0;
	
	// Sor
	w=1.85;
	
	invRe=1.0/Re;
	invdx=1.0/dx;
	invdx2=1.0/dx/dx;
	invdy=1.0/dy;
	invdy2=1.0/dy/dy;	
	invdt=1.0/dt;
	
	cout<<"dx="<<dx<<endl;
	cout<<"dy="<<dy<<endl;
	cout<<"invdx="<<invdx<<endl;
	cout<<"invdy="<<invdy<<endl;
	cout<<"invdx2="<<invdx2<<endl;
	cout<<"invdy2="<<invdy2<<endl;
	cout<<"Re="<<Re<<endl;
	cout<<"invRe="<<invRe<<endl;
	cout<<"dt="<<dt<<endl;

	for(i=0;i<Nx+1;i++){
		for(j=0;j<Ny+2;j++){
			Uold[i][j]=0.0;
			Un[i][j]=0.0;
			U1[i][j]=0.0;
			Ud[i][j]=0.0;
			Uh[i][j]=0.0;
			Unew[i][j]=0.0;			
		}
	}

	for(i=0;i<Nx+2;i++){
		for(j=0;j<Ny+1;j++){
			Vold[i][j]=0.0;
			Vn[i][j]=0.0;
			V1[i][j]=0.0;
			Vd[i][j]=0.0;
			Vh[i][j]=0.0;
			Vnew[i][j]=0.0;			
		}
	}		
	
	for(i=0;i<Nx+2;i++){
		for(j=0;j<Ny+2;j++){
			SL[i][j]=0.0;
			SL1[i][j]=0.0;			
			Vor[i][j]=0.0;
			fi[i][j]=0.0;
			finew[i][j]=0.0;		
			fibuf[i][j]=0.0;
			P[i][j]=0.0;
		}
	}

	for(i=0;i<Nx+2;i++){
		for(j=0;j<Ny+2;j++){
			for(k=0;k<5;k++){
				Pcof[i][j][k]=0.0;
			}
		}
	}

	// TDMA
	
	double *d;
	double *cc,*dd,*x;
	double *a,*b,*c;
	a=new double[Nx+2];
	b=new double[Nx+2];
	c=new double[Nx+2];
	d=new double[Nx+2];
	cc=new double[Nx+2];
	dd=new double[Nx+2];
	x=new double[Nx+2];
	
	for(i=0;i<Nx+2;i++){
		a[i]=0.0;
		b[i]=0.0;
		c[i]=0.0;
		d[i]=0.0;
		cc[i]=0.0;
		dd[i]=0.0;
		x[i]=0.0;	
	}
	
	itime=0;
		
	for(i=1;i<Nx+1;i++){
		for(j=1;j<Ny+1;j++){
			
			Pcof[i][j][1]=invdx2;
			Pcof[i][j][2]=invdy2;
			Pcof[i][j][3]=invdx2;
			Pcof[i][j][4]=invdy2;			
			
			if(i==1){
				Pcof[i][j][3]=0.0;
			}
			if(j==1){
				Pcof[i][j][4]=0.0;
			}
			
		}
	}
	
	for(i=1;i<Nx+1;i++){
		for(j=1;j<Ny+1;j++){
			Pcof[i][j][0]=Pcof[i][j][1]+Pcof[i][j][2]+Pcof[i][j][3]+Pcof[i][j][4];
			Pcof[i][j][0]=1.0/Pcof[i][j][0];
		}
	}


	itime=0;
	///////////////////////////////////////////////
	// BC Old
	///////////////////////////////////////////////
	for(i=0;i<Nx+1;i++){
		Uold[i][0]=-Uold[i][1];
		Uold[i][Ny+1]=2.0-Uold[i][Ny];
	}
	for(j=0;j<Ny+2;j++){
		Uold[0][j]=0.0;
		Uold[Nx][j]=0.0;
	}
		
	for(j=0;j<Ny+1;j++){
		Vold[0][j]=-Vold[1][j];
		Vold[Nx+1][j]=-Vold[Nx][j];
	}
	for(i=0;i<Nx+2;i++){
		Vold[i][0]=0.0;
		Vold[i][Ny]=0.0;
	}
	
	

	while(eps1>0.0000000001 ){
		
		/////////////////////////////////////////////
		//Boundary condition
		/////////////////////////////////////////////
		for(i=0;i<Nx+1;i++){
			Un[i][0]=-Un[i][1];
			Un[i][Ny+1]=2.0-Un[i][Ny];
		}
		for(j=0;j<Ny+2;j++){
			Un[0][j]=0.0;
			Un[Nx][j]=0.0;
		}
		
		for(j=0;j<Ny+1;j++){
			Vn[0][j]=-Vn[1][j];
			Vn[Nx+1][j]=-Vn[Nx][j];
		}
		for(i=0;i<Nx+2;i++){
			Vn[i][0]=0.0;
			Vn[i][Ny]=0.0;
		}
		
		/////////////////////////////////////////////
		// U1
		/////////////////////////////////////////////
		for(i=1;i<Nx;i++){
			
			for(j=1;j<Ny+1;j++){
				d[j]=0.0;
				cc[j]=0.0;
				dd[j]=0.0;
				x[j]=0.0;
				a[j]=0.0;
				b[j]=0.0;
				c[j]=0.0;
			}
			
			ai=-dt*invRe*invdy2*0.5;
			bi=1+2.0*dt*invRe*invdy2*0.5;
			ci=ai;
						
			for(j=1;j<Ny+1;j++){
				a[j]=ai;
				b[j]=bi;
				c[j]=ci;
			}
			
			for(j=1;j<Ny+1;j++){
				convox=invdx*(0.5*(Uold[i+1][j]*Uold[i+1][j]+Uold[i][j]*Uold[i][j])-0.5*(Uold[i-1][j]*Uold[i-1][j]+Uold[i][j]*Uold[i][j]));
				convoy=invdy*(0.25*(Uold[i][j+1]+Uold[i][j])*(Vold[i][j]+Vold[i+1][j])-0.25*(Uold[i][j]+Uold[i][j-1])*(Vold[i][j-1]+Vold[i+1][j-1]));
				Ho=convox+convoy;
			
				convnx=invdx*(0.5*(Un[i+1][j]*Un[i+1][j]+Un[i][j]*Un[i][j])-0.5*(Un[i-1][j]*Un[i-1][j]+Un[i][j]*Un[i][j]));
				convny=invdy*(0.25*(Un[i][j+1]+Un[i][j])*(Vn[i][j]+Vn[i+1][j])-0.25*(Un[i][j]+Un[i][j-1])*(Vn[i][j-1]+Vn[i+1][j-1]));
				Hn=convnx+convny;

				diffx=invdx2*(Un[i+1][j]-2.0*Un[i][j]+Un[i-1][j]);
				diffy=invdy2*(Un[i][j+1]-2.0*Un[i][j]+Un[i][j-1]);
				diffall=diffx+diffy;
			
				d[j]=0.5*dt*(-3.0*Hn+Ho+2.0*invRe*diffall);
			}

			// Forward
			cc[1]=c[1]/b[1];
			for(j=2;j<Ny+1;j++){
				cc[j]=c[j]/(b[j]-a[j]*cc[j-1]);
			}
			
			dd[1]=d[1]/b[1];
			for(j=2;j<Ny+1;j++){
				dd[j]=(d[j]-a[j]*dd[j-1])/(b[j]-a[j]*cc[j-1]);
			}
			
			//Backward
			x[Ny]=dd[Ny];
			for(j=Ny-1;j>=1;j--){
			x[j]=dd[j]-cc[j]*x[j+1];
			}						
			
			for(j=1;j<Ny+1;j++){
				U1[i][j]=x[j];
			}
		
		}
		
		/////////////////////////////////////////////
		// Ud
		/////////////////////////////////////////////
		for(j=1;j<Ny+1;j++){
		
			for(i=1;i<Nx;i++){
				d[i]=0.0;
				cc[i]=0.0;
				dd[i]=0.0;
				x[i]=0.0;
				a[i]=0.0;
				b[i]=0.0;
				c[i]=0.0;
			}

			ai=-dt*invRe*invdx2*0.5;
			bi=1+2.0*dt*invRe*invdx2*0.5;
			ci=ai;
			
			for(i=1;i<Nx;i++){
				a[i]=ai;
				b[i]=bi;
				c[i]=ci;
			}
			
			for(i=1;i<Nx;i++){
				d[i]=U1[i][j];
			}
			
			// Forward
			cc[1]=c[1]/b[1];
			for(i=2;i<Nx;i++){
				cc[i]=c[i]/(b[i]-a[i]*cc[i-1]);
			}
			
			dd[1]=d[1]/b[1];
			for(i=2;i<Nx;i++){
				dd[i]=(d[i]-a[i]*dd[i-1])/(b[i]-a[i]*cc[i-1]);
			}
			
			//Backward
			x[Nx-1]=dd[Nx-1];
			for(i=Nx-2;i>=1;i--){
			x[i]=dd[i]-cc[i]*x[i+1];
			}						
			
			for(i=1;i<Nx;i++){
				Ud[i][j]=x[i];
			}
						
		}
		
		/////////////////////////////////////////////
		// V1
		/////////////////////////////////////////////		
		for(i=1;i<Nx+1;i++){
			
			for(j=1;j<Ny;j++){
				d[j]=0.0;
				cc[j]=0.0;
				dd[j]=0.0;
				x[j]=0.0;
				a[j]=0.0;
				b[j]=0.0;
				c[j]=0.0;
			}
			
			ai=-dt*invRe*invdy2*0.5;
			bi=1+2.0*dt*invRe*invdy2*0.5;
			ci=ai;
						
			for(j=1;j<Ny;j++){
				a[j]=ai;
				b[j]=bi;
				c[j]=ci;
			}
			
			for(j=1;j<Ny;j++){
				convox=invdx*(0.25*(Uold[i][j+1]+Uold[i][j])*(Vold[i+1][j]+Vold[i][j])-0.25*(Uold[i-1][j+1]+Uold[i-1][j])*(Vold[i-1][j]+Vold[i][j]));
				convoy=invdy*(0.5*(Vold[i][j+1]*Vold[i][j+1]+Vold[i][j]*Vold[i][j])-0.5*(Vold[i][j-1]*Vold[i][j-1]+Vold[i][j]*Vold[i][j]));
				Ho=convox+convoy;

				convnx=invdx*(0.25*(Un[i][j+1]+Un[i][j])*(Vn[i+1][j]+Vn[i][j])-0.25*(Un[i-1][j+1]+Un[i-1][j])*(Vn[i-1][j]+Vn[i][j]));
				convny=invdy*(0.5*(Vn[i][j+1]*Vn[i][j+1]+Vn[i][j]*Vn[i][j])-0.5*(Vn[i][j-1]*Vn[i][j-1]+Vn[i][j]*Vn[i][j]));
				Hn=convnx+convny;

				diffx=invdx2*(Vn[i+1][j]-2.0*Vn[i][j]+Vn[i-1][j]);
				diffy=invdy2*(Vn[i][j+1]-2.0*Vn[i][j]+Vn[i][j-1]);
				diffall=diffx+diffy;
			
				d[j]=0.5*dt*(-3.0*Hn+Ho+2.0*invRe*diffall);
			}

			// Forward
			cc[1]=c[1]/b[1];
			for(j=2;j<Ny;j++){
				cc[j]=c[j]/(b[j]-a[j]*cc[j-1]);
			}
			
			dd[1]=d[1]/b[1];
			for(j=2;j<Ny;j++){
				dd[j]=(d[j]-a[j]*dd[j-1])/(b[j]-a[j]*cc[j-1]);
			}
			
			//Backward
			x[Ny-1]=dd[Ny-1];
			for(j=Ny-2;j>=1;j--){
			x[j]=dd[j]-cc[j]*x[j+1];
			}						
			
			for(j=1;j<Ny;j++){
				V1[i][j]=x[j];
			}
		
		}		
		/////////////////////////////////////////////
		// Vd
		/////////////////////////////////////////////
		for(j=1;j<Ny;j++){
		
			for(i=1;i<Nx+1;i++){
				d[i]=0.0;
				cc[i]=0.0;
				dd[i]=0.0;
				x[i]=0.0;
				a[i]=0.0;
				b[i]=0.0;
				c[i]=0.0;
			}

			ai=-dt*invRe*invdx2*0.5;
			bi=1+2.0*dt*invRe*invdx2*0.5;
			ci=ai;
			
			for(i=1;i<Nx+1;i++){
				a[i]=ai;
				b[i]=bi;
				c[i]=ci;
			}
			
			for(i=1;i<Nx+1;i++){
				d[i]=V1[i][j];
			}
			
			// Forward
			cc[1]=c[1]/b[1];
			for(i=2;i<Nx+1;i++){
				cc[i]=c[i]/(b[i]-a[i]*cc[i-1]);
			}
			
			dd[1]=d[1]/b[1];
			for(i=2;i<Nx+1;i++){
				dd[i]=(d[i]-a[i]*dd[i-1])/(b[i]-a[i]*cc[i-1]);
			}
			
			//Backward
			x[Nx]=dd[Nx];
			for(i=Nx-1;i>=1;i--){
			x[i]=dd[i]-cc[i]*x[i+1];
			}						
			
			for(i=1;i<Nx+1;i++){
				Vd[i][j]=x[i];
			}
						
		}	
		
		/////////////////////////////////////////////
		// Uh && Vh
		/////////////////////////////////////////////
		for(i=0;i<Nx+1;i++){
			for(j=0;j<Ny+2;j++){
				Uh[i][j]=Un[i][j]+Ud[i][j];
			}
		}
		for(i=0;i<Nx+2;i++){
			for(j=0;j<Ny+1;j++){
				Vh[i][j]=Vn[i][j]+Vd[i][j];
			}
		}

		/////////////////////////////////////////////
		//Boundary condition
		/////////////////////////////////////////////
		for(i=0;i<Nx+1;i++){
			Uh[i][0]=-Uh[i][1];
			Uh[i][Ny+1]=2.0-Uh[i][Ny];
		}
		for(j=0;j<Ny+2;j++){
			Uh[0][j]=0.0;
			Uh[Nx][j]=0.0;
		}
		
		for(j=0;j<Ny+1;j++){
			Vh[0][j]=-Vh[1][j];
			Vh[Nx+1][j]=-Vh[Nx][j];
		}
		for(i=0;i<Nx+2;i++){
			Vh[i][0]=0.0;
			Vh[i][Ny]=0.0;
		}

		/////////////////////////////////////////////
		// Pseudo-pressure 1
		/////////////////////////////////////////////
		for(i=0;i<Nx+2;i++){
			for(j=0;j<Ny+2;j++){
				fi[i][j]=0.0;
			}
		}
		
		while(eps2>0.00000001){
			
			// BC
			for(i=0;i<Nx+2;i++){
				fi[i][0]=fi[i][1];
				fi[i][Ny+1]=fi[i][Ny];
			}
			for(j=0;j<Ny+2;j++){
				fi[0][j]=fi[1][j];
				fi[Nx+1][j]=fi[Nx][j];	
			}		
		
			for(i=1;i<Nx+1;i++){
				for(j=1;j<Ny+1;j++){
					buf =invdt*(Uh[i][j]-Uh[i-1][j])*invdx;
					buf+=invdt*(Vh[i][j]-Vh[i][j-1])*invdy;
					rhs=fi[i+1][j]*Pcof[i][j][1];
					rhs+=fi[i][j+1]*Pcof[i][j][2];
					rhs+=finew[i-1][j]*Pcof[i][j][3];
					rhs+=finew[i][j-1]*Pcof[i][j][4];
					finew[i][j]=(1.0-w)*fi[i][j]+w*Pcof[i][j][0]*(rhs-buf);
				}
			}
			
			double s=0.0;
			double s1=0.0;
			for(i=1;i<Nx+1;i++){
				for(j=1;j<Ny+1;j++){
					s=s+fabs(finew[i][j]-fi[i][j]);
					s1+=fabs(fi[i][j]);
				}
			}
			eps2=s/s1/(Nx)/(Ny);
		
			for(i=1;i<Nx+1;i++){
				for(j=1;j<Ny+1;j++){
					fi[i][j]=finew[i][j];
				}
			}

		}

		eps2=1.0;
		
		/////////////////////////////////////////
		// New
		/////////////////////////////////////////
		for(i=1;i<Nx;i++){
			for(j=1;j<Ny+1;j++){
				Unew[i][j]=Uh[i][j]-dt*invdx*(finew[i+1][j]-finew[i][j]);
			}
		}
		for(i=1;i<Nx+1;i++){
			for(j=1;j<Ny;j++){
				Vnew[i][j]=Vh[i][j]-dt*invdy*(finew[i][j+1]-finew[i][j]);
			}
		}
		//////////////////////////////////////////
		// Pressure
		//////////////////////////////////////////
		for(i=1;i<Nx+1;i++){
			for(j=1;j<Ny+1;j++){
				P[i][j]=finew[i][j]-0.5*dt*invRe*((finew[i+1][j]-2.0*finew[i][j]+finew[i-1][j])*invdx2+(finew[i][j+1]-2.0*finew[i][j]+finew[i][j-1])*invdy2);
			}
		}
		
		//////////////////////////////////////////
		// Convergence
		//////////////////////////////////////////	
		double s=0.0;
		double s1=0.0;
		
		for(i=1;i<Nx;i++){
			for(j=1;j<Ny+1;j++){
				s=s+fabs(Unew[i][j]-Un[i][j]);
				s1=s1+fabs(Un[i][j]);
			}
		}
		
		for(i=1;i<Nx+1;i++){
			for(j=1;j<Ny;j++){
				s=s+fabs(Vnew[i][j]-Vn[i][j]);
				s1=s1+fabs(Vn[i][j]);				
			}
		}	
		
		//eps1=s/s1/(2*Ny*Nx);
		eps1=s/(2*Ny*Nx);
		
		/////////////////////////////////////////////
		// Update
		/////////////////////////////////////////////
		for(i=0;i<Nx+1;i++){
			for(j=0;j<Ny+2;j++){
				Uold[i][j]=Un[i][j];
			}
		}
		
		for(i=0;i<Nx+2;i++){
			for(j=0;j<Ny+1;j++){
				Vold[i][j]=Vn[i][j];
			}
		}	
		
		for(i=0;i<Nx+1;i++){
			for(j=0;j<Ny+2;j++){
				Un[i][j]=Unew[i][j];
			}
		}
		
		for(i=0;i<Nx+2;i++){
			for(j=0;j<Ny+1;j++){
				Vn[i][j]=Vnew[i][j];
			}
		}	
		cout<<"EPS="<<eps1<<" Itime="<<itime<<endl;
		
		/////////////////////////////////////////////
		// Post process
		/////////////////////////////////////////////
		if(itime%100==3){
		char filename[128];
		FILE *tempfile;
		
		sprintf(filename, "N_ALL%06i.dat", itime);
		tempfile = fopen(filename, "w+");
		fprintf(tempfile, "TITLE = \"Example: 3D Plot\"\n");
		
		fprintf(tempfile, "VARIABLES = \"X\", \"Y\",\"U\",\"V\",\"fi\",\"P\",\"Vort\",\"Stream\"\n");
		fprintf(tempfile, "ZONE I=%i, J=%i, F=POINT\n",Nx,Ny);
		
		for(i=1;i<Nx+1;i++){
			for(j=1;j<Ny+1;j++){
				fprintf(tempfile, "%5.5f\t%5.5f\t%5.5f\t%5.5f\t%5.5f\t%5.5f\t%5.5f\t%5.5f\n", i*dx-0.5*dx, j*dy-0.5*dy, (Unew[i][j]+Unew[i-1][j])*0.5, (Vnew[i][j]+Vnew[i][j-1])*0.5, finew[i][j], P[i][j], -Vor[i][j], SL[i][j]);		
			}
		}
		fclose(tempfile);	
		}
		
		if(itime%100==1){
		char filename[128];
		FILE *tempfile;
		
		sprintf(filename, "U%06i.dat", itime);
		tempfile = fopen(filename, "w+");

		i=65;
			for(j=1;j<Ny+1;j++){
				fprintf(tempfile, "%5.5f\t%5.5f\n", j*dy-0.5*dy, (Unew[i][j]+Unew[i-1][j])*0.5);		
			}
		
		fclose(tempfile);	
		}		
		
		if(itime%100==2){
		char filename[128];
		FILE *tempfile;
		
		sprintf(filename, "V%06i.dat", itime);
		tempfile = fopen(filename, "w+");

		j=65;
		for(i=1;i<Nx+1;i++){
			//for(j=1;j<Ny+1;j++){
				fprintf(tempfile, "%5.5f\t%5.5f\n", i*dx-0.5*dx, (Vnew[i][j]+Vnew[i][j-1])*0.5);		
			//}
		}
		fclose(tempfile);	
		}
		
		//////////////////////////////////////////////////
		// Table1 && Table2 Ghia
		//////////////////////////////////////////////////
		if(itime%100==0){
			char filename[128];
			FILE *tempfile;
			
			sprintf(filename, "Table2_%06i.dat", itime);
			tempfile = fopen(filename, "w+");
	
			j=65;

			fprintf(tempfile, "%3i\t%5.5f\n", 129, (Vnew[129][j]+Vnew[129][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 125, (Vnew[125][j]+Vnew[125][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 124, (Vnew[124][j]+Vnew[124][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 123, (Vnew[123][j]+Vnew[123][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 122, (Vnew[122][j]+Vnew[122][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 117, (Vnew[117][j]+Vnew[117][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 111, (Vnew[111][j]+Vnew[111][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 104, (Vnew[104][j]+Vnew[104][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 65, (Vnew[65][j]+Vnew[65][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 31, (Vnew[31][j]+Vnew[31][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 30, (Vnew[30][j]+Vnew[30][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 21, (Vnew[21][j]+Vnew[21][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 13, (Vnew[13][j]+Vnew[13][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 11, (Vnew[11][j]+Vnew[11][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 10, (Vnew[10][j]+Vnew[10][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 9, (Vnew[9][j]+Vnew[9][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 1, (Vnew[1][j]+Vnew[1][j-1])*0.5);		

			fclose(tempfile);			
			
		}

		if(itime%100==0){
			char filename[128];
			FILE *tempfile;
			
			sprintf(filename, "Table1_%06i.dat", itime);
			tempfile = fopen(filename, "w+");
	
			j=65;

			fprintf(tempfile, "%3i\t%5.5f\n", 129, (Vnew[129][j]+Vnew[129][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 126, (Vnew[126][j]+Vnew[126][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 125, (Vnew[125][j]+Vnew[125][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 124, (Vnew[124][j]+Vnew[124][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 123, (Vnew[123][j]+Vnew[123][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 110, (Vnew[110][j]+Vnew[110][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 95, (Vnew[95][j]+Vnew[95][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 80, (Vnew[80][j]+Vnew[80][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 65, (Vnew[65][j]+Vnew[65][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 59, (Vnew[59][j]+Vnew[59][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 37, (Vnew[37][j]+Vnew[37][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 23, (Vnew[23][j]+Vnew[23][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 14, (Vnew[14][j]+Vnew[14][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 10, (Vnew[10][j]+Vnew[10][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 9, (Vnew[9][j]+Vnew[9][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 8, (Vnew[8][j]+Vnew[8][j-1])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 1, (Vnew[1][j]+Vnew[1][j-1])*0.5);		

			fclose(tempfile);			
			
		}
	
		if(itime%100==0){
			char filename[128];
			FILE *tempfile;
			
			sprintf(filename, "Table_1_%06i.dat", itime);
			tempfile = fopen(filename, "w+");
	
			i=65;

			fprintf(tempfile, "%3i\t%5.5f\n", 129, (Unew[i][129]+Unew[i-1][129])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 126, (Unew[i][126]+Unew[i-1][126])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 125, (Unew[i][125]+Unew[i-1][125])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 124, (Unew[i][124]+Unew[i-1][124])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 123, (Unew[i][123]+Unew[i-1][123])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n", 110, (Unew[i][110]+Unew[i-1][110])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n",  95, (Unew[i][95]+Unew[i-1][95])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n",  80, (Unew[i][80]+Unew[i-1][80])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n",  65, (Unew[i][65]+Unew[i-1][65])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n",  59, (Unew[i][59]+Unew[i-1][59])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n",  37, (Unew[i][37]+Unew[i-1][37])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n",  23, (Unew[i][23]+Unew[i-1][23])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n",  14, (Unew[i][14]+Unew[i-1][14])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n",  10, (Unew[i][10]+Unew[i-1][10])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n",   9, (Unew[i][9]+Unew[i-1][9])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n",   8, (Unew[i][8]+Unew[i-1][8])*0.5);		
			fprintf(tempfile, "%3i\t%5.5f\n",   1, (Unew[i][1]+Unew[i-1][1])*0.5);		

			fclose(tempfile);		
		}
	
	
		if(itime%100==0){
			/////////////////////////////////////////////////
			//Streamline function 
			/////////////////////////////////////////////////

			for(i=1;i<Nx+1;i++){
				for(j=1;j<Ny+1;j++){
					
					double s=0.0;
					
					for(k=1;k<=j;k++){
					s=s+(Unew[i][k]+Unew[i-1][k])*0.5;
					}
					
					SL1[i][j]=s*dy;				
				}
			}
		
			for(i=1;i<Nx+1;i++){
				SL[i][1]=SL1[i][1]/2.0;
				for(j=2;j<Ny+1;j++){
					SL[i][j]=(SL1[i][j]+SL1[i][j-1])/2.0;
				}
			} 

			char filename[128];
			FILE *tempfile;			
			sprintf(filename, "SL%06i.dat", itime);
			tempfile = fopen(filename, "w+");
			fprintf(tempfile, "TITLE = \"Example: 2D Plot\"\n");
			
			//fprintf(tempfile, "VARIABLES = \"X\", \"Y\",\"Z\",\"U\",\"V\",\"W\",\"T\",\"C1\",\"C2\",\"C3\"\n");
			fprintf(tempfile, "VARIABLES = \"X\", \"Y\",\"Streamline\"\n");
			fprintf(tempfile, "ZONE I=%i, J=%i, F=POINT\n",Nx,Ny);
			
			for(i=1;i<Nx+1;i++){
				for(j=1;j<Ny+1;j++){
					fprintf(tempfile, "%5.5f\t%5.5f\t%5.5f\n",i*dx-dx*0.5,j*dy-dy*0.5, SL[i][j]);
				}
			}
			fclose(tempfile);	
			
		
		}
		
		//////////////////////////////////////////////////////////
		// Vorticity
		///////////////////////////////////////////////////////////
		if(itime%100==1){		
			
			for(i=1;i<Nx+1;i++){
				for(j=1;j<Ny+1;j++){
					Vor[i][j]= 0.5*invdx*( 0.5*(Vnew[i+1][j]+Vnew[i+1][j-1]) - 0.5*(Vnew[i-1][j]+Vnew[i-1][j-1]))-0.5*invdy*(0.5*(Unew[i][j+1]+Unew[i-1][j+1])-0.5*(Unew[i][j-1]+Unew[i-1][j-1]));
				}
			}
			char filename[128];
			FILE *tempfile;			
			sprintf(filename, "Vor%06i.dat", itime);
			tempfile = fopen(filename, "w+");
			fprintf(tempfile, "TITLE = \"Example: 2D Plot\"\n");
			
			//fprintf(tempfile, "VARIABLES = \"X\", \"Y\",\"Z\",\"U\",\"V\",\"W\",\"T\",\"C1\",\"C2\",\"C3\"\n");
			fprintf(tempfile, "VARIABLES = \"X\", \"Y\",\"Vorticity\"\n");
			fprintf(tempfile, "ZONE I=%i, J=%i, F=POINT\n",Nx,Ny);
			
			for(i=1;i<Nx+1;i++){
				for(j=1;j<Ny+1;j++){
					fprintf(tempfile, "%5.5f\t%5.5f\t%5.5f\n",i*dx-dx*0.5,j*dy-dy*0.5, -Vor[i][j]);
				}
			}
			fclose(tempfile);						
		}
		
		
		if(itime%100==0){
			slmin=SL[1][1];
			imin=1;
			jmin=1;	
			
			for(i=1;i<Nx+1;i++){
				for(j=1;j<Ny+1;j++){
					if(SL[i][j]<slmin){
						slmin=SL[i][j];
						imin=i;
						jmin=j;
					}
				}
			}
			fprintf(tempfilemin, "%5.5f\t%5.5f\t%5.5f\t%5.5f\n", imin*dx-0.5*dx, jmin*dy-0.5*dy, slmin, -Vor[imin][jmin]);
			//cout<<imin*dx-05*dx<<"\t"<<jmin*dy-0.5*dy<<"\t"<<slmin<<endl;
		}
		
		
		
		if(itime%100==0){
			fprintf(tempfileuc, "%5.5f\t%5.5f\n",itime*dt, Unew[65][65]);			
		}
		itime++;
			
	}

	fclose(tempfileuc);
	fclose(tempfilemin);
	
	
		
		
}
