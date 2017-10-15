#include<iostream>
#include<fstream>
#include<cmath>
#include"Vector.h"
#include"Random64.h"

using namespace std;
/*------------------------------------------------DIMENSIONAL VARIABLES--------------------------------------------------------*/
 
double L=1,m=1,dt=1e-3;                                  /*Independent Variables:Length,Mass,Time*/
double g=2*L*pow(dt,-2);                                 /*Buckingham Group*/
double v=L/dt;                                           /*Buckingham Group*/
double I=m*L*L;                                          /*Buckingham Group*/
double Gamma=pow(L*dt*dt,-0.5);                          /*Buckingham Group*/
double KHertz=m*pow(dt,-2);                              /*Buckingham Group*/
double KSpring=1e-2*m*pow(L*dt*dt*dt*dt,-0.5);           /*Buckingham Group*/
double Kx=KSpring,Ky=KSpring,Kz=KSpring,Kcundall=KSpring;/*Constant of the springs on each axis and constant of the contact force*/
double Lx=100*L,Ly=100*L,Lz=100*L;                       /*Lengths for the animation*/
double GammaSpring=1e-1/dt,GammaTangential=1e-3/dt;      /*Coefficients of friction for the spring and the contact force*/
/*------------------------------------------------ADIMENSIONAL VARIABLES-------------------------------------------------------*/
const double MU=0.4;              /*Coefficient for contact force*/
const int Nx=60,Ny=3,Nz=1;        /*Number of springs in each axis*/
const int N2D=Nx*Ny ,N3D=Nx*Ny*Nz;/*Total number of springs*/
/*------------------------------------------------OMELYAN PEFRL CONSTANTS-------------------------------------------------------*/
const double Xi=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Chi=-0.06626458266981849;
/*------------------------------------------------CLASS DECLARATION------------------------------------------------------------*/
class Body;
class Collider;
class Oscillator;
/*------------------------------------------------CLASS BODY------------------------------------------------------------------*/
class Body{
private:
  vector3D R,V,F,RBreak;/*Translation Variables*/
  vector3D Theta,W,Tau; /*Rotational Variables*/
  double m,r,I;         /*Intrinsic Variables*/
public:
  void Start(double Rx0,double Ry0,double Rz0,double Vx0,double Vy0,double Vz0,
	     double m0,double r0,double x0Break,double y0Break,double z0Break,
	     double Theta0,double Psi0,double Phi0,double Wx0,double Wy0,double Wz0);
  /*Move Translation Variables*/
  void Move_R(double dt, double Constante);
  void Move_V(double dt, double Constante);
  /*Add/Erases Forces*/
  void EraseForce(void);
  void AddForce_x(double F0x);
  void AddForce_y(double F0y);
  void AddForce_z(double F0z);
  void AddForce(vector3D F0);
  /*Move Rotational Variables*/
  void Move_Theta(double dt, double Constante);
  void Move_W(double dt, double Constante);
  /*Add/Erases Torsions*/
  void EraseTorsion(void);
  void AddTorsion(vector3D Tau0);
  /*Draw the Trajectory of the Particle*/
  void Draw(void);
  void Draw3D(void);
  /*Inline Functions*/
  double X(void){return RBreak.x()+R.x();};
  double Y(void){return RBreak.y()+R.y();};
  double Z(void){return RBreak.z()+R.z();};
  double NormF(void){return norma(F);};
  double NormTau(void){return norma(Tau);};
  friend class Collider;
  friend class Oscillator;
};
/*-------------------------------------------------BODY'S FUNCTIONS-------------------------------------------------------------------*/
void Body::Start(double Rx0,double Ry0,double Rz0,double Vx0,double Vy0,double Vz0,
		 double m0,double r0,double x0Break,double y0Break,double z0Break,
		 double Theta0,double Psi0,double Phi0,double Wx0,double Wy0,double Wz0){
  R.cargue(Rx0,Ry0,Rz0);V.cargue(Vx0,Vy0,Vz0);m=m0;r=r0;RBreak.cargue(x0Break,y0Break,z0Break);
  Theta.cargue(Theta0,Psi0,Phi0);W.cargue(Wx0,Wy0,Wz0);I=2.0/5*m*r*r;}
void Body::Move_R(double dt, double Constante){R+=V*(Constante*dt);}
void Body::Move_V(double dt, double Constante){V+=F*(Constante*dt/m);}
void Body::EraseForce(void){F.cargue(0,0,0);}
void Body::AddForce_x(double F0x){vector3D F0; F0.cargue(F0x,0,0);F+=F0;}
void Body::AddForce_y(double F0y){vector3D F0; F0.cargue(0,F0y,0);F+=F0;}
void Body::AddForce_z(double F0z){vector3D F0; F0.cargue(0,0,F0z);F+=F0;}
void Body::AddForce(vector3D F0){F+=F0;}
void Body::Move_Theta(double dt, double Constante){Theta+=W*(Constante*dt);}
void Body::Move_W(double dt, double Constante){W+=Tau*(Constante*dt/I);}
void Body::EraseTorsion(void){Tau.cargue(0,0,0);}
void Body::AddTorsion(vector3D Tau0){Tau+=Tau0;}
void Body::Draw(void){
  cout<<", "<<X()<<"+"<<r<<"*cos(t),"<<Y()<<"+"<<r<<"*sin(t)";
  cout<<", "<<X()<<"+"<<r*cos(Theta.z())/(2*M_PI)<<"*t,"<<Y()<<"+"<<r*sin(Theta.z())/(2*M_PI)<<"*t";}
void Body::Draw3D(void){
  cout<<", "<<X()<<"+"<<r<<"sin(v)*cos(u),"<<Y()<<"+"<<r<<"*sin(v)*sin(u),"<<Z()<<"+"<<r<<"*cos(v)";
  cout<<", "<<X()<<"+"<<r*sin(Theta.y())*cos(Theta.z())/(2*M_PI*M_PI)<<"*u*v,"
      <<Y()<<"+"<<r*sin(Theta.y())*cos(Theta.z())/(2*M_PI*M_PI)<<"*u*v";}
/*------------------------------------------------CLASS COLLIDER--------------------------------------------------------------*/
class Collider{
private:
  vector3D l[N2D+1][N2D+1]; bool InCollision[N2D+1][N2D+1];
public:
  void Start(void);
  void AllForces(Body* Spring,double dt);
  void InteractionForceSpring(Body & Spring1, Body & Spring2);
  void InteractionForceSphere(Body & Spring1, Body & Spring2,vector3D & l,bool & InCollision,double dt);
};
/*--------------------------------------COLLIDER'S FUNCTIONS------------------------------------------------------------------------*/
void Collider::Start(void){int i,j; for(i=0;i<N2D+1;i++){for(j=i+1;j<N2D+1;j++){l[i][j].cargue(0,0,0); InCollision[i][j]=false;}}}
void Collider::AllForces(Body* Spring,double dt){
  int i,j;
  /*Erase all Forces and  Torsions*/
  for(i=0;i<N2D+1;i++){Spring[i].EraseForce();Spring[i].EraseTorsion();}
  /*Add Gravitational Force due by the Earth's Gravitational Field*/
  for(i=N2D;i<N2D+1;i++){Spring[i].AddForce_y(-Spring[i].m*g);}
  /*Add Viscous Force on the springs*/
  for(i=0;i<N2D;i++){Spring[i].AddForce(Spring[i].V*(-GammaSpring*Spring[i].m));}
  /*Add forces on y for extremum springs in y axis*/
  for(i=0;i<Nx;i++){Spring[i].AddForce_y(-Ky*Spring[i].R.y()); Spring[Nx*(Ny-1)+i].AddForce_y(-Ky*Spring[Nx*(Ny-1)+i].R.y());}
  /*Local Interaction with periodic boundary conditions*/
  for(j=0;j<Ny;j++){for(i=0;i<Nx;i++){
      if(j==Ny-1){Spring[i+j*Nx].AddForce_y(0);InteractionForceSpring(Spring[j*Nx+(i%Nx)],Spring[j*Nx+((i+1)%Nx)]);}/*Last x-line*/
      else{InteractionForceSpring(Spring[j*Nx+(i%Nx)],Spring[j*Nx+((i+1)%Nx)]);InteractionForceSpring(Spring[i],Spring[i+Nx]);}}}
  /*Sphere and Spring*/
  for(i=0;i<N2D;i++){for(j=i+1;j<N2D+1;j++){InteractionForceSphere(Spring[i],Spring[j],l[i][j],InCollision[i][j],dt);}}
}

void Collider::InteractionForceSpring(Body & Spring1, Body & Spring2){
  double drx=Spring2.R.x()-Spring1.R.x();
  double dry=Spring2.R.y()-Spring1.R.y();
  double drz=Spring2.R.z()-Spring1.R.z();
  double F1x,F1y,F1z; vector3D F1;
  F1.cargue(Kx*drx,Ky*dry,Kz*drz);
  Spring1.AddForce(F1); Spring2.AddForce(F1*(-1));
}
void Collider::InteractionForceSphere(Body & Spring1, Body & Spring2,vector3D & l,bool & InCollision,double dt){
  vector3D R21,n,V1,V2,Vc,Vcn,Vct,t,Fn,Ft,F2; 
  double R1x,R1y,R1z,R2x,R2y,R2z,r1,r2,d21,s,m1,m2,m12,componentVcn,componentFn,normVct,Ftmax,normFt;
  double ERFF=1e-8;
  R1x=Spring1.X(); R2x=Spring2.X(); R1y=Spring1.Y(); R2y=Spring2.Y(); R1z=Spring1.Z(); R2z=Spring2.Z(); 
  R21.cargue(R2x-R1x,R2y-R1y,R2z-R1z); d21=norma(R21);  s=(Spring1.r+Spring2.r)-d21;
  V1=Spring1.V; V2=Spring2.V; 
  if(s>0){ //If crashing,
    /*Geometry and Contact Dinamycs*/
    m1=Spring1.m;   m2=Spring2.m;   m12=(m1*m2)/(m1+m2);
    r1=Spring1.r;   r2=Spring2.r;
    n=R21/d21;
    /*Contact Velocity and Tangent Vector*/
    Vc=(V2-V1)-(Spring2.W^n)*r2-(Spring1.W^n)*r1;
    componentVcn=Vc*n; Vcn=n*componentVcn; Vct=Vc-Vcn;  normVct=norma(Vct);
    if(normVct<ERFF){t.cargue(0,0,0);} else{t=Vct/normVct;}
    
    /*Normal Forces*/
    /*Hertz's Force*/
    componentFn=KHertz*pow(s,1.5); 
    /*Plastic Disipation*/
    componentFn-=m12*sqrt(s)*Gamma*componentVcn; if(componentFn<0){componentFn=0;}
    Fn=n*componentFn;
    
    /*Tangential Forces*/
    /*Static Force*/
    l+=(Vct*dt);
    Ft=l*(-Kcundall);
    /*Kinetic Force*/
    Ftmax=MU*componentFn; normFt=norma(Ft);
    if(normFt>Ftmax){Ft=l*(-Ftmax/norma(l));}
    /*Viscous Force*/
    Ft-=m2*GammaTangential*Vct;
    
    /*Building Total Force*/
    F2=Fn+Ft;
    Spring2.AddForce(F2);      Spring2.AddTorsion((n*(-r2))^Ft);      
    Spring1.AddForce(F2*(-1)); Spring1.AddTorsion((n*r1)^(Ft*(-1))); 
    
    InCollision=true;
  }
  else if(InCollision==true){l.cargue(0,0,0); InCollision=false;}
}
/*------------------------------------------------CLASS RIGID BODY------------------------------------------------------------*/

/*------------------------------------------------CLASS OSCILLATOR------------------------------------------------------------*/
class Oscillator{
private:
  
public:
  void DrawInteractionSpring(Body & Spring1,Body & Spring2);
};
/*-------------------------------------OSCILLATOR'S FUNCTIONS--------------------------------------------*/
void Oscillator::DrawInteractionSpring(Body & Spring1,Body & Spring2){
  vector3D D,X; double d,a;
  D.cargue(Spring2.X()-Spring1.X(),Spring2.Y()-Spring1.Y(),Spring2.Z()-Spring1.Z()); d=norma(D);
  X.cargue(1,0,0);
  a=(D*X)/d; /*Angle between D and x-axis*/
  /*Cicloid Parametrization
  x=r*(t-sin(t))
  y=r*(1-cos(t))
  Active Rotation Matrix [cos -sin; sin cos]*/
  cout<<", "<<Spring1.X()<<"+"<<d<<"*(t-sin(2*t))*"<<cos(a)<<"-"<<d<<"*(1-cos(2*t))*"<<sin(a)
      <<"," <<Spring1.Y()<<"+"<<d<<"*(t-sin(2*t))*"<<sin(a)<<"+"<<d<<"*(1-cos(2*t))*"<<cos(a);
}
/*-------------------------------------GLOBAL'S FUNCTIONS------------------------------------------------*/
void StartAnimation(double x0Break, double y0Break, double h){
  //cout<<"set terminal gif animate"<<endl; 
  //cout<<"set output 'MySpring2D.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange ["<<-x0Break<<":"<<x0Break*Nx<<"]"<<endl;
  //cout<<"set yrange ["<<-y0Break<<":"<<y0Break*Ny+h<<"]"<<endl;
  //cout<<"set autoscale"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:2*pi]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void StartSquare(void){cout<<"plot 0,0 ";}
void StartFixedMasses(double x0Break,double y0Break,double r0){
  int i;
  for(i=0;i<Nx;i++){cout<<", "<<x0Break*i <<"+"<<r0<<"*cos(t),"<<-y0Break  <<"+"<<r0<<"*sin(t)";} /*Down line*/
  for(i=0;i<Nx;i++){cout<<", "<<x0Break*i <<"+"<<r0<<"*cos(t),"<<Ny*y0Break<<"+"<<r0<<"*sin(t)";} /*Top line*/
  for(i=0;i<Ny;i++){cout<<", "<<-x0Break  <<"+"<<r0<<"*cos(t),"<<y0Break*i <<"+"<<r0<<"*sin(t)";} /*Left line*/
  for(i=0;i<Ny;i++){cout<<", "<<Nx*x0Break<<"+"<<r0<<"*cos(t),"<<y0Break*i <<"+"<<r0<<"*sin(t)";} /*Right line*/
  /*Corners*/
  cout<<", "<<-x0Break  <<"+"<<r0<<"*cos(t),"<<-y0Break  <<"+"<<r0<<"*sin(t)";/*Lower Left*/
  cout<<", "<<Nx*x0Break<<"+"<<r0<<"*cos(t),"<<-y0Break  <<"+"<<r0<<"*sin(t)";/*Lower Right*/
  cout<<", "<<-x0Break  <<"+"<<r0<<"*cos(t),"<<Ny*y0Break<<"+"<<r0<<"*sin(t)";/*Upper Left*/
  cout<<", "<<Nx*x0Break<<"+"<<r0<<"*cos(t),"<<Ny*y0Break<<"+"<<r0<<"*sin(t)";/*Upper Right*/
}
void FinishSquare(void){cout<<endl;}
/*------------------------------------------------INTEGRATOR------------------------------------------------------------------*/
void Integrator(Body *Spring,Collider Hooke,double dt){
  /*Move with Omelyan PEFRL*/
  int i;
  for(i=0;i<N2D+1;i++){Spring[i].Move_R(dt,Xi);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_Theta(dt,Xi);}
  Hooke.AllForces(Spring,dt);
  for(i=0;i<N2D+1;i++){Spring[i].Move_V(dt,(1-2*Lambda)/2);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_W(dt,(1-2*Lambda)/2);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_R(dt,Chi);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_Theta(dt,Chi);}
  Hooke.AllForces(Spring,dt);
  for(i=0;i<N2D+1;i++){Spring[i].Move_V(dt,Lambda);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_W(dt,Lambda);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_R(dt,1-2*(Chi+Xi));}
  for(i=0;i<N2D+1;i++){Spring[i].Move_Theta(dt,1-2*(Chi+Xi));}
  Hooke.AllForces(Spring,dt);
  for(i=0;i<N2D+1;i++){Spring[i].Move_V(dt,Lambda);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_W(dt,Lambda);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_R(dt,Chi);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_Theta(dt,Chi);}
  Hooke.AllForces(Spring,dt);
  for(i=0;i<N2D+1;i++){Spring[i].Move_V(dt,(1-2*Lambda)/2);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_W(dt,(1-2*Lambda)/2);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_R(dt,Xi);}
  for(i=0;i<N2D+1;i++){Spring[i].Move_Theta(dt,Xi);}
}
/*------------------------------------------------MAIN PROGRAM----------------------------------------------------------------*/
int main(void)
{
  int i,j,M=Ny*(Nx-1)+Nx*(Ny-1);
  double t;                  /*Times*/
  double tdrawings,Ndrawings;/*Accountants*/
  Body Spring[N2D+1];        /*Oscillators Masses*/
  Oscillator SpringK[M];     /*Oscillators Springs*/
  Collider Hooke;            /*Collider*/

  double h=50*Ny;
  double x0Break=Lx, y0Break=Ly, z0Break=0;                                          /*Distance of Separation in each axis*/
  double m0=2*m,r0=35*L,Rx0=0,Ry0=0,Rz0=0,Vx0=0,Vy0=0,Vz0=0;                         /*Initial Conditions*/
  double Theta0=0,Phi0=0,Psi0=0,Wx0=0,Wy0=0,Wz0=0;                                   /*Initial Conditions*/
  double m1=1*m,r1=200*L,Rx1=Nx*x0Break/2,Ry1=Ny*y0Break+h,Rz1=0,Vx1=10*v,Vy1=0 ,Vz1=0;/*Initial Conditions Sphere*/
  double Theta1=0,Phi1=0,Psi1=0,Wx1=0,Wy1=0,Wz1=0.4/dt;                                 /*Initial Conditions Sphere*/
  double omegax=sqrt(Kx/m0),omegay=sqrt(Ky/m0),omegaz=sqrt(Kz/m0);                   /*One of Frequencies in each axis*/
  double Tx=2*M_PI/omegax, tmax=15*Tx;                                                /*Time Step and Tmax*/

  /*Start Springs*/
  for(i=0;i<Nx;i++){for(j=0;j<Ny;j++){
      Spring[i+Nx*j].Start(Rx0,Ry0,Rz0,Vx0,Vy0,Vz0,m0,r0,i*x0Break,j*y0Break,z0Break,Theta0,Psi0,Phi0,Wx0,Wy0,Wz0);}}
  /*Start Sphere*/
  Spring[N2D].Start(Rx1,Ry1,Rz1,Vx1,Vy1,Vz1,m1,r1,0,0,0,Theta1,Phi1,Psi1,Wx1,Wy1,Wz1);
  /*Start Collider*/
  Hooke.Start();
  /*Start Animation*/
  StartAnimation(x0Break,y0Break,h);  Ndrawings=1000;
  /*Integration*/
  for(t=tdrawings=0;t<tmax;t+=dt,tdrawings+=dt)
    {
      
      if(tdrawings>tmax/Ndrawings){
	StartSquare();
	//StartFixedMasses(x0Break,y0Break,r0);
	for(i=0;i<N2D+1;i++){Spring[i].Draw();}
	//for(j=0;j<Ny;j++){for(i=0;i<Nx;i++){SpringK[i].DrawInteractionSpring(Spring[i+j*Nx],Spring[i+1+j*Nx]);}}
	FinishSquare();
	tdrawings=0;}
      
      Integrator(Spring,Hooke,dt);
    }
  return 0;
}
