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
const int Nx=20,Ny=6,Nz=3;        /*Number of springs in each axis*/
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
  double m,r,I,I1,I2,I3;/*Intrinsic Variables*/
  double q0,q1,q2,q3;   /*Quaternions*/
  double Nx,Ny,Nz;      /*Torsions in Inertial System*/
public:
  void Start(double Rx0,double Ry0,double Rz0,double Vx0,double Vy0,double Vz0,
	     double m0,double r0,double x0Break,double y0Break,double z0Break,
	     double Theta0,double Phi0,double Psi0,double Wx0,double Wy0,double Wz0,
	     double Ix, double Iy, double Iz);
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
  double WX(void){return W.x();};
  double WY(void){return W.y();};
  double WZ(void){return W.z();};
  double GetTheta(void){return acos((q0*q0-q1*q1-q2*q2+q3*q3)/(q0*q0+q1*q1+q2*q2+q3*q3));};
  double GetPhi(void){return atan((q1*q3+q0*q2)/(q0*q1-q2*q3));};
  double GetPsi(void){return atan((q1*q3-q0*q2)/(q0*q1+q2*q3));};
  /*Rigid Body*/
  void CalculateXYZ(vector3D & R);
  void Rotate(double dt);
  void CalculateTorsion(void);
  friend class Collider;
  friend class Oscillator;
};
/*-------------------------------------------------BODY'S FUNCTIONS-------------------------------------------------------------------*/
void Body::Start(double Rx0,double Ry0,double Rz0,double Vx0,double Vy0,double Vz0,
		 double m0,double r0,double x0Break,double y0Break,double z0Break,
		 double Theta0,double Phi0,double Psi0,double Wx0,double Wy0,double Wz0,
		 double Ix, double Iy, double Iz){
  R.cargue(Rx0,Ry0,Rz0);V.cargue(Vx0,Vy0,Vz0);m=m0;r=r0;RBreak.cargue(x0Break,y0Break,z0Break);
  Theta.cargue(Theta0,Phi0,Psi0);W.cargue(Wx0,Wy0,Wz0);I=2.0/5*m*r*r;I1=Ix;I2=Iy;I3=Iz;
  /*Initial Position (in Quaternions)*/
  q0=cos(0.5*Theta0)*cos(0.5*(Phi0+Psi0));
  q1=sin(0.5*Theta0)*cos(0.5*(Phi0-Psi0));
  q2=sin(0.5*Theta0)*sin(0.5*(Phi0-Psi0));
  q3=cos(0.5*Theta0)*sin(0.5*(Phi0+Psi0));
}
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
  cout<<", "<<X()<<"+"<<r<<"*sin(v)*cos(u),"<<Y()<<"+"<<r<<"*sin(v)*sin(u),"<<Z()<<"+"<<r<<"*cos(v)";
  cout<<", "<<X()<<"+"<<r*sin(Theta.z())*cos(Theta.x())/(2*M_PI*M_PI)<<"*u*v,"
      <<Y()<<"+"<<r*sin(Theta.x())*sin(Theta.z())/(2*M_PI*M_PI)<<"*u*v,"
      <<Z()<<"+"<<r*cos(Theta.z())/(M_PI)<<"*v";
}
void Body::CalculateXYZ(vector3D & R){
  vector3D e1,e2,e3; e1.cargue(1,0,0); e2.cargue(0,1,0); e3.cargue(0,0,1);
  double xi=R*e1,yi=R*e2,zi=R*e3,x,y,z;
  /*Transform Non inertial system to inertial system (Rotate the axis of Body)*/
  x=(q0*q0+q1*q1-q2*q2-q3*q3)*xi+2*(q1*q2-q0*q3)*yi+2*(q1*q3-q0*q2)*zi;
  y=2*(q1*q2+q0*q3)*xi+(q0*q0-q1*q1+q2*q2-q3*q3)*yi+2*(q2*q3-q0*q1)*zi;
  z=2*(q1*q3+q0*q2)*xi+2*(q2*q3+q0*q1)*yi+(q0*q0-q1*q1-q2*q2+q3*q3)*zi;
  R.cargue(x,y,z);
}
void Body::Rotate(double dt){
  double q0old=q0,q1old=q1,q2old=q2,q3old=q3;
  double Wx=W.x(),Wy=W.y(),Wz=W.z();
  double Wxold=Wx,Wyold=Wy,Wzold=Wz;
  double Thetax,Thetay,Thetaz;
  /*Update quaternions from angular velocity*/
  q0+=dt*0.5*(-q1old*Wx-q2old*Wy-q3old*Wz);
  q1+=dt*0.5*( q0old*Wx-q3old*Wy+q2old*Wz);
  q2+=dt*0.5*( q3old*Wx+q0old*Wy-q1old*Wz);
  q3+=dt*0.5*(-q2old*Wx+q1old*Wy+q0old*Wz);
  /*Update Angular velocity from torsion (Euler equation's)*/
  Wx+=dt*((Nx/I1)+Wyold*Wzold*(I2-I3)/I1);
  Wy+=dt*((Ny/I2)+Wxold*Wzold*(I3-I1)/I2);
  Wz+=dt*((Nz/I3)+Wxold*Wyold*(I1-I2)/I3);
  /*Update Angles from Quaternions*/
  Thetax=Theta.x();Thetay=Theta.y();Thetaz=Theta.z();               /*(Theta,Phi,Psi)*/
  Thetax+=acos((q0*q0-q1*q1-q2*q2+q3*q3)/(q0*q0+q1*q1+q2*q2+q3*q3));/*Theta*/
  Thetay+=atan((q1*q3+q0*q2)/(q0*q1-q2*q3));                        /*Phi*/
  Thetaz+=atan((q1*q3-q0*q2)/(q0*q1+q2*q3));                        /*Psi*/
  Theta.cargue(Thetax,Thetay,Thetaz);
}
void Body::CalculateTorsion(void){
  /*Torsion in the system attached to the body. In this case for the Sphere (Matrix A)*/
  double Fx,Fy,Fz;
  Fx=F.x();Fy=F.y();Fz=F.z()-m*g;
  Ny=(q0*q0+q1*q1-q2*q2-q3*q3)*Fx*r +2*(q1*q2+q0*q3)*Fy*r           -2*(q1*q3-q0*q2)*Fz*r;
  Nx=2*(q1*q2-q0*q3)*Fx*r           +(q0*q0-q1*q1+q2*q2-q3*q3)*Fy*r +2*(q2*q3+q0*q1)*Fz*r;
  Nz=2*(q1*q3+q0*q2)*Fx*r           +2*(q2*q3-q0*q1)*Fy*r           +(q0*q0-q1*q1-q2*q2+q3*q3)*Fz*r;
}
/*------------------------------------------------CLASS COLLIDER--------------------------------------------------------------*/
class Collider{
private:
  vector3D l[N3D+1][N3D+1]; bool InCollision[N3D+1][N3D+1];
public:
  void Start(void);
  void AllForces(Body* Spring,double dt);
  void InteractionForceSpring(Body & Spring1, Body & Spring2);
  void InteractionForceSphere(Body & Spring1, Body & Spring2,vector3D & l,bool & InCollision,double dt);
};
/*--------------------------------------COLLIDER'S FUNCTIONS------------------------------------------------------------------------*/
void Collider::Start(void){int i,j; for(i=0;i<N3D+1;i++){for(j=i+1;j<N3D+1;j++){l[i][j].cargue(0,0,0); InCollision[i][j]=false;}}}
void Collider::AllForces(Body* Spring,double dt){
  int i,j,k;
  /*Erase all Forces and  Torsions*/
  for(i=0;i<N3D+1;i++){Spring[i].EraseForce();Spring[i].EraseTorsion();}
  /*Add Gravitational Force due by the Earth's Gravitational Field*/
  for(i=N3D;i<N3D+1;i++){Spring[i].AddForce_z(-Spring[i].m*g);}
  /*Add Viscous Force on the springs*/
  for(i=0;i<N3D;i++){Spring[i].AddForce(Spring[i].V*(-GammaSpring*Spring[i].m));}
  /*Add forces on y for extremum springs in z axis*/
  for(j=0;j<Ny;j++){for(i=0;i<Nx;i++){
      Spring[i+Nx*j].AddForce_z(-Kz*Spring[i+Nx*j].R.z());
      Spring[(Nz-1)*Ny*Nx+i+Nx*j].AddForce_z(-Kz*Spring[(Nz-1)*Ny*Nx+i+Nx*j].R.z());}}
  /*Local Interaction with periodic boundary conditions*/
  for(k=0;k<Nz;k++){for(j=0;j<Ny;j++){for(i=0;i<Nx;i++){
      if(k==Nz-1){/*Last z-line*/
	Spring[k*Nx*Ny+i+j*Nx].AddForce_z(0);
	InteractionForceSpring(Spring[k*Nx*Ny+j*Nx+(i%Nx)],Spring[k*Nx*Ny+j*Nx+((i+1)%Nx)]); /*x-axis*/
	InteractionForceSpring(Spring[k*Nx*Ny+(j%Ny)*Nx+i],Spring[k*Nx*Ny+((j+1)%Ny)*Nx+i]);}/*y-axis*/
      else{
	InteractionForceSpring(Spring[k*Nx*Ny+j*Nx+i],Spring[(k+1)*Nx*Ny+j*Nx+i]);           /*z-axis*/
	InteractionForceSpring(Spring[k*Nx*Ny+j*Nx+(i%Nx)],Spring[k*Nx*Ny+j*Nx+((i+1)%Nx)]); /*x-axis*/
	InteractionForceSpring(Spring[k*Nx*Ny+(j%Ny)*Nx+i],Spring[k*Nx*Ny+((j+1)%Ny)*Nx+i]);}/*y-axis*/
      }}}
  /*Sphere and Spring*/
  for(i=0;i<N3D;i++){for(j=i+1;j<N3D+1;j++){InteractionForceSphere(Spring[i],Spring[j],l[i][j],InCollision[i][j],dt);}}
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
/*------------------------------------------------CLASS OSCILLATOR------------------------------------------------------------*/
class Oscillator{
private:
  
public:
  void DrawInteractionSpring3D(Body & Spring1,Body & Spring2);
};
/*-------------------------------------OSCILLATOR'S FUNCTIONS--------------------------------------------*/
void Oscillator::DrawInteractionSpring3D(Body & Spring1,Body & Spring2){
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
void StartAnimation3D(double x0Break, double y0Break, double z0Break, double h){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'MySpring3DR3.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xlabel 'x'"<<endl;
  cout<<"set ylabel 'y'"<<endl;
  cout<<"set zlabel 'z'"<<endl;
  cout<<"set xrange ["<<-x0Break<<":"<<x0Break*Nx<<"]"<<endl;
  cout<<"set yrange ["<<-y0Break<<":"<<y0Break*Ny<<"]"<<endl;
  cout<<"set zrange ["<<-z0Break<<":"<<z0Break*Nz+h<<"]"<<endl;
  //cout<<"set autoscale"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set urange [0:2*pi]"<<endl;
  cout<<"set vrange [0:pi]"<<endl;
  cout<<"set isosamples 12, 12"<<endl;  
}
void StartSquare3D(void){cout<<"splot 0,0,0 ";}
void StartFixedMasses3D(double x0Break,double y0Break,double z0Break,double r0){
  int i,j;
  for(i=0;i<Nx;i++){for(j=0;j<Nz;j++){
      cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*j<<"+"<<r0<<"*cos(v)";}}/*Front Plane*/
  for(i=0;i<Nx;i++){for(j=0;j<Nz;j++){
      cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*Ny<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*j<<"+"<<r0<<"*cos(v)";}}/*BackPlane*/
  for(i=0;i<Ny;i++){for(j=0;j<Nz;j++){
      cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*j<<"+"<<r0<<"*cos(v)";}}/*Left Plane*/
  for(i=0;i<Ny;i++){for(j=0;j<Nz;j++){
      cout<<", "<<x0Break*Nx<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*j<<"+"<<r0<<"*cos(v)";}}/*RigPlane*/
  for(i=0;i<Nx;i++){for(j=0;j<Ny;j++){
      cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*j<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";}}/*Down Plane*/
  for(i=0;i<Nx;i++){for(j=0;j<Ny;j++){
      cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*j<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*Nz<<"+"<<r0<<"*cos(v)";}}/*TopPlane*/
  /*Lines Corners*/
  for(i=0;i<Nx;i++){
    cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";}/*Down Front*/
  for(i=0;i<Nx;i++){
    cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*Ny<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";}/*Down Back*/
  for(i=0;i<Ny;i++){
    cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";}/*Down Left*/
  for(i=0;i<Ny;i++){
    cout<<", "<<x0Break*Nx<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";}/*Down Right*/
  for(i=0;i<Nx;i++){
    cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*Nz<<"+"<<r0<<"*cos(v)";}/*Top Front*/
  for(i=0;i<Nx;i++){
    cout<<", "<<x0Break*i<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*Ny<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*Nz<<"+"<<r0<<"*cos(v)";}/*Top Back*/
  for(i=0;i<Ny;i++){
    cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*Nz<<"+"<<r0<<"*cos(v)";}/*Top Left*/
  for(i=0;i<Ny;i++){
    cout<<", "<<x0Break*Nx<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*i<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*Nz<<"+"<<r0<<"*cos(v)";}/*Top Right*/
  for(i=0;i<Nz;i++){
    cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*i<<"+"<<r0<<"*cos(v)";}/*Front Left*/
  for(i=0;i<Nz;i++){
    cout<<", "<<x0Break*Nx<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*i<<"+"<<r0<<"*cos(v)";}/*Front Right*/
  for(i=0;i<Nz;i++){
    cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*Ny<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*i<<"+"<<r0<<"*cos(v)";}/*Back Left*/
  for(i=0;i<Nz;i++){
    cout<<", "<<x0Break*Nx<<"+"<<r0<<"*cos(u)*sin(v),"<<y0Break*Ny<<"+"<<r0<<"*sin(u)*sin(v),"<<z0Break*i<<"+"<<r0<<"*cos(v)";}/*Back Right*/
  /*Corners*/
  cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";/*Lower Left Front*/
  cout<<", "<<Nx*x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";/*Lower Right Front*/
  cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<Ny*y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";/*Lower Left Back*/
  cout<<", "<<Nx*x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<Ny*y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<-z0Break<<"+"<<r0<<"*cos(v)";/*Lower RightBack*/
  cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<Nz*z0Break<<"+"<<r0<<"*cos(v)";/*Top Left Front*/
  cout<<", "<<Nx*x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<-y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<Nz*z0Break<<"+"<<r0<<"*cos(v)";/*Top Right Front*/
  cout<<", "<<-x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<Ny*y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<Nz*z0Break<<"+"<<r0<<"*cos(v)";/*Top Left Back*/
  cout<<", "<<Nx*x0Break<<"+"<<r0<<"*cos(u)*sin(v),"<<Ny*y0Break<<"+"<<r0<<"*sin(u)*sin(v),"<<Nz*z0Break<<"+"<<r0<<"*cos(v)";/*Top RightBack*/
}
void FinishSquare(void){cout<<endl;}
/*------------------------------------------------INTEGRATOR------------------------------------------------------------------*/
void Integrator(Body *Spring,Collider Hooke,double dt){
  /*Move with Omelyan PEFRL*/
  int i;
  for(i=0;i<N3D+1;i++){Spring[i].Move_R(dt,Xi);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_Theta(dt,Xi);}
  Hooke.AllForces(Spring,dt);
  for(i=0;i<N3D+1;i++){Spring[i].Move_V(dt,(1-2*Lambda)/2);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_W(dt,(1-2*Lambda)/2);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_R(dt,Chi);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_Theta(dt,Chi);}
  Hooke.AllForces(Spring,dt);
  for(i=0;i<N3D+1;i++){Spring[i].Move_V(dt,Lambda);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_W(dt,Lambda);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_R(dt,1-2*(Chi+Xi));}
  for(i=0;i<N3D+1;i++){Spring[i].Move_Theta(dt,1-2*(Chi+Xi));}
  Hooke.AllForces(Spring,dt);
  for(i=0;i<N3D+1;i++){Spring[i].Move_V(dt,Lambda);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_W(dt,Lambda);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_R(dt,Chi);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_Theta(dt,Chi);}
  Hooke.AllForces(Spring,dt);
  for(i=0;i<N3D+1;i++){Spring[i].Move_V(dt,(1-2*Lambda)/2);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_W(dt,(1-2*Lambda)/2);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_R(dt,Xi);}
  for(i=0;i<N3D+1;i++){Spring[i].Move_Theta(dt,Xi);}
}
void IntegratorRigidBody(Body & Sphere,double dt){
  Sphere.CalculateTorsion();
  Sphere.Rotate(dt);
  vector3D R; R.cargue(Sphere.X(),Sphere.Y(),Sphere.Z());
  Sphere.CalculateXYZ(R);
}
/*------------------------------------------------MAIN PROGRAM----------------------------------------------------------------*/
int main(void)
{
  int i,j,k,M=Nz*Ny*(Nx-1)+Nz*Nx*(Ny-1)+Nx*Ny*(Nz-1);
  double t;                  /*Time*/
  double tdrawings,Ndrawings;/*Accountants*/
  Body Spring[N3D+1];        /*Oscillators Masses*/
  Oscillator SpringK[M];     /*Oscillators Springs*/
  Collider Hooke;            /*Collider*/

  double h=50*Nz;
  double x0Break=Lx, y0Break=Ly, z0Break=Lz;                         /*Distance of Separation in each axis*/
  double m0=2*m,r0=35*L,Rx0=0,Ry0=0,Rz0=0,Vx0=0,Vy0=0,Vz0=0;         /*Initial Conditions*/
  double Theta0=0,Phi0=0,Psi0=0,Wx0=0,Wy0=0,Wz0=0;                   /*Initial Conditions*/
  double m1=1*m,r1=200*L,Vx1=10*v,Vy1=0 ,Vz1=0;                      /*Initial Conditions Sphere*/
  double Rx1=Nx*x0Break/2,Ry1=Ny*y0Break/2,Rz1=Nz*z0Break+h;         /*Initial Conditions Sphere*/
  double Theta1=0,Phi1=0,Psi1=0;                                     /*Initial Conditions Sphere*/
  double Wx1=0.1/dt,Wy1=0,Wz1=0.4/dt,I1=(2.0/5)*m1*r1*r1,I2=I1,I3=I1;/*Initial Conditions Sphere*/
  double omegax=sqrt(Kx/m0),omegay=sqrt(Ky/m0),omegaz=sqrt(Kz/m0);   /*One of Frequencies in each axis*/
  double Tx=2*M_PI/omegax, tmax=25*Tx;                               /*Time Step and Tmax*/

  /*Start Springs*/
  for(k=0;k<Nz;k++){for(j=0;j<Ny;j++){for(i=0;i<Nx;i++){
	Spring[i+Nx*j+Nx*Ny*k].Start(Rx0,Ry0,Rz0,Vx0,Vy0,Vz0,m0,r0,i*x0Break,j*y0Break,k*z0Break,Theta0,Psi0,Phi0,Wx0,Wy0,Wz0,0,0,0);}}}
  /*Start Sphere*/
  Spring[N3D].Start(Rx1,Ry1,Rz1,Vx1,Vy1,Vz1,m1,r1,0,0,0,Theta1,Psi1,Phi1,Wx1,Wy1,Wz1,I1,I2,I3);
  /*Start Collider*/
  Hooke.Start();
  /*Start Animation*/
  StartAnimation3D(x0Break,y0Break,z0Break,2*h);  Ndrawings=1000;
  /*Integration*/
  for(t=tdrawings=0;t<tmax;t+=dt,tdrawings+=dt)
    {
      
      if(tdrawings>tmax/Ndrawings){
	StartSquare3D();
	//StartFixedMasses3D(x0Break,y0Break,z0Break,r0);
	for(i=0;i<N3D+1;i++){Spring[i].Draw3D();}
	//for(j=0;j<Ny;j++){for(i=0;i<Nx;i++){SpringK[i].DrawInteractionSpring3D(Spring[i+j*Nx],Spring[i+1+j*Nx]);}}
	FinishSquare();
	tdrawings=0;}
      
      Integrator(Spring,Hooke,dt);
      IntegratorRigidBody(Spring[N3D],dt);
      //cout<<Spring[N3D].Z()<<endl;
    }
  return 0;
}
