#include<iostream>
#include<fstream>
#include<cmath>
#include"Vector.h"
#include"Random64.h"

using namespace std;

const double g=0.5, K=1;
const double Kx=2,Ky=2,Kz=1;               /*Spring's Constant in each axis*/
const double Lx=100,Ly=100,Lz=100;         /*Lenght for parametrize the sides of the box*/
const double Gamma=20, Kcundall=10, MU=0.4;/*Contact Force*/
const double Beta=5,Gamma_t=0.5;           /*Viscous Force over Spring*/
const int Nx=80,Ny=4,Nz=1;                 /*Number of springs in each axis*/
const int N2D=Nx*Ny ,N3D=Nx*Ny*Nz;         /*Total number of springs*/

const double Xi=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Chi=-0.06626458266981849;

class Body;
class Collider;
class Oscillator;
/*--------------------------------------------------CLASS BODY------------------------------------------------------------------------*/
class Body{
private:
  vector3D W,Tau; double m,r,Theta,I;
  vector3D R,V,F,RBreak;
public:
  void Start(double Rx0,double Ry0,double Rz0,double Vx0,double Vy0,double Vz0,
	     double m0,double r0,double x0Break,double y0Break,double z0Break,
	     double Theta0,double W0);
  /*Move the Position and Velocity in each component*/
  void Move_R(double dt, double Constante);
  void Move_V(double dt, double Constante);
  /*Inline Functions*/
  double X(void){return RBreak.x()+R.x();};
  double Y(void){return RBreak.y()+R.y();};
  double Z(void){return RBreak.z()+R.z();};
  double NormF(void){return norma(F);};
  /*Add/Erases Forces*/
  void EraseForce(void);
  void AddForce_x(double F0x);
  void AddForce_y(double F0y);
  void AddForce_z(double F0z);
  /*Move Rotational Variables*/
  void Move_Theta(double dt, double Constante);
  void Move_W(double dt, double Constante);
  /*Add/Erases Torsions*/
  void EraseTorsion(void);
  void AddForce(vector3D F0);
  void AddTorsion(vector3D Tau0);
  /*Draw the Trajectory of the Particle*/
  void Draw(void);
  void Draw3D(void);
  friend class Collider;
  friend class Oscillator;
};
/*-------------------------------------------------BODY'S FUNCTIONS-------------------------------------------------------------------*/
void Body::Start(double Rx0,double Ry0,double Rz0,double Vx0,double Vy0,double Vz0,
		 double m0,double r0,double x0Break,double y0Break,double z0Break,
		 double Theta0,double W0)
{R.cargue(Rx0,Ry0,Rz0);V.cargue(Vx0,Vy0,Vz0);m=m0;r=r0;RBreak.cargue(x0Break,y0Break,z0Break);
  Theta=Theta0;W.cargue(0,0,W0);I=2.0/5*m*r*r;}
void Body::Move_R(double dt, double Constante){R+=V*(Constante*dt);}
void Body::Move_V(double dt, double Constante){V+=F*(Constante*dt/m);}
void Body::EraseForce(void){F.cargue(0,0,0);}
void Body::EraseTorsion(void){Tau.cargue(0,0,0);}
void Body::AddForce_x(double F0x){vector3D F0; F0.cargue(F0x,0,0);F+=F0;}
void Body::AddForce_y(double F0y){vector3D F0; F0.cargue(0,F0y,0);F+=F0;}
void Body::AddForce_z(double F0z){vector3D F0; F0.cargue(0,0,F0z);F+=F0;}
void Body::Move_Theta(double dt, double Constante){Theta+=W.z()*(Constante*dt);}
void Body::Move_W(double dt, double Constante){W+=Tau*(Constante*dt/I);}
void Body::AddForce(vector3D F0){F+=F0;}
void Body::AddTorsion(vector3D Tau0){Tau+=Tau0;}
void Body::Draw(void){cout<<", "<<X()<<"+"<<r<<"*cos(t),"<<Y()<<"+"<<r<<"*sin(t), "
			  <<X()<<"+"<<r*cos(Theta)/(2*M_PI)<<"*t,"<<Y()<<"+"<<r*sin(Theta)/(2*M_PI)<<"*t";}
void Body::Draw3D(void){cout<<", "<<X()<<"+"<<r<<"sin(v)*cos(t),"<<Y()<<"+"<<r<<"*sin(v)*sin(u),"<<Z()<<"+"<<r<<"*cos(v)";}
/*--------------------------------------CLASS COLLIDER-------------------------------------------------------------------------------*/
class Collider{
private:
  vector3D l[N2D+1][N2D+1]; bool InCollision[N2D+1][N2D+1];
public:
  void Start(void);
  void AllForces(Body* Spring,double dt);
  void InteractionForce(Body & Spring1, Body & Spring2);
  void InteractionForce2(Body & Spring1, Body & Spring2,vector3D & l,bool & InCollision,double dt);
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
  for(i=0;i<N2D;i++){Spring[i].AddForce(Spring[i].V*(-Beta*Spring[i].m));}
  /*Add forces on y for extremum springs in y axis*/
  for(i=0;i<Nx;i++){Spring[i].AddForce_y(-Ky*Spring[i].R.y()); Spring[Nx*(Ny-1)+i].AddForce_y(-Ky*Spring[Nx*(Ny-1)+i].R.y());}
  /*Local Interaction*/
  for(j=0;j<Ny;j++){for(i=0;i<Nx;i++){
      if(j==Ny-1){Spring[i+j*Nx].AddForce_y(0);InteractionForce(Spring[j*Nx+(i%Nx)],Spring[j*Nx+((i+1)%Nx)]);}/*Last x-line*/
      else{InteractionForce(Spring[j*Nx+(i%Nx)],Spring[j*Nx+((i+1)%Nx)]);InteractionForce(Spring[i],Spring[i+Nx]);}}}
  /*Sphere and Spring*/
  for(i=0;i<N2D;i++){InteractionForce2(Spring[i],Spring[N2D],l[i][N2D],InCollision[i][N2D],dt);}
}

void Collider::InteractionForce(Body & Spring1, Body & Spring2){
  double drx=Spring2.R.x()-Spring1.R.x();
  double dry=Spring2.R.y()-Spring1.R.y();
  double drz=Spring2.R.z()-Spring1.R.z();
  double F1x,F1y,F1z; vector3D F1;
  F1.cargue(Kx*drx,Ky*dry,Kz*drz);
  Spring1.AddForce(F1); Spring2.AddForce(F1*(-1));
}
void Collider::InteractionForce2(Body & Spring1, Body & Spring2,vector3D & l,bool & InCollision,double dt){
  vector3D R21,n,V1,V2,Vc,Vcn,Vct,t,Fn,Ft,F2; 
  double R1x,R1y,R1z,R2x,R2y,R2z,r1,r2,d21,s,m1,m2,m12,componentVcn,componentFn,normVct,Ftmax,normFt;
  double ERFF=1e-8;
  R1x=Spring1.X(); R2x=Spring2.X();
  R1y=Spring1.Y(); R2y=Spring2.Y();
  R1z=Spring1.Z(); R2z=Spring2.Z(); 
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
    componentFn=K*s; 
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
    Ft-=m12*Gamma_t*Vct;
    
    /*Building Total Force*/
    F2=Fn+Ft;
    Spring2.AddForce(F2);      Spring2.AddTorsion((n*(-r2))^Ft);      
    Spring1.AddForce(F2*(-1)); Spring1.AddTorsion((n*r1)^(Ft*(-1))); 
    
    InCollision=true;
    //cout<<Spring2.NormF()<<"\t"<<Spring2.Tau.x()<<"\t"<<Spring2.Tau.y()<<"\t"<<Spring2.Tau.z()<<endl;
  }
  else if(InCollision==true){l.cargue(0,0,0); InCollision=false;}
}
/*-------------------------------------CLASS OSCILLATOR--------------------------------------------------*/
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
void StartAnimation(double x0Break){
  //cout<<"set terminal gif animate"<<endl; 
  //cout<<"set output 'MySpring2D.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange ["<<-x0Break<<":"<<x0Break*Nx<<"]"<<endl;
  //cout<<"set yrange [-900:600]"<<endl;
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
/*-------------------------------------MAIN PROGRAM------------------------------------------------------*/

int main(void)
{
  int i,j,M=Ny*(Nx-1)+Nx*(Ny-1);
  double t,dt=1e-1;          /*Times*/
  double tdrawings,Ndrawings;/*Accountants*/
  Body Spring[N2D+1];        /*Oscillators Masses*/
  Oscillator SpringK[M];     /*Oscillators Springs*/
  Collider Hooke;            /*Collider*/

  double x0Break=Lx/2, y0Break=Ly/2, z0Break=0;                                      /*Distance of Separation in each axis*/
  double m0=1 ,r0=18,Rx0=0   ,Ry0=0  ,Rz0=0,Vx0=0,Vy0=0 ,Vz0=0,Theta0=0,W0=0; /*Initial Conditions*/
  double m1=30,r1=80,Rx1=Nx*x0Break/2,Ry1=(Ny+2)*y0Break,Rz1=0,Vx1=0,Vy1=0 ,Vz1=0,Theta1=0,W1=4; /*Initial Conditions*/
  double omegax=sqrt(Kx/m0),omegay=sqrt(Ky/m0),omegaz=sqrt(Kz/m0);                   /*One of Frequencies in each axis*/
  double Tx=2*M_PI/omegax, tmax=5*Tx;                                                /*Time Step and Tmax*/

  /*Start Springs*/
  for(i=0;i<Nx;i++){for(j=0;j<Ny;j++){Spring[i+Nx*j].Start(Rx0,Ry0,Rz0,Vx0,Vy0,Vz0,m0,r0,i*x0Break,j*y0Break,z0Break,Theta0,W0);}}
  /*Start Sphere*/
  Spring[N2D].Start(Rx1,Ry1,Rz1,Vx1,Vy1,Vz1,m1,r1,0,0,0,Theta1,W1);
  /*Start Collider*/
  Hooke.Start();
  /*Start Animation*/
  StartAnimation(x0Break);  Ndrawings=1000;
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
      
      //for(i=0;i<N2D+1;i++){cout<<"\t"<<i<<"\t"<<Spring[i].Y()<<"\t"<<Spring[i].NormF()<<endl;}
      
      //Move with Omelyan PEFRL
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
  return 0;
}
