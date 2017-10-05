#include<iostream>
#include<fstream>
#include<cmath>
#include"Vector.h"
#include"Random64.h"

using namespace std;

const double g=10, K=1e4;
const double Kx=1,Ky=1,Kz=1;               /*Spring's Constant in each axis*/
const double Lx=100,Ly=100,Lz=100;         /*Lenght for parametrize the sides of the box*/
const int Nx=20,Ny=20,Nz=1;                /*Number of springs in each axis*/
const int N2D=Nx*Ny ,N3D=Nx*Ny*Nz;         /*Total number of springs*/
const double Gamma=50, Kcundall=10, MU=0.4;/*Contact Force*/

const double Xi=0.1786178958448091;
const double Lambda=-0.2123418310626054;
const double Chi=-0.06626458266981849;

class Body;
class Collider;
class Oscillator;
/*--------------------------------------------------CLASS BODY------------------------------------------------------------------------*/
class Body{
private:
  vector3D W,Tau; double Theta,I;
  double Rx,Ry,Rz,Vx,Vy,Vz,Fx,Fy,Fz,m,r,xBreak,yBreak,zBreak; 
public:
  void Start(double Rx0,double Ry0,double Rz0,double Vx0,double Vy0,double Vz0,
	     double m0,double r0,double x0Break,double y0Break,double z0Break,
	     double Theta0,double W0);
  /*Move the Position and Velocity in each component*/
  void Move_x(double dt, double Constante);
  void Move_y(double dt, double Constante);
  void Move_z(double dt, double Constante);
  void Move_Vx(double dt, double Constante);
  void Move_Vy(double dt, double Constante);
  void Move_Vz(double dt, double Constante);
  /*Inline Functions*/
  double Getx(void){return Rx;};
  double Gety(void){return Ry;};
  double Getz(void){return Rz;};
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
{Rx=Rx0; Ry=Ry0; Rz=Rz0; Vx=Vx0; Vy=Vy0; Vz=Vz0; m=m0; r=r0; xBreak=x0Break; yBreak=y0Break; zBreak=z0Break;
  Theta=Theta0;W.cargue(0,0,W0);I=2.0/5*m*r*r;}
void Body::Move_x(double dt, double Constante){Rx+=Vx*(Constante*dt);}
void Body::Move_y(double dt, double Constante){Ry+=Vy*(Constante*dt);}
void Body::Move_z(double dt, double Constante){Rz+=Vz*(Constante*dt);}
void Body::Move_Vx(double dt, double Constante){Vx+=Fx*(Constante*dt/m);}
void Body::Move_Vy(double dt, double Constante){Vy+=Fy*(Constante*dt/m);}
void Body::Move_Vz(double dt, double Constante){Vz+=Fz*(Constante*dt/m);}
void Body::EraseForce(void){Fx=0;Fy=0;Fz=0;}
void Body::AddForce_x(double F0x){Fx+=F0x;}
void Body::AddForce_y(double F0y){Fy+=F0y;}
void Body::AddForce_z(double F0z){Fz+=F0z;}
void Body::Draw(void){cout<<", "<<xBreak+Rx<<"+"<<r<<"*cos(t),"<<yBreak+Ry<<"+"<<r<<"*sin(t), "
			  <<Rx+xBreak<<"+"<<r*cos(Theta)/(2*M_PI)<<"*t,"<<Ry+yBreak<<"+"<<r*sin(Theta)/(2*M_PI)<<"*t";}
void Body::Draw3D(void){cout<<", "<<Rx<<"+"<<r<<"sin(v)*cos(t),"<<Ry<<"+"<<r<<"*sin(v)*sin(u),"<<Rz<<"+"<<r<<"*cos(v)";}
void Body::Move_Theta(double dt, double Constante){Theta+=W.z()*(Constante*dt);}
void Body::Move_W(double dt, double Constante){W+=Tau*(Constante*dt/I);}
void Body::EraseTorsion(void){Tau.cargue(0,0,0);}
void Body::AddForce(vector3D F0){AddForce_x(F0.x());AddForce_y(F0.y());AddForce_z(F0.z());}
void Body::AddTorsion(vector3D Tau0){Tau+=Tau0;}
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
  for(i=0;i<N2D+1;i++){Spring[i].AddForce_y(-Spring[i].m*g);}
  /*Add forces on x for extremum springs in x axis*/
  for(i=0;i<Nx;i++){for(j=0;j<Ny;j++){
      Spring[0+j*Ny].AddForce_x(-Kx*Spring[0+j*Ny].Rx); Spring[Nx-1+j*Ny].AddForce_x(-Kx*Spring[Nx-1+j*Ny].Rx);}}
  /*Add forces on y for extremum springs in y axis*/
  for(i=0;i<Nx;i++){Spring[i].AddForce_y(-Ky*Spring[i].Ry); Spring[Nx*(Ny-1)+i].AddForce_y(-Ky*Spring[Nx*(Ny-1)+i].Ry);}
  /*Local Interaction*/
  for(i=0;i<N2D;i++){
    if(i%Nx==Nx-1){Spring[i].AddForce_x(0);InteractionForce(Spring[i],Spring[i+Nx]);}       /*Last y-line*/
    else if((N2D-i)<=Ny){Spring[i].AddForce_y(0);InteractionForce(Spring[i],Spring[i+1]);}  /*Last x-line*/
    else if(i==N2D-1){Spring[i].AddForce_x(0);Spring[i].AddForce_y(0);}                     /*Last Spring*/
    else{InteractionForce(Spring[i],Spring[i+1]);InteractionForce(Spring[i],Spring[i+Nx]);}}
  
}

void Collider::InteractionForce(Body & Spring1, Body & Spring2){
  double drx=Spring2.Rx-Spring1.Rx;
  double dry=Spring2.Ry-Spring1.Ry;
  double drz=Spring2.Rz-Spring1.Rz;
  double F1x,F1y,F1z; vector3D F1;
  F1x=Kx*drx;F1y=Ky*dry;F1z=Kz*drz; F1.cargue(F1x,F1y,F1z);
  Spring1.AddForce(F1); Spring2.AddForce(F1*(-1));
}
void Collider::InteractionForce2(Body & Spring1, Body & Spring2,vector3D & l,bool & InCollision,double dt){
  vector3D R21,n,V1,V2,Vc,Vcn,Vct,t,Fn,Ft,F2; 
  double R1x,R1y,R1z,R2x,R2y,R2z,r1,r2,d21,s,m1,m2,m12,V1x,V1y,V1z,V2x,V2y,V2z,componentVcn,componentFn,normVct,Ftmax,normFt;
  double ERFF=1e-8;
  R1x=Spring1.Rx; R2x=Spring2.Rx; V1x=Spring1.Vx; V2x=Spring2.Vx;
  R1y=Spring1.Ry; R2y=Spring2.Ry; V1y=Spring1.Vy; V2y=Spring2.Vy;
  R1z=Spring1.Rz; R2z=Spring2.Rz; V1z=Spring1.Vz; V2z=Spring2.Vz;
  R21.cargue(R2x-R1x,R2y-R1y,R2z-R1z); d21=norma(R21);  s=(Spring1.r+Spring2.r)-d21;
  V1.cargue(V1x,V1y,V1z); V2.cargue(V2x,V2y,V2z); 
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
    componentFn=K*pow(s,1.5); 
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

    /*Building Total Force*/
    F2=Fn+Ft;
    Spring2.AddForce(F2);      Spring2.AddTorsion((n*(-r2))^Ft);      
    Spring1.AddForce(F2*(-1)); Spring1.AddTorsion((n*r1)^(Ft*(-1))); 
    
    InCollision=true;
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
  D.cargue(Spring2.Rx-Spring1.Rx,Spring2.Ry-Spring1.Ry,Spring2.Rz-Spring1.Rz); d=norma(D);
  X.cargue(1,0,0);
  a=(D*X)/d; /*Angle between D and x-axis*/
  /*Cicloid Parametrization
  x=r*(t-sin(t))
  y=r*(1-cos(t))
  Active Rotation Matrix [cos -sin; sin cos]*/
  cout<<", "<<Spring1.Rx<<"+"<<d<<"*(t-sin(2*t))*"<<cos(a)<<"-"<<d<<"*(1-cos(2*t))*"<<sin(a)
      <<"," <<Spring1.Ry<<"+"<<d<<"*(t-sin(2*t))*"<<sin(a)<<"+"<<d<<"*(1-cos(2*t))*"<<cos(a);
}
/*-------------------------------------GLOBAL'S FUNCTIONS------------------------------------------------*/
void StartAnimation(void){
  //cout<<"set terminal gif animate"<<endl; 
  //cout<<"set output 'MySpring2D.gif'"<<endl;
  cout<<"unset key"<<endl;
  //cout<<"set xrange [-100:1100]"<<endl;
  //cout<<"set yrange [-900:600]"<<endl;
  cout<<"set autoscale"<<endl;
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
  
  double m0=1 ,r0=5 ,Rx0=Lx/8   ,Ry0=Lx/10  ,Rz0=0,Vx0=0,Vy0=0 ,Vz0=0,Theta0=0,W0=0; /*Initial Conditions*/
  double m1=10,r1=50,Rx1=Lx*Nx/4,Ry1=3*Ly*Ny/5,Rz1=0,Vx1=0,Vy1=0 ,Vz1=0,Theta1=0,W1=4; /*Initial Conditions*/
  double x0Break=Lx/2, y0Break=Ly/2, z0Break=0;                                      /*Distance of Separation in each axis*/
  double omegax=sqrt(Kx/m0),omegay=sqrt(Ky/m0),omegaz=sqrt(Kz/m0);                   /*One of Frequencies in each axis*/
  double Tx=2*M_PI/omegax, tmax=2*Tx;                                                /*Time Step and Tmax*/

  /*Start Springs*/
  for(i=0;i<Nx;i++){for(j=0;j<Ny;j++){Spring[i+Nx*j].Start(Rx0,Ry0,Rz0,Vx0,Vy0,Vz0,m0,r0,i*x0Break,j*y0Break,z0Break,Theta0,W0);}}
  /*Start Sphere*/
  Spring[N2D].Start(Rx1,Ry1,Rz1,Vx1,Vy1,Vz1,m1,r1,0,0,0,Theta1,W1);
  /*Start Collider*/
  Hooke.Start();
  /*Start Animation*/
  StartAnimation();  Ndrawings=1000;
  /*Integration*/
  for(t=tdrawings=0;t<tmax;t+=dt,tdrawings+=dt)
    {
      
      if(tdrawings>tmax/Ndrawings){
	StartSquare();
	//StartFixedMasses(x0Break,y0Break,r0);
	for(i=0;i<N2D+1;i++){Spring[i].Draw();}
	FinishSquare();
	tdrawings=0;}
      /*
      for(i=0;i<N2D;i++){cout<<t<<"\t"<<i<<"\t"<<Spring[i].Getx()<<"\t"<<Spring[i].Gety()<<endl;}
      */
      //Move with Omelyan PEFRL
      for(i=0;i<N2D+1;i++){Spring[i].Move_x(dt,Xi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_y(dt,Xi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_z(dt,Xi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Theta(dt,Xi);}
      Hooke.AllForces(Spring,dt);
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vx(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vy(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vz(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_W(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_x(dt,Chi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_y(dt,Chi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_z(dt,Chi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Theta(dt,Chi);}
      Hooke.AllForces(Spring,dt);
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vx(dt,Lambda);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vy(dt,Lambda);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vz(dt,Lambda);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_W(dt,Lambda);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_x(dt,1-2*(Chi+Xi));}
      for(i=0;i<N2D+1;i++){Spring[i].Move_y(dt,1-2*(Chi+Xi));}
      for(i=0;i<N2D+1;i++){Spring[i].Move_z(dt,1-2*(Chi+Xi));}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Theta(dt,1-2*(Chi+Xi));}
      Hooke.AllForces(Spring,dt);
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vx(dt,Lambda);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vy(dt,Lambda);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vz(dt,Lambda);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_W(dt,Lambda);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_x(dt,Chi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_y(dt,Chi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_z(dt,Chi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Theta(dt,Chi);}
      Hooke.AllForces(Spring,dt);
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vx(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vy(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Vz(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_W(dt,(1-2*Lambda)/2);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_x(dt,Xi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_y(dt,Xi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_z(dt,Xi);}
      for(i=0;i<N2D+1;i++){Spring[i].Move_Theta(dt,Xi);}
    }
  return 0;
}
