#include<wavefunction.h>
#include<hamop.h>
#include<fluid.h>
#include<grid.h>

#define THRESH 1e-6
#define OOS 1.0/6.0
#define TOT 2.0/3.0
#define FOT 5.0/3.0
#define SQRTT sqrt(2.0)

#define complex complex<double>

void wavefunction::nullify()
{
  for (long i = 0; i < wf_dim; i++)
    start[i] = complex(0.0,0.0); 
}


// ------------------ for other overloaded version see below ------------------
void wavefunction::init(grid g, int inittype, double width, double k, double time, FILE* filename, int output_of_interest)
{
  long xindex, yindex, zindex, index_i,index,i;
  double x,y,z;
  double offs=0.0;
  double realpart, imagpart;
  long ctrl_id;
  
  if (inittype==99)
    { 
      for (i=0; i<output_of_interest; i++) 
	{
	  cout << "Scanning output no. " << i << "... - ignoring ..." << endl;
	  for (zindex=0; zindex<g.ngps_z(); zindex++)
	    {
	      for (yindex=0; yindex<g.ngps_y(); yindex++)
		{
		  for (xindex=0; xindex<g.ngps_x(); xindex++)
		    {      
		      ctrl_id=fscanf(filename,"%lf %lf",&realpart,&imagpart);
		    };
		};
	    };
	};
      cout << "Now I'm storing!" << endl;
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      for (xindex=0; xindex<g.ngps_x(); xindex++)
		{      
		  index=g.index(xindex,yindex,zindex);
		  ctrl_id=fscanf(filename,"%lf %lf",&realpart,&imagpart);
		  //		  printf("%e %e\n",realpart,imagpart);
		  start[index]=complex(realpart,imagpart);
		};
	    };
	};
    };
}

  
void wavefunction::init(grid g, int inittype, double width, double k, double time)
{
  long xindex, yindex, zindex, index_i;
  double x,y,z;
  complex st;
  double offs=0.0;
  
  st=width+0.5*complex(0.0,1.0)*time/width;

  srand((long)(width));

  for (xindex=0; xindex<g.ngps_x(); xindex++)
    { 
      x=g.x(xindex);
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  y=g.y(yindex);
	  for (zindex=0; zindex<g.ngps_z(); zindex++)
	    {
	      z=g.z(zindex);
	      index_i=g.index(xindex,yindex,zindex);
	      switch (inittype)
		{
		case 0 :
		  start[index_i]=cos(3.1415926/width*x)
		    *cos(3.1415926/width*y)
		    *sin(2.0*3.1415926/width*z);
		  break;
		case 1 :
		  start[index_i]=exp(-x*x-y*y);
		  break;
		case 2 :
		  start[index_i]=exp(-0.5*x*x*width)*exp(-0.5*y*y*width)*z*exp(-0.5*z*z*width); // not normalized
		  break;
		case 3 :
		  start[index_i]=rand()/(RAND_MAX+1.0)-0.5;
		  break;
		case 4 :
		  start[index_i]=cos(3.1415926/width*x);
		  break;
		case 5 :
		  if (yindex==0)
		    {
		      //start[index_i]=xindex; 
		      start[index_i]=g.r(xindex)*exp(-g.r(xindex)*width);
		    }
		  else
		    {
		      start[index_i]=complex(0.0,0.0);
		    };
		  break;
		case 6 :
		  start[index_i]=g.r(xindex)*exp(-g.r(xindex));

		  break;
		case 7 : 
		  start[index_i]=exp(-((x+offs-k*time)*(x+offs-k*time)
				       +y*y+z*z)/(4.0*st*width)
				     +complex(0.0,1.0)*k
				     *(x+offs-0.5*k*time))
		    *pow(2.0*M_PI*st*st,(double)(-0.25*g.dimens()));
		  break;
		case 8 :
		  start[index_i]=cos(3.1415926/width*x)*cos(3.1415926/width*y);
		  break;
		case 9 :
		  start[index_i]=cos(3.1415926/width*x)*sin(2.0*3.1415926/width*y);
		  break;
		case 10 :
		  start[index_i]=exp(-((x+offs-k*time)*(x+offs-k*time)
				       +0.5*y*y+0.5*z*z)/(4.0*st*width)
				     +complex(0.0,1.0)*k
				     *(x+offs-0.5*k*time))*y  // <----
		    *pow(2.0*M_PI*st*st,(double)(-0.25*g.dimens()));
		  break;	     	
		case 11 :
		  if (yindex==0)
		    {
		      //start[index_i]=xindex; 
		      start[index_i]=g.r(xindex)*exp(-(g.r(xindex)-4.0)*(g.r(xindex)-4.0)*width);
		    }
		  else
		    {
		      start[index_i]=complex(0.0,0.0);
		    };
		  break;
		case 12 :
		  start[index_i]=exp(-width*((sqrt(x*x+y*y)-2.6)*(sqrt(x*x+y*y)-2.6)));
		  //start[index_i]=exp(-width*((x-6.71)*(x-6.71)+y*y)) + exp(-width*((x+6.71)*(x+6.71)+y*y));
		  break;
		case 20 :
		   if (yindex==1)
		    {
	  	      start[index_i]=0.25/sqrt(2.0*3.141592)*sqrt(8.0)*2.0*g.r(xindex)*g.r(xindex)*exp(-g.r(xindex));
		    }
		  else
		    {
		      start[index_i]=complex(0.0,0.0);
		    };
		  break;

		default :
		  start[index_i]=complex(0.0,0.0);
		};
	    };
	};
    };

}

wavefunction &wavefunction::operator=(const wavefunction &v)
{
    if (this != &v)
    {
       delete[] start;
       wf_dim  = v.wf_dim;
       start = new complex[wf_dim];
       for (long i = 0; i < wf_dim; i++)
           start[i] = v.start[i];
    }
    return *this;
}

wavefunction &wavefunction::operator=(const fluid &v)
{
       delete[] start;
       wf_dim  = v.wf_size();
       start = new complex[wf_dim];
       for (long i = 0; i < wf_dim; i++)
           start[i] = v[i];
   
    return *this;
}

wavefunction operator + (const wavefunction &v, const wavefunction &w )
{
  
  wavefunction temp(v.wf_size());
  for(long i=0; i<v.wf_size(); i++)
    {
      temp[i]=v[i]+w[i];
    };
    
  return temp;

}

wavefunction operator - (const wavefunction &v, const wavefunction &w )
{
  
  wavefunction temp(v.wf_size());
  for(long i=0; i<v.wf_size(); i++)
    {
      temp[i]=v[i]-w[i];
    };
    
  return temp;

}


wavefunction operator / (const wavefunction &v, const wavefunction &w )
{
  
  wavefunction temp(v.wf_size());
  for(long i=0; i<v.wf_size(); i++)
    {
      temp[i]=v[i]/w[i];
    };
    
  return temp;

}


wavefunction operator + (const wavefunction &v, const fluid &w )
{
   wavefunction temp(v.wf_size());
  for(long i=0; i<v.wf_size(); i++)
    {
      temp[i]=v[i]+w[i];
    };
  return temp;
}

wavefunction operator + (const fluid &w, const wavefunction &v )
{
  wavefunction temp(v.wf_size());
  for(long i=0; i<v.wf_size(); i++)
    {
      temp[i]=v[i]+w[i];
    };
  return temp;

}


wavefunction operator + (const wavefunction &v, double z )
{
  wavefunction temp(v.wf_size());
  for(long i=0; i<v.wf_size(); i++)
    {
      temp[i]=v[i]+z;
    };
  return temp;

}


wavefunction& wavefunction::operator *= (double z)
{
  for(long i=0; i<wf_dim; i++)
    start[i]=start[i]*z;
  return *this;
}

wavefunction& wavefunction::operator *= (complex z)
{
  for(long i=0; i<wf_dim; i++)
    start[i]=start[i]*z;
  return *this;
}

wavefunction operator * (double z, const wavefunction &v)
{
  wavefunction temp=v;
  return temp *= z;
}

wavefunction operator * (complex z, const wavefunction &v)
{
  wavefunction temp=v;
  return temp *= z;
}

wavefunction operator * (const wavefunction &v, double z)
{
  wavefunction temp=v;
  return temp *= z;
}

wavefunction operator * (const wavefunction &v, complex z)
{
  wavefunction temp=v;
  return temp *= z;
}

complex operator * (const wavefunction &v, const wavefunction &w )
{
  complex result(0.0,0.0);
  for(long i=0; i<v.wf_size(); i++)
    {
      result+=conj(v[i])*w[i];
    };
  return result;
}


ostream& operator<<(ostream& os, const wavefunction& v)
{  
  for(long i = 0; i < v.wf_size(); i++)
    {   os << real(v[i]) << " " << imag(v[i]) << endl;
    }
  return os;
}


istream& operator>>(istream& is, wavefunction& v)
{  
  double tmpre, tmpim;
  for(long i = 0; i < v.wf_size(); i++)
    {   
      is >> tmpre >> tmpim;
      v[i]=complex(tmpre,tmpim);
    }
  return is;
}



void wavefunction::regrid(grid g, grid g_small, const wavefunction &v)
{
  long xindex,yindex,zindex,index,index_small,xshift,yshift;
  
  xshift=g.offs_x()-g_small.offs_x();
  yshift=g.offs_y()-g_small.offs_y();

  if ((xshift+g_small.ngps_x()<=g.ngps_x()) && (yshift+g_small.ngps_y()<=g.ngps_y()))
    {

      for (xindex=0; xindex<g_small.ngps_x(); xindex++)
	{
	  for (yindex=0; yindex<g_small.ngps_y(); yindex++)
	    {
	      for (zindex=0; zindex<g_small.ngps_z(); zindex++)
		{
		  index=g.index(xindex+xshift,yindex+yshift,zindex);
		  index_small=g_small.index(xindex,yindex,zindex);
		  start[index]=v[index_small];
		};
	    };
	};
    }
  else
    {
      cout << "Ooops! There's a problem in regridding data! " << endl;
    };


}


// overloaded below
void wavefunction::propagate(complex timestep, 
			     double time, 
			     grid g, 
			     hamop hamil, 
			     int me, 
			     int vecpotflag,
			     const wavefunction &staticpot,
			     const wavefunction &tddftpot,
			     double charge)
{

  switch (g.dimens())
    {



    case 15:
      do_cn_step_x_muller(timestep,time,g,hamil,me,staticpot,tddftpot,0,0);
      break;


    };

}



void wavefunction::propagate(complex timestep, 
			     double time, 
			     grid g, 
			     hamop hamil, 
			     int me, 
			     int vecpotflag,
			     const wavefunction &staticpot_x,
 			     const wavefunction &staticpot_y,
			     const wavefunction &staticpot_xy,
			     double charge)
{

  switch (g.dimens())
    {



    case 16 :

      do_cn_step_xy_muller(timestep,time,g,hamil,me,staticpot_x,staticpot_y,staticpot_xy,0,0);
      break;

    };

}


void wavefunction::do_cn_step_x(complex timestep, double time, grid g, hamop hamil, int me, int vecpotflag, long yindex, long zindex, const fluid &wf_one, const wavefunction &wf_two)
{
  long xindex;
  long index, index_xp, index_xm;
  double x, y, z, halfvecpotvecpot, field;
  wavefunction rhs(g.ngps_x()*g.ngps_y()*g.ngps_z());
  double const_diag_part;
  complex vpx;
  wavefunction a(g.ngps_x());
  wavefunction b(g.ngps_x());
  wavefunction c(g.ngps_x());

  const_diag_part=1.0/(g.delt_x()*g.delt_x());


  y=g.y(yindex);
  z=g.z(zindex);

  vpx=-0.5/(g.delt_x()*g.delt_x())-0.5*complex(0.0,1.0)*hamil.vecpot_x(time,me)/g.delt_x();
		  
  //  halfvecpotvecpot=0.5*hamil.vecpot_x(time,me)*hamil.vecpot_x(time,me);
  halfvecpotvecpot=0.0;
  field=hamil.field(time,me);

  for (xindex=0; xindex<g.ngps_x(); xindex++)
    { 
      x=g.x(xindex);
      a[xindex]=-vpx*0.5*complex(0.0,1.0)*timestep;
      c[xindex]=-conj(vpx)*0.5*complex(0.0,1.0)*timestep;
      b[xindex]=1.0-(const_diag_part
		     +hamil.scalarpotx(x,y,z,time,me)
		     +hamil.dftpot(g,x,y,z,time,me,wf_one,wf_two)
		     - complex(0.0,1.0)*hamil.imagpot(xindex,yindex,zindex,time,g)
		     +x*field
		     +halfvecpotvecpot )
	*0.5*complex(0.0,1.0)*timestep;
    };
  

  xindex=0;
  index=g.index(xindex,yindex,zindex);
  index_xp=g.index(xindex+1,yindex,zindex);

  rhs[index]=a[xindex]*start[index_xp]+b[xindex]*start[index];

  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
    { 
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);

      rhs[index]=a[xindex]*start[index_xp]+b[xindex]*start[index]+c[xindex]*start[index_xm];
    };

  xindex=g.ngps_x()-1;
  index=g.index(xindex,yindex,zindex);
  index_xm=g.index(xindex-1,yindex,zindex);

  rhs[index]=c[xindex]*start[index_xm]+b[xindex]*start[index];

  for (xindex=0; xindex<g.ngps_x(); xindex++)
    { 
      x=g.x(xindex);
      a[xindex]=vpx*0.5*complex(0.0,1.0)*timestep;
      c[xindex]=conj(vpx)*0.5*complex(0.0,1.0)*timestep;
      b[xindex]=1.0+(const_diag_part
		     +hamil.scalarpotx(x,y,z,time,me)
		     +hamil.dftpot(g,x,y,z,time,me,wf_one,wf_two)
		     - complex(0.0,1.0)*hamil.imagpot(xindex,yindex,zindex,time,g)
		     +x*field
		     +halfvecpotvecpot )
	*0.5*complex(0.0,1.0)*timestep;
    };


  solve(a,b,c,rhs,g.ngps_x());

  
}




void wavefunction::do_cn_step_x_muller(complex timestep, double time, grid g, hamop hamil, int me, const wavefunction &staticpot, const wavefunction &tddftpot, long yindex, long zindex)
{
  long xindex;
  long index, index_xp, index_xm;
  double x, y, z, halfvecpotvecpot, vecpot;
  wavefunction rhsone(g.ngps_x()*g.ngps_y()*g.ngps_z());
  double oneoverhsquare;
  wavefunction a(g.ngps_x());
  wavefunction b(g.ngps_x());
  wavefunction c(g.ngps_x());	
  wavefunction imagitimestephalfoversixstaticpot(g.ngps_x());
  wavefunction fiveimagitimestephalfoverthreestaticpot(g.ngps_x());
  wavefunction imagitimestephalfoversixtddftpot(g.ngps_x());
  wavefunction fiveimagitimestephalfoverthreetddftpot(g.ngps_x());
  complex aa,bb,cc,imagitimestephalf,imagitimestephalfoversix,fiveimagitimestephalfoverthree;
  double aaa,bbb,ccc;
  double b_upperleft, b_lowerright;
  double vecpotwithprefactor;
  double lambda=sqrt(3.0)-2.0;
  double llambda=-sqrt(3.0)+2.0;
  complex imagi=complex(0.0,1.0);

  oneoverhsquare=1.0/(g.delt_x()*g.delt_x());
  imagitimestephalf=imagi*timestep/2.0;
  imagitimestephalfoversix=imagitimestephalf/6.0;
  fiveimagitimestephalfoverthree=5.0*imagitimestephalf/3.0;

  vecpotwithprefactor=real(timestep/(8.0*g.delt_x())*hamil.vecpot_x(time,me));

  imagitimestephalfoversixstaticpot=imagitimestephalfoversix*staticpot;
  fiveimagitimestephalfoverthreestaticpot=fiveimagitimestephalfoverthree*staticpot;
  imagitimestephalfoversixtddftpot=imagitimestephalfoversix*tddftpot;
  fiveimagitimestephalfoverthreetddftpot=fiveimagitimestephalfoverthree*tddftpot;

  y=g.y(yindex);
  z=g.z(zindex);

  // The matrix S_-(tau/2)
  aaa=OOS-vecpotwithprefactor;
  ccc=OOS+vecpotwithprefactor;
  bbb=TOT;

  
  // Calculate the rhs vector S_-(tau/2) *this
  xindex=0;
  index=g.index(xindex,yindex,zindex);
  index_xp=g.index(xindex+1,yindex,zindex);
  
  rhsone[index]=aaa*start[index_xp]
    +((4.0+lambda)/6.0-lambda*vecpotwithprefactor)*start[index];

  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
    { 
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);

      rhsone[index]=aaa*start[index_xp]+bbb*start[index]+ccc*start[index_xm];
    };

  xindex=g.ngps_x()-1;
  index=g.index(xindex,yindex,zindex);
  index_xm=g.index(xindex-1,yindex,zindex);

  rhsone[index]=ccc*start[index_xm]
    +((4.0+lambda)/6.0-llambda*vecpotwithprefactor)*start[index];


  // The matrix  S_+(tau/2)
  aaa=OOS+vecpotwithprefactor;
  ccc=(OOS-vecpotwithprefactor);
  bbb=TOT;

  b_upperleft=((4.0+lambda)/6.0+lambda*vecpotwithprefactor);
  b_lowerright=((4.0+lambda)/6.0+llambda*vecpotwithprefactor);

  solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,rhsone,g.ngps_x());

  // The constant part of the matrix L_-(tau/2)
  aa=-OOS-0.5*imagitimestephalf*oneoverhsquare;
  cc=aa;
  bb=-FOT+0.5*imagitimestephalf*(2.0*oneoverhsquare);

  // Calculate the rhs vector L_-(tau/2) *this
  xindex=0;
  index=g.index(xindex,yindex,zindex);
  index_xp=g.index(xindex+1,yindex,zindex);
  rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot[index_xp]+0.5*imagitimestephalfoversixtddftpot[index_xp])*start[index_xp]
    +(bb+0.5*fiveimagitimestephalfoverthreestaticpot[index]+0.5*fiveimagitimestephalfoverthreetddftpot[index])*start[index];

  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
    { 
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);
      rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot[index_xp]+0.5*imagitimestephalfoversixtddftpot[index_xp])*start[index_xp]
	+(bb+0.5*fiveimagitimestephalfoverthreestaticpot[index]+0.5*fiveimagitimestephalfoverthreetddftpot[index])*start[index]
	+(cc+0.5*imagitimestephalfoversixstaticpot[index_xm]+0.5*imagitimestephalfoversixtddftpot[index_xm])*start[index_xm];
    };


  xindex=g.ngps_x()-1;
  index=g.index(xindex,yindex,zindex);
  index_xm=g.index(xindex-1,yindex,zindex);
  rhsone[index]=(bb+0.5*fiveimagitimestephalfoverthreestaticpot[index]+0.5*fiveimagitimestephalfoverthreetddftpot[index])*start[index]
    +(cc+0.5*imagitimestephalfoversixstaticpot[index_xm]+0.5*imagitimestephalfoversixtddftpot[index_xm])*start[index_xm];


  // The matrix L_+(tau/2)
  aa=-OOS+0.5*imagitimestephalf*oneoverhsquare;
  cc=-OOS+0.5*imagitimestephalf*oneoverhsquare;
  bb=-FOT-0.5*imagitimestephalf*(2.0*oneoverhsquare);  

  xindex=0;
  x=g.x(xindex);
  index=g.index(xindex,yindex,zindex);
  index_xp=g.index(xindex+1,yindex,zindex);
  a[xindex]=aa-0.5*imagitimestephalfoversixstaticpot[index_xp]-0.5*imagitimestephalfoversixtddftpot[index_xp];
  c[xindex]=1.0; // not used
  b[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot[index]-0.5*fiveimagitimestephalfoverthreetddftpot[index];

  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
    { 
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);

      a[xindex]=aa-0.5*imagitimestephalfoversixstaticpot[index_xp]-0.5*imagitimestephalfoversixtddftpot[index_xp];
      c[xindex]=cc-0.5*imagitimestephalfoversixstaticpot[index_xm]-0.5*imagitimestephalfoversixtddftpot[index_xm];
      b[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot[index]-0.5*fiveimagitimestephalfoverthreetddftpot[index];

    };  

  xindex=g.ngps_x()-1;
  x=g.x(xindex);
  index=g.index(xindex,yindex,zindex);
  index_xm=g.index(xindex-1,yindex,zindex);
  a[xindex]=aa-0.5*imagitimestephalfoversix*hamil.scalarpotx(x+g.delt_x(),y,z,time,me)-0.5*imagitimestephalfoversix*(2.0*tddftpot[xindex-1]-tddftpot[xindex-2]);
  c[xindex]=cc-0.5*imagitimestephalfoversixstaticpot[index_xm]-0.5*imagitimestephalfoversixtddftpot[index_xm];
  b[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot[index]-0.5*fiveimagitimestephalfoverthreetddftpot[index];
  

  solve(a,b,c,rhsone,g.ngps_x());



  // The constant part of the matrix L_-(tau/2)
  aa=-OOS-0.5*imagitimestephalf*oneoverhsquare;
  cc=aa;
  bb=-FOT+0.5*imagitimestephalf*(2.0*oneoverhsquare);

  // Calculate the rhs vector L_-(tau/2) *this
  xindex=0;
  index=g.index(xindex,yindex,zindex);
  index_xp=g.index(xindex+1,yindex,zindex);
  rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot[index_xp]+0.5*imagitimestephalfoversixtddftpot[index_xp])*start[index_xp]
    +(bb+0.5*fiveimagitimestephalfoverthreestaticpot[index]+0.5*fiveimagitimestephalfoverthreetddftpot[index])*start[index];

  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
    { 
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);
      rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot[index_xp]+0.5*imagitimestephalfoversixtddftpot[index_xp])*start[index_xp]
	+(bb+0.5*fiveimagitimestephalfoverthreestaticpot[index]+0.5*fiveimagitimestephalfoverthreetddftpot[index])*start[index]
	+(cc+0.5*imagitimestephalfoversixstaticpot[index_xm]+0.5*imagitimestephalfoversixtddftpot[index_xm])*start[index_xm];
    };


  xindex=g.ngps_x()-1;
  index=g.index(xindex,yindex,zindex);
  index_xm=g.index(xindex-1,yindex,zindex);
  rhsone[index]=(bb+0.5*fiveimagitimestephalfoverthreestaticpot[index]+0.5*fiveimagitimestephalfoverthreetddftpot[index])*start[index]
    +(cc+0.5*imagitimestephalfoversixstaticpot[index_xm]+0.5*imagitimestephalfoversixtddftpot[index_xm])*start[index_xm];


  // The matrix L_+(tau/2)
  aa=-OOS+0.5*imagitimestephalf*oneoverhsquare;
  cc=-OOS+0.5*imagitimestephalf*oneoverhsquare;
  bb=-FOT-0.5*imagitimestephalf*(2.0*oneoverhsquare);  

  xindex=0;
  x=g.x(xindex);
  index=g.index(xindex,yindex,zindex);
  index_xp=g.index(xindex+1,yindex,zindex);
  a[xindex]=aa-0.5*imagitimestephalfoversixstaticpot[index_xp]-0.5*imagitimestephalfoversixtddftpot[index_xp];
  c[xindex]=1.0; // not used
  b[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot[index]-0.5*fiveimagitimestephalfoverthreetddftpot[index];

  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
    { 
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);

      a[xindex]=aa-0.5*imagitimestephalfoversixstaticpot[index_xp]-0.5*imagitimestephalfoversixtddftpot[index_xp];
      c[xindex]=cc-0.5*imagitimestephalfoversixstaticpot[index_xm]-0.5*imagitimestephalfoversixtddftpot[index_xm];
      b[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot[index]-0.5*fiveimagitimestephalfoverthreetddftpot[index];

    };  

  xindex=g.ngps_x()-1;
  x=g.x(xindex);
  index=g.index(xindex,yindex,zindex);
  index_xm=g.index(xindex-1,yindex,zindex);
  a[xindex]=aa-0.5*imagitimestephalfoversix*hamil.scalarpotx(x+g.delt_x(),y,z,time,me)-0.5*imagitimestephalfoversix*(2.0*tddftpot[xindex-1]-tddftpot[xindex-2]);
  c[xindex]=cc-0.5*imagitimestephalfoversixstaticpot[index_xm]-0.5*imagitimestephalfoversixtddftpot[index_xm];
  b[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot[index]-0.5*fiveimagitimestephalfoverthreetddftpot[index];
  

  solve(a,b,c,rhsone,g.ngps_x());



  // The matrix S_-(tau/2)
  vecpotwithprefactor=real(timestep/(8.0*g.delt_x())*hamil.vecpot_x(real(time+timestep),me));
  aaa=OOS-vecpotwithprefactor;
  ccc=OOS+vecpotwithprefactor;
  bbb=TOT;

  
  // Calculate the rhs vector S_-(tau/2) *this
  xindex=0;
  index=g.index(xindex,yindex,zindex);
  index_xp=g.index(xindex+1,yindex,zindex);
  
  rhsone[index]=aaa*start[index_xp]
    +((4.0+lambda)/6.0-lambda*vecpotwithprefactor)*start[index];

  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
    { 
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);

      rhsone[index]=aaa*start[index_xp]+bbb*start[index]+ccc*start[index_xm];
    };

  xindex=g.ngps_x()-1;
  index=g.index(xindex,yindex,zindex);
  index_xm=g.index(xindex-1,yindex,zindex);

  rhsone[index]=ccc*start[index_xm]
    +((4.0+lambda)/6.0-llambda*vecpotwithprefactor)*start[index];

  // The matrix  S_+(tau/2)
  aaa=OOS+vecpotwithprefactor;
  ccc=(OOS-vecpotwithprefactor);
  bbb=TOT;

  b_upperleft=((4.0+lambda)/6.0+lambda*vecpotwithprefactor);
  b_lowerright=((4.0+lambda)/6.0+llambda*vecpotwithprefactor);

  solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,rhsone,g.ngps_x());


}





void wavefunction::do_cn_step_xy_muller(complex timestep, double time, grid g, hamop hamil, int me, const wavefunction &staticpot_x, const wavefunction &staticpot_y, const wavefunction &staticpot_xy, long yindex, long zindex)
{
  long xindex;
  long index, index_xp, index_xm, index_yp, index_ym;
  double x, y, z, halfvecpotvecpot, vecpot;
  wavefunction rhsone(g.ngps_x()*g.ngps_y()*g.ngps_z());
  wavefunction rhsone_x(g.ngps_x());
  wavefunction rhstwo_x(g.ngps_x());
  wavefunction rhsone_y(g.ngps_y());
  wavefunction rhstwo_y(g.ngps_y());
  double oneoverhsquare;
  wavefunction ax(g.ngps_x());
  wavefunction bx(g.ngps_x());
  wavefunction cx(g.ngps_x());	
  wavefunction ay(g.ngps_y());
  wavefunction by(g.ngps_y());
  wavefunction cy(g.ngps_y());
  wavefunction imagitimestepstaticpot_xy(g.ngps_x()*g.ngps_y()*g.ngps_z());
  wavefunction imagitimestephalfoversixstaticpot_x(g.ngps_x());
  wavefunction fiveimagitimestephalfoverthreestaticpot_x(g.ngps_x());
  wavefunction imagitimestephalfoversixstaticpot_y(g.ngps_y());
  wavefunction fiveimagitimestephalfoverthreestaticpot_y(g.ngps_y());

  complex aa,bb,cc,imagitimestep,imagitimestephalf,imagitimestephalfoversix,fiveimagitimestephalfoverthree;
  double aaa,bbb,ccc;
  double b_upperleft, b_lowerright;
  double vecpotwithprefactor;
  double lambda=sqrt(3.0)-2.0;
  double llambda=-sqrt(3.0)+2.0;
  complex imagi=complex(0.0,1.0);

  z=g.z(zindex);

  imagitimestep=imagi*timestep;
  imagitimestephalf=imagi*timestep/2.0;
  imagitimestephalfoversix=imagitimestephalf/6.0;
  fiveimagitimestephalfoverthree=5.0*imagitimestephalf/3.0;

  imagitimestephalfoversixstaticpot_x=imagitimestephalfoversix*staticpot_x;
  fiveimagitimestephalfoverthreestaticpot_x=fiveimagitimestephalfoverthree*staticpot_x;
  imagitimestephalfoversixstaticpot_y=imagitimestephalfoversix*staticpot_y;
  fiveimagitimestephalfoverthreestaticpot_y=fiveimagitimestephalfoverthree*staticpot_y;
  imagitimestepstaticpot_xy=imagitimestep*staticpot_xy;


  // ============================= S_x ========================================

  vecpotwithprefactor=real(timestep/(8.0*g.delt_x())*hamil.vecpot_x(time,me));


  aaa=OOS-vecpotwithprefactor;
  ccc=OOS+vecpotwithprefactor;
  bbb=TOT;


  // Calculate the rhs vector S_-x *this

  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      xindex=0;
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
  
      rhsone[index]=aaa*start[index_xp]
	+((4.0+lambda)/6.0-lambda*vecpotwithprefactor)*start[index];
      
      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  
	  rhsone[index]=aaa*start[index_xp]+bbb*start[index]+ccc*start[index_xm];
	};
      
      xindex=g.ngps_x()-1;
      index=g.index(xindex,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);
      
      rhsone[index]=ccc*start[index_xm]
	+((4.0+lambda)/6.0-llambda*vecpotwithprefactor)*start[index];
    }

  // The matrix  S_+x
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhsone_x[xindex]=rhsone[index];
	};

      aaa=OOS+vecpotwithprefactor;
      ccc=(OOS-vecpotwithprefactor);
      bbb=TOT;

      b_upperleft=((4.0+lambda)/6.0+lambda*vecpotwithprefactor);
      b_lowerright=((4.0+lambda)/6.0+llambda*vecpotwithprefactor);

      rhstwo_x.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,rhsone_x,g.ngps_x());

      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  start[index]=rhstwo_x[xindex];
	};
    };


  // ============================= S_y ========================================


  vecpotwithprefactor=real(timestep/(8.0*g.delt_y())*hamil.vecpot_y(time,me));


  aaa=OOS-vecpotwithprefactor;
  ccc=OOS+vecpotwithprefactor;
  bbb=TOT;


  // Calculate the rhs vector S_-y *this

  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      yindex=0;
      index=g.index(xindex,yindex,zindex);
      index_yp=g.index(xindex,yindex+1,zindex);
  
      rhsone[index]=aaa*start[index_yp]
	+((4.0+lambda)/6.0-lambda*vecpotwithprefactor)*start[index];
      
      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_yp=g.index(xindex,yindex+1,zindex);
	  index_ym=g.index(xindex,yindex-1,zindex);
	  
	  rhsone[index]=aaa*start[index_yp]+bbb*start[index]+ccc*start[index_ym];
	};
      
      yindex=g.ngps_y()-1;
      index=g.index(xindex,yindex,zindex);
      index_ym=g.index(xindex,yindex-1,zindex);
      
      rhsone[index]=ccc*start[index_ym]
	+((4.0+lambda)/6.0-llambda*vecpotwithprefactor)*start[index];
    }


  // The matrix  S_+y
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhsone_y[yindex]=rhsone[index];
	};

      aaa=OOS+vecpotwithprefactor;
      ccc=(OOS-vecpotwithprefactor);
      bbb=TOT;

      b_upperleft=((4.0+lambda)/6.0+lambda*vecpotwithprefactor);
      b_lowerright=((4.0+lambda)/6.0+llambda*vecpotwithprefactor);

      rhstwo_y.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,rhsone_y,g.ngps_y());

      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  start[index]=rhstwo_y[yindex];
	};
    };


  // ==================================== W_x =================================

  oneoverhsquare=1.0/(g.delt_x()*g.delt_x());

  aa=-OOS-0.5*imagitimestephalf*oneoverhsquare;
  cc=aa;
  bb=-FOT+0.5*imagitimestephalf*(2.0*oneoverhsquare);

  // ---------- W_-x
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      xindex=0;
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot_x[xindex+1])*start[index_xp]
	+(bb+0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex])*start[index];
      
      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot_x[xindex+1])*start[index_xp]
	    +(bb+0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex])*start[index]
	    +(cc+0.5*imagitimestephalfoversixstaticpot_x[xindex-1])*start[index_xm];
	};
      
      
      xindex=g.ngps_x()-1;
      index=g.index(xindex,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);
      rhsone[index]=(bb+0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex])*start[index]
	+(cc+0.5*imagitimestephalfoversixstaticpot_x[xindex-1])*start[index_xm];

    };


  
  // --------------- W_+x
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      y=g.y(yindex);
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhsone_x[xindex]=rhsone[index];
	};

      aa=-OOS+0.5*imagitimestephalf*oneoverhsquare;
      cc=-OOS+0.5*imagitimestephalf*oneoverhsquare;
      bb=-FOT-0.5*imagitimestephalf*(2.0*oneoverhsquare);  

      xindex=0;
      x=g.x(xindex);
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      ax[xindex]=aa-0.5*imagitimestephalfoversixstaticpot_x[xindex+1];
      cx[xindex]=1.0; // not used
      bx[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex];
      
      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  
	  ax[xindex]=aa-0.5*imagitimestephalfoversixstaticpot_x[xindex+1];
	  cx[xindex]=cc-0.5*imagitimestephalfoversixstaticpot_x[xindex-1];
	  bx[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex];
	  
	};  
      
      xindex=g.ngps_x()-1;
      x=g.x(xindex);
      index=g.index(xindex,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);
      ax[xindex]=aa-0.5*imagitimestephalfoversix*hamil.scalarpotx(x+g.delt_x(),y,z,time,me);
      cx[xindex]=cc-0.5*imagitimestephalfoversixstaticpot_x[xindex-1];
      bx[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex];
  
      
      rhstwo_x.solve(ax,bx,cx,rhsone_x,g.ngps_x());

      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  start[index]=rhstwo_x[xindex];
	};

    };



  // ==================================== W_y =================================

  oneoverhsquare=1.0/(g.delt_y()*g.delt_y());

  aa=-OOS-0.5*imagitimestephalf*oneoverhsquare;
  cc=aa;
  bb=-FOT+0.5*imagitimestephalf*(2.0*oneoverhsquare);

  // ---------- W_-y
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      yindex=0;
      index=g.index(xindex,yindex,zindex);
      index_yp=g.index(xindex,yindex+1,zindex);
      rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot_y[yindex+1])*start[index_yp]
	+(bb+0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex])*start[index];
      
      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_yp=g.index(xindex,yindex+1,zindex);
	  index_ym=g.index(xindex,yindex-1,zindex);
	  rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot_y[yindex+1])*start[index_yp]
	    +(bb+0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex])*start[index]
	    +(cc+0.5*imagitimestephalfoversixstaticpot_y[yindex-1])*start[index_ym];
	};
      
      
      yindex=g.ngps_y()-1;
      index=g.index(xindex,yindex,zindex);
      index_ym=g.index(xindex,yindex-1,zindex);
      rhsone[index]=(bb+0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex])*start[index]
	+(cc+0.5*imagitimestephalfoversixstaticpot_y[yindex-1])*start[index_ym];

    };


  
  // --------------- W_+y
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      x=g.x(xindex);
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhsone_y[yindex]=rhsone[index];
	};

      aa=-OOS+0.5*imagitimestephalf*oneoverhsquare;
      cc=-OOS+0.5*imagitimestephalf*oneoverhsquare;
      bb=-FOT-0.5*imagitimestephalf*(2.0*oneoverhsquare);  

      yindex=0;
      y=g.y(yindex);
      index=g.index(xindex,yindex,zindex);
      index_yp=g.index(xindex,yindex+1,zindex);
      ay[yindex]=aa-0.5*imagitimestephalfoversixstaticpot_y[yindex+1];
      cy[yindex]=1.0; // not used
      by[yindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex];
      
      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_yp=g.index(xindex,yindex+1,zindex);
	  index_ym=g.index(xindex,yindex-1,zindex);
	  
	  ay[yindex]=aa-0.5*imagitimestephalfoversixstaticpot_y[yindex+1];
	  cy[yindex]=cc-0.5*imagitimestephalfoversixstaticpot_y[yindex-1];
	  by[yindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex];
	  
	};  
      
      yindex=g.ngps_y()-1;
      y=g.y(yindex);
      index=g.index(xindex,yindex,zindex);
      index_ym=g.index(xindex,yindex-1,zindex);
      ay[yindex]=aa-0.5*imagitimestephalfoversix*hamil.scalarpot(x,y+g.delt_y(),z,time,me);
      cy[yindex]=cc-0.5*imagitimestephalfoversixstaticpot_y[yindex-1];
      by[yindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex];
  
      
      rhstwo_y.solve(ay,by,cy,rhsone_y,g.ngps_y());

      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  start[index]=rhstwo_y[yindex];
	};

    };

  // ============================ exp(V(xy)) =================================


      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      start[index]=exp(-imagitimestepstaticpot_xy[index])*start[index];
	    };
	};

  // ==================================== W_x =================================

  oneoverhsquare=1.0/(g.delt_x()*g.delt_x());

  aa=-OOS-0.5*imagitimestephalf*oneoverhsquare;
  cc=aa;
  bb=-FOT+0.5*imagitimestephalf*(2.0*oneoverhsquare);

  // ---------- W_-x
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      xindex=0;
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot_x[xindex+1])*start[index_xp]
	+(bb+0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex])*start[index];
      
      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot_x[xindex+1])*start[index_xp]
	    +(bb+0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex])*start[index]
	    +(cc+0.5*imagitimestephalfoversixstaticpot_x[xindex-1])*start[index_xm];
	};
      
      
      xindex=g.ngps_x()-1;
      index=g.index(xindex,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);
      rhsone[index]=(bb+0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex])*start[index]
	+(cc+0.5*imagitimestephalfoversixstaticpot_x[xindex-1])*start[index_xm];

    };


  
  // --------------- W_+x
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      y=g.y(yindex);
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhsone_x[xindex]=rhsone[index];
	};

      aa=-OOS+0.5*imagitimestephalf*oneoverhsquare;
      cc=-OOS+0.5*imagitimestephalf*oneoverhsquare;
      bb=-FOT-0.5*imagitimestephalf*(2.0*oneoverhsquare);  

      xindex=0;
      x=g.x(xindex);
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      ax[xindex]=aa-0.5*imagitimestephalfoversixstaticpot_x[xindex+1];
      cx[xindex]=1.0; // not used
      bx[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex];
      
      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  
	  ax[xindex]=aa-0.5*imagitimestephalfoversixstaticpot_x[xindex+1];
	  cx[xindex]=cc-0.5*imagitimestephalfoversixstaticpot_x[xindex-1];
	  bx[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex];
	  
	};  
      
      xindex=g.ngps_x()-1;
      x=g.x(xindex);
      index=g.index(xindex,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);
      ax[xindex]=aa-0.5*imagitimestephalfoversix*hamil.scalarpotx(x+g.delt_x(),y,z,time,me);
      cx[xindex]=cc-0.5*imagitimestephalfoversixstaticpot_x[xindex-1];
      bx[xindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_x[xindex];
  
      
      rhstwo_x.solve(ax,bx,cx,rhsone_x,g.ngps_x());

      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  start[index]=rhstwo_x[xindex];
	};

    };


  // ==================================== W_y =================================

  oneoverhsquare=1.0/(g.delt_y()*g.delt_y());

  aa=-OOS-0.5*imagitimestephalf*oneoverhsquare;
  cc=aa;
  bb=-FOT+0.5*imagitimestephalf*(2.0*oneoverhsquare);

  // ---------- W_-y
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      yindex=0;
      index=g.index(xindex,yindex,zindex);
      index_yp=g.index(xindex,yindex+1,zindex);
      rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot_y[yindex+1])*start[index_yp]
	+(bb+0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex])*start[index];
      
      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_yp=g.index(xindex,yindex+1,zindex);
	  index_ym=g.index(xindex,yindex-1,zindex);
	  rhsone[index]=(aa+0.5*imagitimestephalfoversixstaticpot_y[yindex+1])*start[index_yp]
	    +(bb+0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex])*start[index]
	    +(cc+0.5*imagitimestephalfoversixstaticpot_y[yindex-1])*start[index_ym];
	};
      
      
      yindex=g.ngps_y()-1;
      index=g.index(xindex,yindex,zindex);
      index_ym=g.index(xindex,yindex-1,zindex);
      rhsone[index]=(bb+0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex])*start[index]
	+(cc+0.5*imagitimestephalfoversixstaticpot_y[yindex-1])*start[index_ym];

    };


  
  // --------------- W_+y
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      x=g.x(xindex);
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhsone_y[yindex]=rhsone[index];
	};

      aa=-OOS+0.5*imagitimestephalf*oneoverhsquare;
      cc=-OOS+0.5*imagitimestephalf*oneoverhsquare;
      bb=-FOT-0.5*imagitimestephalf*(2.0*oneoverhsquare);  

      yindex=0;
      y=g.y(yindex);
      index=g.index(xindex,yindex,zindex);
      index_yp=g.index(xindex,yindex+1,zindex);
      ay[yindex]=aa-0.5*imagitimestephalfoversixstaticpot_y[yindex+1];
      cy[yindex]=1.0; // not used
      by[yindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex];
      
      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_yp=g.index(xindex,yindex+1,zindex);
	  index_ym=g.index(xindex,yindex-1,zindex);
	  
	  ay[yindex]=aa-0.5*imagitimestephalfoversixstaticpot_y[yindex+1];
	  cy[yindex]=cc-0.5*imagitimestephalfoversixstaticpot_y[yindex-1];
	  by[yindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex];
	  
	};  
      
      yindex=g.ngps_y()-1;
      y=g.y(yindex);
      index=g.index(xindex,yindex,zindex);
      index_ym=g.index(xindex,yindex-1,zindex);
      ay[yindex]=aa-0.5*imagitimestephalfoversix*hamil.scalarpot(x,y+g.delt_y(),z,time,me);
      cy[yindex]=cc-0.5*imagitimestephalfoversixstaticpot_y[yindex-1];
      by[yindex]=bb-0.5*fiveimagitimestephalfoverthreestaticpot_y[yindex];
  
      
      rhstwo_y.solve(ay,by,cy,rhsone_y,g.ngps_y());

      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  start[index]=rhstwo_y[yindex];
	};

    };

  // ============================= S_x ========================================

  vecpotwithprefactor=real(timestep/(8.0*g.delt_x())*hamil.vecpot_x(time,me));


  aaa=OOS-vecpotwithprefactor;
  ccc=OOS+vecpotwithprefactor;
  bbb=TOT;


  // Calculate the rhs vector S_-x *this

  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      xindex=0;
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
  
      rhsone[index]=aaa*start[index_xp]
	+((4.0+lambda)/6.0-lambda*vecpotwithprefactor)*start[index];
      
      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  
	  rhsone[index]=aaa*start[index_xp]+bbb*start[index]+ccc*start[index_xm];
	};
      
      xindex=g.ngps_x()-1;
      index=g.index(xindex,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);
      
      rhsone[index]=ccc*start[index_xm]
	+((4.0+lambda)/6.0-llambda*vecpotwithprefactor)*start[index];
    }

  // The matrix  S_+x
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhsone_x[xindex]=rhsone[index];
	};

      aaa=OOS+vecpotwithprefactor;
      ccc=(OOS-vecpotwithprefactor);
      bbb=TOT;

      b_upperleft=((4.0+lambda)/6.0+lambda*vecpotwithprefactor);
      b_lowerright=((4.0+lambda)/6.0+llambda*vecpotwithprefactor);

      rhstwo_x.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,rhsone_x,g.ngps_x());

      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  start[index]=rhstwo_x[xindex];
	};
    };


  // ============================= S_y ========================================


  vecpotwithprefactor=real(timestep/(8.0*g.delt_y())*hamil.vecpot_y(time,me));


  aaa=OOS-vecpotwithprefactor;
  ccc=OOS+vecpotwithprefactor;
  bbb=TOT;


  // Calculate the rhs vector S_-y *this

  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      yindex=0;
      index=g.index(xindex,yindex,zindex);
      index_yp=g.index(xindex,yindex+1,zindex);
  
      rhsone[index]=aaa*start[index_yp]
	+((4.0+lambda)/6.0-lambda*vecpotwithprefactor)*start[index];
      
      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_yp=g.index(xindex,yindex+1,zindex);
	  index_ym=g.index(xindex,yindex-1,zindex);
	  
	  rhsone[index]=aaa*start[index_yp]+bbb*start[index]+ccc*start[index_ym];
	};
      
      yindex=g.ngps_y()-1;
      index=g.index(xindex,yindex,zindex);
      index_ym=g.index(xindex,yindex-1,zindex);
      
      rhsone[index]=ccc*start[index_ym]
	+((4.0+lambda)/6.0-llambda*vecpotwithprefactor)*start[index];
    }


  // The matrix  S_+y
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhsone_y[yindex]=rhsone[index];
	};

      aaa=OOS+vecpotwithprefactor;
      ccc=(OOS-vecpotwithprefactor);
      bbb=TOT;

      b_upperleft=((4.0+lambda)/6.0+lambda*vecpotwithprefactor);
      b_lowerright=((4.0+lambda)/6.0+llambda*vecpotwithprefactor);

      rhstwo_y.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,rhsone_y,g.ngps_y());

      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  start[index]=rhstwo_y[yindex];
	};
    };



}



void wavefunction::do_muller(complex timestep, 
			     double time, 
			     grid g, 
			     hamop hamil, 
			     const wavefunction &staticpot, 
			     int me, 
			     long noofgridpoints[], 
			     double charge)
{

  double Eff_Charge = charge;

  long xindex, yindex;
  long index, index_xp, index_xm, index_lp, max_noofgridpoints;
  double r;
  wavefunction rhsone(g.ngps_x()*g.ngps_y()*g.ngps_z());
  double wfsq;
  wavefunction aa(g.ngps_x());
  wavefunction bb(g.ngps_x());
  wavefunction cc(g.ngps_x());	
  wavefunction tmp_plus(g.ngps_x());	
  wavefunction tmp_minus(g.ngps_x());	
  wavefunction tmp_plus_two(g.ngps_x());	
  wavefunction tmp_minus_two(g.ngps_x());	
  complex imagi(0.0,1.0);
  complex Delta_two_upperleft, M_two_upperleft;
  double c, cl, cnl, tl, ctilde, fnl, fnl_nom, fnlfnl;
  double aaa,bbb,ccc,Mtilde;
  double b_upperleft, b_lowerright;
  double factor;
  double vecpotwithprefactor;
  double lambda=sqrt(3.0)-2.0;
  double llambda=-sqrt(3.0)+2.0;
  double ul,ur,ll,lr,maxpsi;
  complex aaaa, bbbb, cccc;

  complex halfimagitimestep=0.5*timestep*imagi;
  complex halfimagitimestepOOS=halfimagitimestep*OOS;
  complex halfimagitimestepFOT=halfimagitimestep*FOT;

  // determine the grid size (in radial direction) for each l
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {  
      if (noofgridpoints[yindex]<g.ngps_x())
	{
	  maxpsi=0.0;
	  for (xindex=0; xindex<noofgridpoints[yindex]; xindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      wfsq=real(start[index]*conj(start[index]));
	      if (wfsq>maxpsi) maxpsi=wfsq;
	    };
	  if (wfsq>THRESH*maxpsi) 
	    {
	      noofgridpoints[yindex]+=50;
	      if (noofgridpoints[yindex]>g.ngps_x()) noofgridpoints[yindex]=g.ngps_x();
	    };
	};
    };

  cout << "noofgridpoints at l=0, l=lmax/2, l=lmax: " << noofgridpoints[0] << " " 
       << noofgridpoints[(long)(0.5*g.ngps_y())] << " " << noofgridpoints[g.ngps_y()-1] << "\n";




  // calculate one l-block
  for (yindex=0; yindex<g.ngps_y()-1; yindex++)
    {
      cnl=(yindex+1)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
      c=real(0.25*timestep*hamil.vecpot_z(time,me)*cnl/(2.0*g.delt_x()));
      // calculate W_nl Psi
      tl=(yindex+1)*(yindex+1)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
      fnl_nom=real(0.25*timestep*hamil.vecpot_z(time,me)*tl);
      max_noofgridpoints=noofgridpoints[yindex+1]>noofgridpoints[yindex] ? noofgridpoints[yindex+1] : noofgridpoints[yindex];
      for (xindex=0; xindex<max_noofgridpoints; xindex++)
	{      
	  index=g.index(xindex,yindex,0);
	  index_lp=g.index(xindex,yindex+1,0);
	  fnl=fnl_nom/g.r(xindex);
	  fnlfnl=fnl*fnl;
	  factor=1.0/(2.0*(1.0+fnlfnl));
	  ul=1.0+2.0*fnl-fnlfnl;
	  ur=1.0-2.0*fnl-fnlfnl;
	  tmp_plus[xindex]=(ul*start[index]+ur*start[index_lp])*factor;
	  tmp_minus[xindex]=(-ur*start[index]+ul*start[index_lp])*factor;
	};
      // calculate Y_l^- tmp
      xindex=0;
      Mtilde=TOT+OOS*lambda;
      ctilde=lambda*c;
      tmp_plus_two[xindex]=(Mtilde-ctilde)*tmp_plus[xindex] + (OOS-c)*tmp_plus[xindex+1];
      tmp_minus_two[xindex]=(Mtilde+ctilde)*tmp_minus[xindex] + (OOS+c)*tmp_minus[xindex+1];
      ul=OOS+c;
      ur=OOS-c;
      for (xindex=1; xindex<max_noofgridpoints-1; xindex++)
	{
	  tmp_plus_two[xindex]=ul*tmp_plus[xindex-1] + TOT*tmp_plus[xindex] + ur*tmp_plus[xindex+1];
	  tmp_minus_two[xindex]=ur*tmp_minus[xindex-1] + TOT*tmp_minus[xindex] + ul*tmp_minus[xindex+1];
	};
      xindex=max_noofgridpoints-1;
      Mtilde=TOT+OOS*lambda;
      ctilde=llambda*c;
      tmp_plus_two[xindex]=(OOS+c)*tmp_plus[xindex-1] + (Mtilde-ctilde)*tmp_plus[xindex];
      tmp_minus_two[xindex]=(OOS-c)*tmp_minus[xindex-1] + (Mtilde+ctilde)*tmp_minus[xindex];
      // solve for Y_l^+ tmp_new = tmp_old
      Mtilde=TOT+OOS*lambda;
      ctilde=lambda*c;
      b_upperleft=(Mtilde+ctilde);
      Mtilde=TOT+OOS*lambda;
      ctilde=llambda*c;
      b_lowerright=(Mtilde+ctilde);
      aaa=OOS+c;
      bbb=TOT;
      ccc=OOS-c;
      tmp_plus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_plus_two,max_noofgridpoints);
      Mtilde=TOT+OOS*lambda;
      ctilde=lambda*c;
      b_upperleft=(Mtilde-ctilde);
      Mtilde=TOT+OOS*lambda;
      ctilde=llambda*c;
      b_lowerright=(Mtilde-ctilde);
      aaa=OOS-c;
      bbb=TOT;
      ccc=OOS+c;
      tmp_minus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_minus_two,max_noofgridpoints);
      // update the wavefunction
      for (xindex=0; xindex<max_noofgridpoints; xindex++)
	{
	  index=g.index(xindex,yindex,0);
	  index_lp=g.index(xindex,yindex+1,0);
	  start[index]=(tmp_plus[xindex]-tmp_minus[xindex]);
	  start[index_lp]=(tmp_plus[xindex]+tmp_minus[xindex]);
	};
    };

  // the spatial part
  // The constant part of the matrix L_-(tau)
  aaaa=-OOS-halfimagitimestep/(g.delt_x()*g.delt_x());
  cccc=aaaa;
  bbbb=-FOT+imagi*timestep/(g.delt_x()*g.delt_x());
  Delta_two_upperleft=-2.0/(g.delt_x()*g.delt_x())*(1.0-Eff_Charge*g.delt_x()/(12.0-10.0*Eff_Charge*g.delt_x()));
  M_two_upperleft=-2.0*(1.0+g.delt_x()*g.delt_x()/12.0*Delta_two_upperleft);
  
  // Calculate the rhs vector L_-(tau) *this
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    { 
      xindex=0;
      index=g.index(xindex,yindex,0);
      index_xp=g.index(xindex+1,yindex,0);
      if (yindex==0)
	{
	  rhsone[index]=(aaaa+halfimagitimestepOOS*staticpot[index_xp])*start[index_xp]
	    +(M_two_upperleft*(1.0-halfimagitimestep*staticpot[index])-imagi*Delta_two_upperleft*0.5*timestep)*start[index];
	}
      else
	{
	  rhsone[index]=(aaaa+halfimagitimestepOOS*staticpot[index_xp])*start[index_xp]
	    +(bbbb+halfimagitimestepFOT*staticpot[index])*start[index];
	};
	  
      max_noofgridpoints=noofgridpoints[yindex];
      for (xindex=1; xindex<max_noofgridpoints-1; xindex++)
	{ 
	  index=g.index(xindex,yindex,0);
	  index_xp=g.index(xindex+1,yindex,0);
	  index_xm=g.index(xindex-1,yindex,0);
	  rhsone[index]=(aaaa+halfimagitimestepOOS*staticpot[index_xp])*start[index_xp]
	    +(bbbb+halfimagitimestepFOT*staticpot[index])*start[index]
	    +(cccc+halfimagitimestepOOS*staticpot[index_xm])*start[index_xm];
	};
      

      xindex=max_noofgridpoints-1;
      index=g.index(xindex,yindex,0);
      index_xm=g.index(xindex-1,yindex,0);
      rhsone[index]=(bbbb+halfimagitimestepFOT*staticpot[index])*start[index]
	+(cccc+halfimagitimestepOOS*staticpot[index_xm])*start[index_xm];
    };

  // The matrix L_+(tau)
  aaaa=-OOS+halfimagitimestep/(g.delt_x()*g.delt_x());
  cccc=aaaa;
  bbbb=-FOT-imagi*timestep/(g.delt_x()*g.delt_x());

  for (yindex=0; yindex<g.ngps_y(); yindex++)
    { 
      xindex=0;
      index=g.index(xindex,yindex,0);
      index_xp=g.index(xindex+1,yindex,0);
      aa[xindex]=(aaaa-halfimagitimestep*OOS*staticpot[index_xp]);
      cc[xindex]=1.0; // not used

      if (yindex==0)
	{
	  bb[xindex]=(M_two_upperleft*(1.0+halfimagitimestep*staticpot[index])+imagi*Delta_two_upperleft*0.5*timestep);
	}
      else
	{
	  bb[xindex]=(bbbb-halfimagitimestepFOT*staticpot[index]);
	};

      tmp_plus[xindex]=rhsone[index];

      max_noofgridpoints=noofgridpoints[yindex];
      for (xindex=1; xindex<max_noofgridpoints-1; xindex++)
	{ 
	  index=g.index(xindex,yindex,0);
	  index_xp=g.index(xindex+1,yindex,0);
	  index_xm=g.index(xindex-1,yindex,0);
	  
	  aa[xindex]=(aaaa-halfimagitimestepOOS*staticpot[index_xp]);
	  cc[xindex]=(cccc-halfimagitimestepOOS*staticpot[index_xm]);
	  bb[xindex]=(bbbb-halfimagitimestepFOT*staticpot[index]);
	  tmp_plus[xindex]=rhsone[index];
	};  
      
      xindex=max_noofgridpoints-1;
      r=g.r(xindex);
      index=g.index(xindex,yindex,0);
      index_xm=g.index(xindex-1,yindex,0);
      aa[xindex]=aaaa-halfimagitimestepOOS*(hamil.scalarpotx(r+g.delt_x(),yindex,0.0,time,me)
					   -imagi*hamil.imagpot(xindex,yindex,0,time,g));
      cc[xindex]=(cccc-halfimagitimestepOOS*staticpot[index_xm]);
      bb[xindex]=(bbbb-halfimagitimestepFOT*staticpot[index]);
      tmp_plus[xindex]=rhsone[index];
      
      tmp_minus.solve(aa,bb,cc,tmp_plus,max_noofgridpoints);
      
      for (xindex=0; xindex<max_noofgridpoints; xindex++)
	{ 
	  index=g.index(xindex,yindex,0);
	  start[index]=tmp_minus[xindex];
	};  
    };


  // calculate one l-block
  for (yindex=g.ngps_y()-2; yindex>=0; yindex--)
    {
      cnl=(yindex+1)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
      c=real(0.25*timestep*hamil.vecpot_z(time,me)*cnl/(2.0*g.delt_x()));
      // calculate B_l Psi
      max_noofgridpoints=noofgridpoints[yindex+1]>noofgridpoints[yindex] ? noofgridpoints[yindex+1] : noofgridpoints[yindex];
      for (xindex=0; xindex<max_noofgridpoints; xindex++)
	{      
	  index=g.index(xindex,yindex,0);
	  index_lp=g.index(xindex,yindex+1,0);
	  tmp_plus[xindex]=(start[index]+start[index_lp]);
	  tmp_minus[xindex]=(-start[index]+start[index_lp]);
	};
      // calculate Y_l^- tmp
      xindex=0;
      Mtilde=TOT+OOS*lambda;
      ctilde=lambda*c;
      tmp_plus_two[xindex]=(Mtilde-ctilde)*tmp_plus[xindex] + (OOS-c)*tmp_plus[xindex+1];
      tmp_minus_two[xindex]=(Mtilde+ctilde)*tmp_minus[xindex] + (OOS+c)*tmp_minus[xindex+1];
      ul=OOS+c;
      ur=OOS-c;
      for (xindex=1; xindex<max_noofgridpoints-1; xindex++)
	{
	  tmp_plus_two[xindex]=ul*tmp_plus[xindex-1] + TOT*tmp_plus[xindex] + ur*tmp_plus[xindex+1];
	  tmp_minus_two[xindex]=ur*tmp_minus[xindex-1] + TOT*tmp_minus[xindex] + ul*tmp_minus[xindex+1];
	};
      xindex=max_noofgridpoints-1;
      Mtilde=TOT+OOS*lambda;
      ctilde=llambda*c;
      tmp_plus_two[xindex]=(OOS+c)*tmp_plus[xindex-1] + (Mtilde-ctilde)*tmp_plus[xindex];
      tmp_minus_two[xindex]=(OOS-c)*tmp_minus[xindex-1] + (Mtilde+ctilde)*tmp_minus[xindex];
      // solve for Y_l^+ tmp_new = tmp_old
      Mtilde=TOT+OOS*lambda;
      ctilde=lambda*c;
      b_upperleft=(Mtilde+ctilde);
      Mtilde=TOT+OOS*lambda;
      ctilde=llambda*c;
      b_lowerright=(Mtilde+ctilde);
      aaa=OOS+c;
      bbb=TOT;
      ccc=OOS-c;
      tmp_plus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_plus_two,max_noofgridpoints);
      Mtilde=TOT+OOS*lambda;
      ctilde=lambda*c;
      b_upperleft=(Mtilde-ctilde);
      Mtilde=TOT+OOS*lambda;
      ctilde=llambda*c;
      b_lowerright=(Mtilde-ctilde);
      aaa=OOS-c;
      bbb=TOT;
      ccc=OOS+c;
      tmp_minus.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,tmp_minus_two,max_noofgridpoints);
      // calculate S_nl tmp
      tl=(yindex+1.0)*(yindex+1.0)/sqrt((2.0*yindex+1)*(2.0*yindex+3));
      fnl_nom=real(0.25*timestep*hamil.vecpot_z(real(time+timestep),me)*tl);
      for (xindex=0; xindex<max_noofgridpoints; xindex++)
	{
	  index=g.index(xindex,yindex,0);
	  index_lp=g.index(xindex,yindex+1,0);
	  fnl=fnl_nom/g.r(xindex);
	  fnlfnl=fnl*fnl;
	  ul=1.0-2.0*fnl-fnlfnl;
	  ur=-1.0-2.0*fnl+fnlfnl;
	  factor=1.0/(2.0*(1.0+fnlfnl));
	  start[index]=(ul*tmp_plus[xindex]+ur*tmp_minus[xindex])*factor;
	  start[index_lp]=(-ur*tmp_plus[xindex]+ul*tmp_minus[xindex])*factor;
	};
    };
  
}



void wavefunction::do_cn_step_xy(complex timestep, double time, grid g, hamop hamil, int me, int vecpotflag, long zindex, const fluid &wf_one, const wavefunction &wf_two)
{
  long xindex, yindex;
  long index, index_xp, index_xm, index_ym, index_yp;
  double x, y, z;
  wavefunction rhs(g.ngps_x()*g.ngps_y()*g.ngps_z());
  double const_diag_part, field, halfvecpotxvecpotx, halfvecpotyvecpoty;
  complex vpx, vpy;
  long max_gpsx_gpsy;

  if (g.ngps_x()>g.ngps_y())
    {
      max_gpsx_gpsy=g.ngps_x();
    }
  else
    {
      max_gpsx_gpsy=g.ngps_y();
    };

  wavefunction a(max_gpsx_gpsy);
  wavefunction b(max_gpsx_gpsy);
  wavefunction c(max_gpsx_gpsy);
  wavefunction wf_of_interest(max_gpsx_gpsy);
  wavefunction solveresult(max_gpsx_gpsy);

  z=g.z(zindex);

  // *** calculate the rhs which (y stuff applied to the wf) ***
  const_diag_part=1.0/(g.delt_y()*g.delt_y());
  vpy=-0.5/(g.delt_y()*g.delt_y())-0.5*complex(0.0,1.0)*hamil.vecpot_y(time,me)/g.delt_y();

  halfvecpotxvecpotx=0.5*hamil.vecpot_x(time,me)*hamil.vecpot_x(time,me);
  halfvecpotyvecpoty=0.5*hamil.vecpot_y(time,me)*hamil.vecpot_y(time,me);

  field=hamil.field(time,me);

  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      x=g.x(xindex);
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{ 
	  y=g.y(yindex);
	  a[yindex]=-vpy*0.5*complex(0.0,1.0)*timestep;
	  c[yindex]=-conj(vpy)*0.5*complex(0.0,1.0)*timestep;
	  b[yindex]=1.0-(const_diag_part
			 +hamil.scalarpoty(x,y,z,time,me)
			 +0.5*hamil.dftpot(g,x,y,z,time,me,wf_one,wf_two)
			 - complex(0.0,1.0)*0.5*hamil.imagpot(xindex,yindex,zindex,time,g)
			 +y*field
			 +0.5*halfvecpotyvecpoty )
	    *0.5*complex(0.0,1.0)*timestep;
	};
  

      yindex=0;
      index=g.index(xindex,yindex,zindex);
      index_yp=g.index(xindex,yindex+1,zindex);

      rhs[index]=a[yindex]*start[index_yp]+b[yindex]*start[index];

      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_yp=g.index(xindex,yindex+1,zindex);
	  index_ym=g.index(xindex,yindex-1,zindex);

	  rhs[index]=a[yindex]*start[index_yp]+b[yindex]*start[index]+c[yindex]*start[index_ym];
	};

      yindex=g.ngps_y()-1;
      index=g.index(xindex,yindex,zindex);
      index_ym=g.index(xindex,yindex-1,zindex);

      rhs[index]=c[yindex]*start[index_ym]+b[yindex]*start[index];

    };

  // *** do the CN-step (inversion with respect to x)
  const_diag_part=1.0/(g.delt_x()*g.delt_x());
  vpx=-0.5/(g.delt_x()*g.delt_x())-0.5*complex(0.0,1.0)*hamil.vecpot_x(time,me)/g.delt_x();
  
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    { 
      y=g.y(yindex);
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  x=g.x(xindex);
	  index=g.index(xindex,yindex,zindex);
	  wf_of_interest[xindex]=rhs[index];
	  a[xindex]=vpx*0.5*complex(0.0,1.0)*timestep;
	  c[xindex]=conj(vpx)*0.5*complex(0.0,1.0)*timestep;
	  b[xindex]=1.0+(const_diag_part
			 +hamil.scalarpotx(x,y,z,time,me)
			 +0.5*hamil.dftpot(g,x,y,z,time,me,wf_one,wf_two)
			 - complex(0.0,1.0)*0.5*hamil.imagpot(xindex,yindex,zindex,time,g)
			 +x*field
			 +halfvecpotxvecpotx )
	    *0.5*complex(0.0,1.0)*timestep;
	};
      solveresult.solve(a,b,c,wf_of_interest,g.ngps_x());
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
          start[index]=solveresult[xindex];
	};      
    };







  // *** calculate the rhs which (x stuff applied to the wf) ***
  const_diag_part=1.0/(g.delt_x()*g.delt_x());
  vpx=-0.5/(g.delt_x()*g.delt_x())-0.5*complex(0.0,1.0)*hamil.vecpot_x(time,me)/g.delt_x();

  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      y=g.y(yindex);
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  x=g.x(xindex);
	  a[xindex]=-vpx*0.5*complex(0.0,1.0)*timestep;
	  c[xindex]=-conj(vpx)*0.5*complex(0.0,1.0)*timestep;
	  b[xindex]=1.0-(const_diag_part
			 +hamil.scalarpotx(x,y,z,time,me)
			 +0.5*hamil.dftpot(g,x,y,z,time,me,wf_one,wf_two)
			 - complex(0.0,1.0)*0.5*hamil.imagpot(xindex,yindex,zindex,time,g)
			 +x*field
			 +halfvecpotxvecpotx )
	    *0.5*complex(0.0,1.0)*timestep;
	};
      
      
      xindex=0;
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
      
      rhs[index]=a[xindex]*start[index_xp]+b[xindex]*start[index];
      
      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  
	  rhs[index]=a[xindex]*start[index_xp]+b[xindex]*start[index]+c[xindex]*start[index_xm];
	};

      xindex=g.ngps_x()-1;
      index=g.index(xindex,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);

      rhs[index]=c[xindex]*start[index_xm]+b[xindex]*start[index];

    };


  // *** do the CN-step (inversion with respect to y)
  const_diag_part=1.0/(g.delt_y()*g.delt_y());
  vpy=-0.5/(g.delt_y()*g.delt_y())-0.5*complex(0.0,1.0)*hamil.vecpot_y(time,me)/g.delt_y();
  
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    { 
      x=g.x(xindex);
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  y=g.y(yindex); 
	  index=g.index(xindex,yindex,zindex);
	  wf_of_interest[yindex]=rhs[index];
	  a[yindex]=vpy*0.5*complex(0.0,1.0)*timestep;
	  c[yindex]=conj(vpy)*0.5*complex(0.0,1.0)*timestep;
	  b[yindex]=1.0+(const_diag_part
			 +hamil.scalarpoty(x,y,z,time,me)
			 +0.5*hamil.dftpot(g,x,y,z,time,me,wf_one,wf_two)
			 - complex(0.0,1.0)*0.5*hamil.imagpot(xindex,yindex,zindex,time,g)
			 +y*field
			 +halfvecpotyvecpoty )
	    *0.5*complex(0.0,1.0)*timestep;
	};
      solveresult.solve(a,b,c,wf_of_interest,g.ngps_y());
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
          start[index]=solveresult[yindex];
	};      
    };



  
}



void wavefunction::dump_to_file(grid g,FILE* os, int stepwidth)
{

  long xindex,yindex,zindex,i,counter, counter_ii, counter_iii;
  double u, r, rho, z, legpolnew, legpolold, legpololder;
  complex summingres;
  


  counter=0;
  counter_ii=0;
  counter_iii=0;

  switch (g.dimens())
    {
    case 2 :  case 6 : case 15 : case 16 : case 99 :
      counter=0;
      for (xindex=0; xindex<g.ngps_x(); xindex+=stepwidth)
	{ 
	  counter++;
	  counter_ii=0;
	  for (yindex=0; yindex<g.ngps_y(); yindex+=stepwidth)
	    {
	      counter_ii++;
	      i=g.index(xindex,yindex,0);
	      fprintf(os,"%.14le %.14le\n",real(start[i]),imag(start[i]));
	    };
	};
      break;
    };
}


void wavefunction::dump_dens1d_to_file(grid g, FILE* os, int stepwidth)
{

  wavefunction density(g.ngps_x());
  long index, xindex, yindex, zindex;
  long index_xp, index_xm;


  switch (g.dimens())
  {
  case 16:
  zindex=0;
  density.nullify();


  //density
  for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      density[xindex]=density[xindex]+2.0*real(conj(start[index])*start[index])*g.delt_y();
	    };
	};


  //write to file
  for (xindex=0; xindex<g.ngps_x(); xindex+=stepwidth)
    { 
      fprintf(os,"%.14le\n",real(density[xindex]));
    };
  break;

  };


};











void wavefunction::dump_wf1d_to_file(grid g, FILE* os, int stepwidth)
{

  wavefunction rhs(g.ngps_x()*g.ngps_y()*g.ngps_z());
  wavefunction rhs_x(g.ngps_x());
  wavefunction result_x(g.ngps_x());
  wavefunction current(g.ngps_x());
  wavefunction density(g.ngps_x());
  wavefunction theta(g.ngps_x());
  wavefunction wf1d(g.ngps_x());
  long index, xindex, yindex, zindex;
  long index_xp, index_xm;
  long xprimeindex;
  double aaa, ccc, bbb;
  double b_upperleft, b_lowerright;
  complex imagi(0.0,1.0);
  double fcteps=1.0E-308;
  double lambda=sqrt(3.0)-2.0;
  double llambda=-sqrt(3.0)+2.0;
  double thetaoffset;

  switch (g.dimens())
  {
  case 16:
  zindex=0;
  current.nullify();
  density.nullify();
  theta.nullify();


  //density
  for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      density[xindex]=density[xindex]+2.0*real(conj(start[index])*start[index])*g.delt_y();
	    };
	};


  //rhs
  aaa=1.0/(2.0*g.delt_x());
  bbb=0.0;
  ccc=-aaa;


  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      xindex=0;
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
  
      rhs[index]=aaa*start[index_xp]+lambda/(2.0*g.delt_x())*start[index];

      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  
	  rhs[index]=aaa*start[index_xp]+bbb*start[index]+ccc*start[index_xm];
	};
      
      xindex=g.ngps_x()-1;
      index=g.index(xindex,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);
      
      rhs[index]=llambda/(2.0*g.delt_x())*start[index]+ccc*start[index_xm];
    }

  //lhs
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhs_x[xindex]=rhs[index];
	};

      aaa=OOS;
      ccc=OOS;
      bbb=OOS*4.0;

      b_upperleft=(4.0+lambda)/6.0;
      b_lowerright=(4.0+lambda)/6.0;

      result_x.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,rhs_x,g.ngps_x());

      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhs[index]=result_x[xindex];
	};
    };
  
  //current
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  current[xindex]=current[xindex]-2.0*real(conj(start[index])*rhs[index]*imagi*g.delt_y());
	};
    };
 
  //theta
  thetaoffset=0.0;
  xindex=g.offs_x();
  for (xprimeindex=0; xprimeindex<=xindex; xprimeindex++)
    {
      thetaoffset=thetaoffset+real(current[xprimeindex])/(real(density[xprimeindex])+fcteps)*g.delt_x();
    };

  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      for (xprimeindex=0; xprimeindex<=xindex; xprimeindex++)
	{
	  theta[xindex]=theta[xindex]+real(current[xprimeindex])/(real(density[xprimeindex])+fcteps)*g.delt_x();
	};
      theta[xindex]=theta[xindex]-thetaoffset;
    };


  //wf1d
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      wf1d[xindex]=sqrt(0.5*real(density[xindex]))*exp(imagi*real(theta[xindex]));
    };  
      
  //write to file
  for (xindex=0; xindex<g.ngps_x(); xindex+=stepwidth)
    { 
      fprintf(os,"%.14le %.14le\n",real(wf1d[xindex]),imag(wf1d[xindex]));
    };
  break;


};
}


wavefunction wavefunction::arccos(grid g, double deltat)
{
  wavefunction rhs(g.ngps_x());
  long xindex;
  complex imagi=complex(0.0,1.0);

  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      //      rhs[xindex]=-(1.0/deltat)*asin(imag(start[xindex]));
      rhs[xindex]=imagi*log(start[xindex])/deltat;
    };


  return rhs;


};

wavefunction wavefunction::dens1d(grid g)
{

  wavefunction rhs(g.ngps_x()*g.ngps_y()*g.ngps_z());
  wavefunction rhs_x(g.ngps_x());
  wavefunction result_x(g.ngps_x());
  wavefunction current(g.ngps_x());
  wavefunction density(g.ngps_x());
  wavefunction theta(g.ngps_x());
  wavefunction wf1d(g.ngps_x());
  long index, xindex, yindex, zindex;
  long index_xp, index_xm;
  long xprimeindex;
  double aaa, ccc, bbb;
  double b_upperleft, b_lowerright;
  complex imagi(0.0,1.0);
  double fcteps=1.0E-308;
  double lambda=sqrt(3.0)-2.0;
  double llambda=-sqrt(3.0)+2.0;
  double thetaoffset;

  switch (g.dimens())
  {
  case 16:
  zindex=0;
  current.nullify();
  density.nullify();
  theta.nullify();


  //density
  for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      density[xindex]=density[xindex]+2.0*real(conj(start[index])*start[index])*g.delt_y();
	    };
	};


  //rhs
  aaa=1.0/(2.0*g.delt_x());
  bbb=0.0;
  ccc=-aaa;


  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      xindex=0;
      index=g.index(xindex,yindex,zindex);
      index_xp=g.index(xindex+1,yindex,zindex);
  
      rhs[index]=aaa*start[index_xp]+lambda/(2.0*g.delt_x())*start[index];

      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
	  
	  rhs[index]=aaa*start[index_xp]+bbb*start[index]+ccc*start[index_xm];
	};
      
      xindex=g.ngps_x()-1;
      index=g.index(xindex,yindex,zindex);
      index_xm=g.index(xindex-1,yindex,zindex);
      
      rhs[index]=llambda/(2.0*g.delt_x())*start[index]+ccc*start[index_xm];
    }

  //lhs
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhs_x[xindex]=rhs[index];
	};

      aaa=OOS;
      ccc=OOS;
      bbb=OOS*4.0;

      b_upperleft=(4.0+lambda)/6.0;
      b_lowerright=(4.0+lambda)/6.0;

      result_x.solve_toep(aaa,bbb,b_upperleft,b_lowerright,ccc,rhs_x,g.ngps_x());

      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  rhs[index]=result_x[xindex];
	};
    };
  
  //current
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  index=g.index(xindex,yindex,zindex);
	  current[xindex]=current[xindex]-2.0*real(conj(start[index])*rhs[index]*imagi*g.delt_y());
	};
    };
 
  //theta
  thetaoffset=0.0;
  xindex=g.offs_x();
  for (xprimeindex=0; xprimeindex<=xindex; xprimeindex++)
    {
      thetaoffset=thetaoffset+real(current[xprimeindex])/(real(density[xprimeindex])+fcteps)*g.delt_x();
    };

  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      for (xprimeindex=0; xprimeindex<=xindex; xprimeindex++)
	{
	  theta[xindex]=theta[xindex]+real(current[xprimeindex])/(real(density[xprimeindex])+fcteps)*g.delt_x();
	};
      theta[xindex]=theta[xindex]-thetaoffset;
    };


  //wf1d
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      wf1d[xindex]=sqrt(0.5*real(density[xindex]))*exp(imagi*real(theta[xindex]));
    };  
      
  return wf1d;
  break;


};
}






void  wavefunction::calculate_fixed_potential_array(grid g, hamop hamil, double time, int me)
{  
  double x,y,z;			
  long index,xindex,yindex,zindex;
  complex imagi(0.0,1.0);

  switch (g.dimens())
    {
    case 15 :
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  x=g.x(xindex);
	  yindex=0;
	  zindex=0;
	  index=g.index(xindex,yindex,zindex);
	  start[index]=hamil.scalarpot(x,y,z,time,me)-imagi*hamil.imagpot(xindex,yindex,zindex,time,g);
	};
      break;
    case 16 :
      zindex=0;
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  x=g.x(xindex);
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    { 
	      y=g.y(yindex);
	      index=g.index(xindex,yindex,zindex);
	      start[index]=hamil.scalarpot(x,y,z,time,me)-imagi*hamil.imagpot(xindex,yindex,zindex,time,g);
	    };
	};
      break;
      
    case 14 : case 24 :
      
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  x=g.r(xindex);
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    { 
	      index=g.index(xindex,yindex,0);
	      start[index]=hamil.scalarpot(x,yindex,0.0,0.0,me)+0.5*yindex*(yindex+1)/(x*x)
		-imagi*hamil.imagpot(xindex,yindex,0,time,g);
	    };
	};
      break;

    };

}

complex wavefunction::energy(double time, grid g, hamop hamil,int me, const double masses[], const wavefunction &staticpot_x, const wavefunction &staticpot_y, const wavefunction &staticpot_xy, double charge)
{
  double Eff_Charge = charge;

  complex result(0.0,0.0);
  complex result_x(0.0,0.0);
  complex result_y(0.0,0.0);
  complex result_xy(0.0,0.0);
  double const_diag_part;
  long xindex, yindex, zindex;
  double x,y,z,r,field,halfvecpotvecpot;
  complex vpx, vpy, vpz;
  long index, index_xp, index_xm, index_yp, index_ym, index_zp, index_zm;
  double b_upperleft, M_two_upperleft, Delta_two_upperleft;
  wavefunction tmp(g.ngps_x());	
  wavefunction tmp_two(g.ngps_x());	
  wavefunction rhsone_y(g.ngps_y());	
  wavefunction rhstwo_y(g.ngps_y());	
  wavefunction rhsone(g.ngps_x()*g.ngps_y()*g.ngps_z());
  wavefunction staticpot(g.ngps_x()*g.ngps_y()*g.ngps_z());
  double aaaa, bbbb, cccc;
  double aaa, bbb, ccc, hh;

  switch (g.dimens())
    {
    case 3 :
      const_diag_part=1.0/(g.delt_x()*g.delt_x())+1.0/(g.delt_y()*g.delt_y())
	+1.0/(g.delt_z()*g.delt_z());

      vpx=-0.5/(g.delt_x()*g.delt_x())
	-0.5*complex(0.0,1.0)*hamil.vecpot_x(time,me)/g.delt_x();
      vpy=-0.5/(g.delt_y()*g.delt_y())
	-0.5*complex(0.0,1.0)*hamil.vecpot_y(time,me)/g.delt_y();
      vpz=-0.5/(g.delt_z()*g.delt_z())
	-0.5*complex(0.0,1.0)*hamil.vecpot_z(time,me)/g.delt_z();
      
      halfvecpotvecpot=0.5*(
			       hamil.vecpot_z(time,me)*hamil.vecpot_z(time,me)
			       +hamil.vecpot_y(time,me)*hamil.vecpot_y(time,me)
			       +hamil.vecpot_x(time,me)*hamil.vecpot_x(time,me)
			       );

      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  x=g.x(xindex);
	  for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	    {
	      y=g.y(yindex);
	      for (zindex=1; zindex<g.ngps_z()-1; zindex++)
		{
		  z=g.z(zindex);
		  
		  index=g.index(xindex,yindex,zindex);
		  index_xp=g.index(xindex+1,yindex,zindex);
		  index_xm=g.index(xindex-1,yindex,zindex);
		  index_yp=g.index(xindex,yindex+1,zindex);
		  index_ym=g.index(xindex,yindex-1,zindex);
		  index_zp=g.index(xindex,yindex,zindex+1);
		  index_zm=g.index(xindex,yindex,zindex-1);
		  
		  result=result+conj(start[index])*
		    (  (const_diag_part
			+hamil.scalarpot(x,y,z,time,me)
			+(x+y+z)*hamil.field(time,me) 
			+ halfvecpotvecpot
			)*start[index]
		       +vpx*start[index_xp]+conj(vpx)*start[index_xm]
		       +vpy*start[index_yp]+conj(vpy)*start[index_ym]
		       +vpz*start[index_zp]+conj(vpz)*start[index_zm] )
		    *g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
      break;
    case 2 : case 6 :
      const_diag_part=1.0/(g.delt_x()*g.delt_x())+1.0/(g.delt_y()*g.delt_y());
      vpx=-0.5/(g.delt_x()*g.delt_x());
      vpy=-0.5/(g.delt_y()*g.delt_y());
      z=0.0;
      zindex=0;


      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  x=g.x(xindex);
	  for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	    {
	      y=g.y(yindex);
			  
	      index=g.index(xindex,yindex,zindex);
	      index_xp=g.index(xindex+1,yindex,zindex);
	      index_xm=g.index(xindex-1,yindex,zindex);
	      index_yp=g.index(xindex,yindex+1,zindex);
	      index_ym=g.index(xindex,yindex-1,zindex);
	      
	      result=result+conj(start[index])*
		(  (const_diag_part
		    +hamil.scalarpot(x,y,z,time,me)
		    )*start[index]
		   +vpx*start[index_xp]+conj(vpx)*start[index_xm]
		   +vpy*start[index_yp]+conj(vpy)*start[index_ym])
		*g.delt_x()*g.delt_y();
	    };
	};
    break;
    case 16 :
      // ------------------------ x-part ----------------
      zindex=0;
      hh=g.delt_x()*g.delt_x();
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  y=g.y(yindex);
	  
	  xindex=0;
	  x=g.x(xindex);
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);

	  rhsone[index]=(-2.0/hh-5.0/3.0*staticpot_x[xindex])*start[index] 
	    + (1.0/hh-1.0/6.0*staticpot_x[xindex+1])*start[index_xp];

	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    {
	      x=g.x(xindex);
	      index=g.index(xindex,yindex,zindex);
	      index_xp=g.index(xindex+1,yindex,zindex);
	      index_xm=g.index(xindex-1,yindex,zindex);

	      rhsone[index]=(1.0/hh-1.0/6.0*staticpot_x[xindex-1])*start[index_xm]
		+ (-2.0/hh-5.0/3.0*staticpot_x[xindex])*start[index] 
		+ (1.0/hh-1.0/6.0*staticpot_x[xindex+1])*start[index_xp];
	    };

	  xindex=g.ngps_x()-1;
	  x=g.x(xindex);
	  index=g.index(xindex,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);

	  rhsone[index]=(1.0/hh-1.0/6.0*staticpot_x[xindex-1])*start[index_xm]
	    + (-2.0/hh-5.0/3.0*staticpot_x[xindex])*start[index]; 

	};	  

      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      tmp[xindex]=rhsone[index];
	    };
	  
	  aaa=-OOS;
	  ccc=aaa;
	  bbb=-FOT;

	  tmp_two.solve_toep(aaa,bbb,bbb,bbb,ccc,tmp,g.ngps_x());
	  
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      rhsone[index]=tmp_two[xindex];
	    };
	};

      result_x=(*this)*rhsone;

      // ------------------------ y-part ----------------
      hh=g.delt_y()*g.delt_y();
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  x=g.x(xindex);
	  
	  yindex=0;
	  y=g.y(yindex);
	  index=g.index(xindex,yindex,zindex);
	  index_yp=g.index(xindex,yindex+1,zindex);

	  rhsone[index]=(-2.0/hh-5.0/3.0*staticpot_y[yindex])*start[index] 
	    + (1.0/hh-1.0/6.0*staticpot_y[yindex+1])*start[index_yp];

	  for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	    {
	      y=g.y(yindex);
	      index=g.index(xindex,yindex,zindex);
	      index_yp=g.index(xindex,yindex+1,zindex);
	      index_ym=g.index(xindex,yindex-1,zindex);

	      rhsone[index]=(1.0/hh-1.0/6.0*staticpot_y[yindex-1])*start[index_ym]
		+ (-2.0/hh-5.0/3.0*staticpot_y[yindex])*start[index] 
		+ (1.0/hh-1.0/6.0*staticpot_y[yindex+1])*start[index_yp];
	    };

	  yindex=g.ngps_y()-1;
	  y=g.y(yindex);
	  index=g.index(xindex,yindex,zindex);
	  index_ym=g.index(xindex,yindex-1,zindex);

	  rhsone[index]=(1.0/hh-1.0/6.0*staticpot_y[yindex-1])*start[index_ym]
	    + (-2.0/hh-5.0/3.0*staticpot_y[yindex])*start[index]; 

	};	  

      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      rhsone_y[yindex]=rhsone[index];
	    };
	  
	  aaa=-OOS;
	  ccc=aaa;
	  bbb=-FOT;

	  rhstwo_y.solve_toep(aaa,bbb,bbb,bbb,ccc,rhsone_y,g.ngps_y());
	  
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      rhsone[index]=rhstwo_y[yindex];
	    };
	};

      result_y=(*this)*rhsone;

      // ------------------------ xy-part ----------------

     for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,zindex);
	      rhsone[index]=staticpot_xy[index]*start[index];
	    };
	};

     result_xy=(*this)*rhsone;

      result=(result_x+result_y+result_xy)*g.delt_x()*g.delt_y();

    break;

    case 99 :
      const_diag_part=1.0/(g.delt_x()*g.delt_x()*masses[0])+1.0/(g.delt_y()*g.delt_y()*masses[1]);
      z=0.0;
      zindex=0;
      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  x=g.x(xindex);
	  for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	    {
	      y=g.y(yindex);
		  
	      vpx=-0.5/(g.delt_x()*g.delt_x()*masses[0])
		-0.5*complex(0.0,1.0)*y*hamil.vecpot_x(time,me)/(g.delt_x()*masses[0]);
	      vpy=-0.5/(g.delt_y()*g.delt_y()*masses[1]);
		  
		  index=g.index(xindex,yindex,zindex);
		  index_xp=g.index(xindex+1,yindex,zindex);
		  index_xm=g.index(xindex-1,yindex,zindex);
		  index_yp=g.index(xindex,yindex+1,zindex);
		  index_ym=g.index(xindex,yindex-1,zindex);
		  
		  result=result+conj(start[index])*
		    (  (const_diag_part
			+hamil.scalarpot(x,y,z,time,me)
			+0.5*y*y*hamil.vecpot_x(time,me)*hamil.vecpot_x(time,me)
			)*start[index]
		       +vpx*start[index_xp]+conj(vpx)*start[index_xm]
		       +vpy*start[index_yp]+conj(vpy)*start[index_ym])
		    *g.delt_x()*g.delt_y();
	    };
	};
    break;
    case 1 :  case 5 : case 15 :  
      const_diag_part=1.0/(g.delt_x()*g.delt_x());
      y=0.0;
      z=0.0;
      yindex=0;
      zindex=0;
      halfvecpotvecpot=0.5*hamil.vecpot_x(time,me)*hamil.vecpot_x(time,me);
      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  x=g.x(xindex);
		  
	  vpx=-0.5/(g.delt_x()*g.delt_x())
	    -0.5*complex(0.0,1.0)*hamil.vecpot_x(time,me)/g.delt_x();
		  
	  index=g.index(xindex,yindex,zindex);
	  index_xp=g.index(xindex+1,yindex,zindex);
	  index_xm=g.index(xindex-1,yindex,zindex);
		  

	  result=result+conj(start[index])*
	    (  
	     ( const_diag_part
	       +hamil.scalarpot(x,y,z,time,me)
	       //	       +x*hamil.field(time,me)
	       +halfvecpotvecpot
	       )*start[index]
	     +vpx*start[index_xp]+conj(vpx)*start[index_xm]
	       )
	    *g.delt_x();
	};
      break;
    case 4 : case 7 :
      const_diag_part=1.0/(g.delt_x()*g.delt_x());
      r=g.r(0);

      field=hamil.field(0.0,me);

      yindex=0;
      index=g.index(0,0,0);
      index_xp=g.index(1,0,0);
      index_yp=g.index(0,1,0);

      result=result+conj(start[index])/(2.0*yindex+1.0)*
	(  (const_diag_part
	    +hamil.scalarpot(r,0.0,0.0,time,me)
	    +(yindex*(yindex+1.0))/(2.0*r*r) )*start[index]
	   -0.5/(g.delt_x()*g.delt_x())*start[index_xp] 
	   +field*(1.0*(yindex+1.0)/(2.0*yindex+3.0)*r*start[index_yp] ))
	    *g.delt_x()*12.56637061;


      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	{
	  index=g.index(0,yindex,0);
	  index_xp=g.index(1,yindex,0);
	  index_yp=g.index(0,yindex+1,0);
	  index_ym=g.index(0,yindex-1,0);
	  
	  result=result+conj(start[index])/(2.0*yindex+1.0)*
	    (  (const_diag_part
		+hamil.scalarpot(r,0.0,0.0,time,me)
		+(yindex*(yindex+1.0))/(2.0*r*r) )*start[index]
	       -0.5/(g.delt_x()*g.delt_x())*start[index_xp] 
	       +field*(1.0*yindex/(2.0*yindex-1.0)*r*start[index_ym]
						  +1.0*(yindex+1.0)/(2.0*yindex+3.0)*r*start[index_yp] )
	       )
	    *g.delt_x()*12.56637061;
	};

      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{ 
	  r=g.r(xindex);
	  yindex=0;
	  index=g.index(xindex,yindex,0);
	  index_xp=g.index(xindex+1,yindex,0);
	  index_xm=g.index(xindex-1,yindex,0);
	  index_yp=g.index(xindex,yindex+1,0);
		  
	  result=result+conj(start[index])/(2.0*yindex+1.0)*
	    (  (const_diag_part
		+hamil.scalarpot(r,0.0,0.0,time,me)
		+(yindex*(yindex+1.0))/(2.0*r*r) )*start[index]
	       -0.5/(g.delt_x()*g.delt_x())*start[index_xp]
	       -0.5/(g.delt_x()*g.delt_x())*start[index_xm] 
	       +field*(1.0*(yindex+1.0)/(2.0*yindex+3.0)*r*start[index_yp]))
	    *g.delt_x()*12.56637061;

	  for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      index_xp=g.index(xindex+1,yindex,0);
	      index_xm=g.index(xindex-1,yindex,0);
	      index_yp=g.index(xindex,yindex+1,0);
	      index_ym=g.index(xindex,yindex-1,0);
		  
	      result=result+conj(start[index])/(2.0*yindex+1.0)*
		(  (const_diag_part
		    +hamil.scalarpot(r,0.0,0.0,time,me)
		    +(yindex*(yindex+1.0))/(2.0*r*r) )*start[index]
		   -0.5/(g.delt_x()*g.delt_x())*start[index_xp]
		   -0.5/(g.delt_x()*g.delt_x())*start[index_xm] 
		   +field*(1.0*yindex/(2.0*yindex-1.0)*r*start[index_ym]
						  +1.0*(yindex+1.0)/(2.0*yindex+3.0)*r*start[index_yp] )
		   )		
		*g.delt_x()*12.56637061;
	    };
	};
      break;
    case 14 : case 24 :
      aaaa=1.0/(g.delt_x()*g.delt_x());
      cccc=aaaa;
      bbbb=-2.0/(g.delt_x()*g.delt_x());
      Delta_two_upperleft=-2.0/(g.delt_x()*g.delt_x())*(1.0-Eff_Charge*g.delt_x()/(12.0-10.0*Eff_Charge*g.delt_x()));
      M_two_upperleft=-2.0*(1.0+g.delt_x()*g.delt_x()/12.0*Delta_two_upperleft);
  
      // apply (Delta_2 + M_2 V_eff,l) to Phi
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{ 
	  xindex=0;
	  index=g.index(xindex,yindex,0);
	  index_xp=g.index(xindex+1,yindex,0);
	  if (yindex==0)
	    {
	      rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		+(M_two_upperleft*(staticpot[index])+Delta_two_upperleft)*start[index];
	    }
	  else
	    {
	      rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		+(bbbb-FOT*(staticpot[index]))*start[index];
	    };
	  
	  for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	    { 
	      index=g.index(xindex,yindex,0);
	      index_xp=g.index(xindex+1,yindex,0);
	      index_xm=g.index(xindex-1,yindex,0);
	      rhsone[index]=(aaaa-(staticpot[index_xp])/6.0)*start[index_xp]
		+(bbbb-FOT*(staticpot[index]))*start[index]
		+(cccc-(staticpot[index_xm])/6.0)*start[index_xm];
	    };
	  
	  
	  xindex=g.ngps_x()-1;
	  index=g.index(xindex,yindex,0);
	  index_xm=g.index(xindex-1,yindex,0);
	  rhsone[index]=(bbbb-FOT*(staticpot[index]))*start[index]
	    +(cccc-(staticpot[index_xm])/6.0)*start[index_xm];
	};
      
      // solve M_2^-1 times that stuff
      aaaa=-OOS;
      cccc=aaaa;
      bbbb=-FOT;
      
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{ 
	  xindex=0;
	  index=g.index(xindex,yindex,0);
	  index_xp=g.index(xindex+1,yindex,0);

	  if (yindex==0)
	    {
	      b_upperleft=M_two_upperleft;
	    }
	  else
	    {
	      b_upperleft=-FOT;
	    };

	  
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      index=g.index(xindex,yindex,0);
	      tmp_two[xindex]=rhsone[index];
	    };  
	  
	  tmp.solve_toep(aaaa,bbbb,b_upperleft,bbbb,cccc,tmp_two,g.ngps_x());
	  
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      index=g.index(xindex,yindex,0);
	      rhsone[index]=tmp[xindex];
	    };  
      
	};

      // and finally do the integration 
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      result=result+conj(start[index])*rhsone[index]*g.delt_x();
	    };
	};


      break;
    };

  return result;
}
  

void  wavefunction::setzeroboundary(grid g)
{
  long xindex, yindex, zindex;

  zindex=0;

  xindex=0;
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      start[g.index(xindex,yindex,zindex)]=complex(0.0,0.0);
    };
  xindex=g.ngps_x()-1;
  for (yindex=0; yindex<g.ngps_y(); yindex++)
    {
      start[g.index(xindex,yindex,zindex)]=complex(0.0,0.0);
    };

  yindex=0;
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      start[g.index(xindex,yindex,zindex)]=complex(0.0,0.0);
    };
  yindex=g.ngps_y()-1;
  for (xindex=0; xindex<g.ngps_x(); xindex++)
    {
      start[g.index(xindex,yindex,zindex)]=complex(0.0,0.0);
    };

}

double wavefunction::norm(grid g)
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index,indexm,indexmm;

  switch (g.dimens())
    {
    case 14 : case 24 :
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      result=result+conj(start[index])*start[index]*g.delt_x();
	    };
	};
      break;
    case 4 : case 7 :
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      result=result+conj(start[index])*start[index]*g.delt_x()/(2.0*yindex+1)*12.56637061;
	    };
	};
      break;
    default :
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  for (zindex=0; zindex<g.ngps_z(); zindex++)
	    {
	      for (yindex=0; yindex<g.ngps_y(); yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+conj(start[index])*start[index]*g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
    };

  return real(result);

}


double wavefunction::non_ionized(grid g, long box)
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index;

  switch (g.dimens())
    {
    case 3 : 
      for (xindex=g.offs_x()-box; xindex<g.offs_x()+box; xindex++)
	{ 
	  for (zindex=g.offs_z()-box; zindex<g.offs_z()+box; zindex++)
	    {
	      for (yindex=g.offs_y()-box; yindex<g.offs_y()+box; yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+conj(start[index])*start[index]
		    *g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
      break;
    case 2 : case 6 : case 99 : case 16 :
      for (xindex=g.offs_x()-box; xindex<g.offs_x()+box; xindex++)
	{ 
	  for (yindex=g.offs_y()-box; yindex<g.offs_y()+box; yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      
	      result=result+conj(start[index])*start[index]
		*g.delt_x()*g.delt_y();
	    };
	};
      break;
    case 1 : case 5 : case 15 :
      for (xindex=g.offs_x()-box; xindex<g.offs_x()+box; xindex++)
	{
	  result=result+conj(start[xindex])*start[xindex]
	    *g.delt_x();
	};
      break;
    case 14 : case 24 :
      for (xindex=0; xindex<box; xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      result=result+conj(start[index])*start[index]*g.delt_x();
	    };
	};
      break;
    case 4 : case 7 :
      for (xindex=0; xindex<box; xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      result=result+conj(start[index])*start[index]*g.delt_x()/(2.0*yindex+1)*12.56637061;
	    };
	};
      break;
    default : result=0.0;
    };

  return real(result);

}

double wavefunction::sing_ionized(grid g, long box)
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index;

  switch (g.dimens())
    {
    case 3 :
      for (xindex=g.offs_x()-box; xindex<g.offs_x()+box; xindex++)
	{ 
	  for (zindex=g.offs_z()-box; zindex<g.offs_z()+box; zindex++)
	    {
	      for (yindex=0; yindex<g.offs_y()-box; yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+conj(start[index])*start[index]
		    *g.delt_x()*g.delt_y()*g.delt_z();
		};
	      for (yindex=g.offs_y()+box; yindex<g.ngps_y(); yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+conj(start[index])*start[index]
		    *g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
      for (xindex=g.offs_x()-box; xindex<g.offs_x()+box; xindex++)
	{ 
	  for (yindex=g.offs_y()-box; yindex<g.offs_y()+box; yindex++)
	    {
	      for (zindex=0; zindex<g.offs_z()-box; zindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+conj(start[index])*start[index]
		    *g.delt_x()*g.delt_y()*g.delt_z();
		};
	      for (zindex=g.offs_z()+box; zindex<g.ngps_z(); zindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+conj(start[index])*start[index]
		    *g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
      for (zindex=g.offs_z()-box; zindex<g.offs_z()+box; zindex++)
	{ 
	  for (yindex=g.offs_y()-box; yindex<g.offs_y()+box; yindex++)
	    {
	      for (xindex=0; xindex<g.offs_x()-box; xindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+conj(start[index])*start[index]
		    *g.delt_x()*g.delt_y()*g.delt_z();
		};
	      for (xindex=g.offs_x()+box; xindex<g.ngps_x(); xindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+conj(start[index])*start[index]
		    *g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
      break;
    case 2 : case 6 : case 16 : case 99 :
      for (xindex=g.offs_x()-box; xindex<g.offs_x()+box; xindex++)
	{ 
	  for (yindex=0; yindex<g.offs_y()-box; yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      
	      result=result+conj(start[index])*start[index]
		*g.delt_x()*g.delt_y();
	    };
	  for (yindex=g.offs_y()+box; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      
	      result=result+conj(start[index])*start[index]
		*g.delt_x()*g.delt_y();
	    };
	};
      for (yindex=g.offs_y()-box; yindex<g.offs_y()+box; yindex++)
	{
	  for (xindex=0; xindex<g.offs_x()-box; xindex++)
	    {
	      index=g.index(xindex,yindex,0);
		  
	      result=result+conj(start[index])*start[index]
		*g.delt_x()*g.delt_y();
	    };
	  for (xindex=g.offs_x()+box; xindex<g.ngps_x(); xindex++)
	    {
	      index=g.index(xindex,yindex,0);
		  
	      result=result+conj(start[index])*start[index]
		*g.delt_x()*g.delt_y();
	    };
	};
      break;
    case 1 : case 5 : case 15 :
      for (xindex=0; xindex<g.offs_x()-box; xindex++)
	{
	  result=result+conj(start[xindex])*start[xindex]
	    *g.delt_x();
	};
      for (xindex=g.offs_x()+box; xindex<g.ngps_x(); xindex++)
	{
	  result=result+conj(start[xindex])*start[xindex]
	    *g.delt_x();
	};     
      break;
    case 4 : case 7 : 
      for (xindex=box; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      result=result+conj(start[index])*start[index]*g.delt_x()/(2.0*yindex+1)*12.56637061;
	    };
	};
      break;
    case 14 : case 24 :
      for (xindex=box; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      result=result+conj(start[index])*start[index]*g.delt_x();
	    };
	};
      break;
    default : result=0.0;
    };


  return real(result);

}


double wavefunction::doub_ionized(grid g, long box)
{
  complex  result=(0.0,0.0);
  long xindex, yindex, zindex;
  long index;

switch (g.dimens())
    {
    case 3:
     for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  for (yindex=g.offs_y()-box; yindex<g.offs_y()+box; yindex++)
	    {
	      for (zindex=0; zindex<g.ngps_z(); zindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+conj(start[index])*start[index]*g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{ 
	  for (zindex=g.offs_z()-box; zindex<g.offs_z()+box; zindex++)
	    {
	      for (xindex=0; xindex<g.ngps_x(); xindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+conj(start[index])*start[index]*g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
      for (zindex=0; zindex<g.ngps_z(); zindex++)
	{ 
	  for (xindex=g.offs_x()-box; xindex<g.offs_x()+box; xindex++)
	    {
	      for (yindex=0; yindex<g.ngps_y(); yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+conj(start[index])*start[index]*g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
      result=result-3.0*non_ionized(g,box)-2.0*sing_ionized(g,box); 
    break;
    
    case 16 :
      for (xindex=0; xindex<g.offs_x()-box; xindex++)
	{ 
	  for (yindex=0; yindex<g.offs_y()-box; yindex++)
	    {
	      index=g.index(xindex,yindex,0);
  	      
	      result=result+conj(start[index])*start[index]
		*g.delt_x()*g.delt_y();
	    };
	  for (yindex=g.offs_y()+box; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      
	      result=result+conj(start[index])*start[index]
		*g.delt_x()*g.delt_y();
	    };
	};
      for (xindex=g.offs_x()+box; xindex<g.ngps_x(); xindex++)
	{
	  for (yindex=0; yindex<g.offs_y()-box; yindex++)
	    {
	      index=g.index(xindex,yindex,0);
		  
	      result=result+conj(start[index])*start[index]
		*g.delt_x()*g.delt_y();
	    };
	  for (yindex=g.offs_y()+box; yindex<g.ngps_y(); yindex++)
	    {
	      index=g.index(xindex,yindex,0);
		  
	      result=result+conj(start[index])*start[index]
		*g.delt_x()*g.delt_y();
	    };
	};
  
    break;

    default: result=0.0;
    };

  return real(result);

}

double wavefunction::expect_x(grid g)
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index;
  double x;

  switch (g.dimens())
    {
    case 4 : case 7 : case 14 :  case 24 :
      break;
    default :      
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{ 
	  x=g.x(xindex);
	  for (zindex=0; zindex<g.ngps_z(); zindex++)
	    {
	      for (yindex=0; yindex<g.ngps_y(); yindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+x*conj(start[index])*start[index]*g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
    };

  return real(result);

}




double wavefunction::expect_y(grid g)
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index;
  double y;

  switch (g.dimens())
    {
    case 4 : case 7 : case 14 : case 24 :
      break;
    default :
      for (yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  y=g.y(yindex);
	  for (xindex=0; xindex<g.ngps_x(); xindex++)
	    { 
	      for (zindex=0; zindex<g.ngps_z(); zindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+y*conj(start[index])*start[index]*g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
    };

  return real(result);

}


double wavefunction::expect_dipole(grid g)
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index;
  double x,y;

  switch (g.dimens())
    {
    case 4 : case 7 : case 14 : case 24 :
      break;
    default :
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  x=g.x(xindex);
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    { 
	      y=g.y(yindex); 
	      for (zindex=0; zindex<g.ngps_z(); zindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  result=result+(x+y)*conj(start[index])*start[index]*g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
    };

  return real(result);

}



complex wavefunction::expect_transition_dipole(grid g, const wavefunction &a)
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index;
  double x,y;

  switch (g.dimens())
    {
    case 4 : case 7 : case 14 : case 24 :
      break;
    default :
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  x=g.x(xindex);
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    { 
	      y=g.y(yindex); 
	      for (zindex=0; zindex<g.ngps_z(); zindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  result=result+(x+y)*conj(start[index])*a[index]*g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
    };

  return result;

}



double wavefunction::expect_dipole_acceleration(grid g, double eps)
{
  complex result(0.0,0.0);
  long xindex, yindex, zindex;
  long index;
  double x, y;

  switch (g.dimens())
    {
    case 16:
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  x=g.x(xindex);
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	     y=g.y(yindex); 
	     for (zindex=0; zindex<g.ngps_z(); zindex++)
		{
		  index=g.index(xindex,yindex,zindex);
		  
		  result=result+(2*x/pow(x*x+eps,3/2)+2*y/pow(y*y+eps,3/2))*conj(start[index])*start[index]*g.delt_x()*g.delt_y()*g.delt_z();
		};
	    };
	};
    break;  
    default: break;
    };

  return real(result);

}



void wavefunction::solve(const wavefunction &a, const wavefunction &b, 
			 const wavefunction &c, const wavefunction &psi, long dimens)
{

  wavefunction e(dimens), f(dimens);
  long i;


  e[0]=-b[0]/a[0];
  f[0]=psi[0]/a[0];
  for (i=1; i<dimens; i++)
    {
      e[i]=(-c[i]/e[i-1] - b[i])/a[i];
      f[i]=(psi[i] + c[i]/e[i-1]*f[i-1])/a[i];
    };

  start[dimens-1]=-f[dimens-1]/e[dimens-1];
  for (i=dimens-2; i>=0; i--)
    {
      start[i]=(start[i+1]-f[i])/e[i];
    };
  

}

// this is for pseudo-Toeplitz structured matrices, i.e., a, c is always the same,
// b also, apart from the corner elements. Thus b_upperleft and b_lowerright are extra arguments.
// version with a, b, c, double !!! 
void wavefunction::solve_toep(double a, double b, double b_upperleft, double b_lowerright, 
			 double c, const wavefunction &psi, long dimens)
{

  wavefunction f(dimens);
  fluid e(dimens);
  long i;
  double covera, bovera;

  covera=c/a;
  bovera=b/a;


  e[0]=-b_upperleft/a;
  f[0]=psi[0]/a;
  for (i=1; i<dimens-1; i++)
    {
      e[i]=-covera/e[i-1] - bovera;
      f[i]=psi[i]/a + covera/e[i-1]*f[i-1];
    };
  e[dimens-1]=-covera/e[dimens-2] - b_lowerright/a;
  f[dimens-1]=psi[dimens-1]/a + covera/e[dimens-2]*f[dimens-2];


  start[dimens-1]=-f[dimens-1]/e[dimens-1];
  for (i=dimens-2; i>=0; i--)
    {
      start[i]=(start[i+1]-f[i])/e[i];
    };
  

}


void wavefunction::apply(const wavefunction &a, const wavefunction &b, 
			 const wavefunction &c, const wavefunction &psi)
{

  long dimens=psi.wf_size();
  long i;


  start[0]=a[0]*psi[1]+b[0]*psi[0];

  for (i=1; i<dimens-1; i++)
    {
      start[i]=a[i]*psi[i+1]+b[i]*psi[i]+c[i]*psi[i-1];
    };
  
  start[dimens-1]=b[dimens-1]*psi[dimens-1]+c[dimens-1]*psi[dimens-2];


}



void  wavefunction::calculate_fixed_potential_array_x(grid g, hamop hamil, double time, int me)
{  
  double x,y,z;			
  long xindex,yindex,zindex;
  complex imagi(0.0,1.0);

  switch (g.dimens())
    {
    case 16:
      yindex=0;
      zindex=0;
      y=g.y(yindex);
      z=g.z(zindex);
      for(xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  x=g.x(xindex);
	  start[xindex]=hamil.scalarpotx(x,y,z,time,me)-imagi*hamil.imagpotx(xindex,yindex,zindex,time,g);
	};
      break;
    }
}

void  wavefunction::calculate_fixed_potential_array_y(grid g, hamop hamil, double time, int me)
{  
  double x,y,z;			
  long xindex,yindex,zindex;
  complex imagi(0.0,1.0);

  switch (g.dimens())
    {
    case 16:
      xindex=0;
      yindex=0;
      x=g.x(xindex);
      z=g.z(zindex);
      for(yindex=0; yindex<g.ngps_y(); yindex++)
	{
	  y=g.y(yindex);
	  start[yindex]=hamil.scalarpoty(x,y,z,time,me)-imagi*hamil.imagpoty(xindex,yindex,zindex,time,g);
	};
      break;
    }
}


void wavefunction::calculate_fixed_potential_array_xy(grid g, hamop hamil, double time, int me)
{  
  double x,y,z;			
  long index,xindex,yindex,zindex;

  switch (g.dimens())
    {
    case 16 :
      zindex=0;
      z=g.z(zindex);
      for (xindex=0; xindex<g.ngps_x(); xindex++)
	{
	  x=g.x(xindex);
	  for (yindex=0; yindex<g.ngps_y(); yindex++)
	    {
	      y=g.y(yindex);
	      index=g.index(xindex,yindex,zindex);
	      start[index]=hamil.interactionpotxy(x,y,z,time,me);
	    };
	};
      break;
    };
}


