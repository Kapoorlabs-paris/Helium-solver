#include<hamop.h>
#include<wavefunction.h>
#include<fluid.h>


hamop::hamop(grid g, 
	     double (*fpx)(double, int),
	     double (*fpy)(double, int),
	     double (*fpz)(double, int),
	     double (*fpsx)(double, double, double, double, int),
	     double (*fpsy)(double, double, double, double, int),
	     double (*fpsz)(double, double, double, double, int),
	     double (*fpixy)(double, double, double, double, int),
	     double (*fpimx)(long, long, long, double, grid),
	     double (*fpimy)(long, long, long, double, grid), 
	     double (*fpf)(double, int),
	     double (*dft)(grid, double, double, double, double, int, const fluid &, const wavefunction &)
	    )
{
  delta_x=g.delt_x();
  delta_y=g.delt_y();
  delta_z=g.delt_z();

  hamopvecpotx=fpx;
  hamopvecpoty=fpy;
  hamopvecpotz=fpz;
  hamopscalarpotx=fpsx;
  hamopscalarpoty=fpsy;
  hamopscalarpotz=fpsz;
  hamopinteractionpotxy=fpixy;
  hamopimagpotx=fpimx;
  hamopimagpoty=fpimy;
  hamopfield=fpf;
  hamopdftpot=dft;
}

double hamop::vecpot_x(double time, int me)
{
  return hamopvecpotx(time,me);

}  

double hamop::vecpot_y(double time, int me)
{
  return hamopvecpoty(time,me);

}  

double hamop::vecpot_z(double time, int me)
{
  return hamopvecpotz(time,me);

}  

double hamop::scalarpot(double x, double y, double z, double time, int me)
{
  return hamopscalarpotx(x,y,z,time,me)+hamopscalarpoty(x,y,z,time,me)+hamopscalarpotz(x,y,z,time,me);
}  

double hamop::scalarpotx(double x, double y, double z, double time, int me)
{
  return hamopscalarpotx(x,y,z,time,me);
}  
double hamop::scalarpoty(double x, double y, double z, double time, int me)
{
  return hamopscalarpoty(x,y,z,time,me);
}  
double hamop::scalarpotz(double x, double y, double z, double time, int me)
{
  return hamopscalarpotz(x,y,z,time,me);
}

double hamop::interactionpotxy(double x, double y, double z, double time, int me)
{
  return hamopinteractionpotxy(x,y,z,time,me);
}

double hamop::imagpot(long xindex, long yindex, long zindex, double time, grid g)
{
  return hamopimagpotx(xindex,yindex,zindex,time,g)+hamopimagpoty(xindex,yindex,zindex,time,g);
}  

double hamop::imagpotx(long xindex, long yindex, long zindex, double time, grid g)
{
  return hamopimagpotx(xindex,yindex,zindex,time,g);
}  

double hamop::imagpoty(long xindex, long yindex, long zindex, double time, grid g)
{
  return hamopimagpoty(xindex,yindex,zindex,time,g);
}  


double hamop::field(double time, int me)
{
  return hamopfield(time,me);

}  

double hamop::dftpot(grid g, double x, double y, double z, double time, int me, 
		     const fluid &wf_one, const wavefunction &wf_two)
{
  return hamopdftpot(g,x,y,z,time,me,wf_one,wf_two);
}






