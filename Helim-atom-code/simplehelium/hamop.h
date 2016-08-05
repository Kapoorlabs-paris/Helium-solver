#ifndef hamop_h
#define hamop_h hamop_h
#include<complex>
#include<iostream>
#include<grid.h>

using namespace std;

class wavefunction;
class fluid;

class hamop
{
 public:
  hamop(grid g, 
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
       );
  double vecpot_x(double time, int me);  
  double vecpot_y(double time, int me);  
  double vecpot_z(double time, int me);
  double scalarpot(double x, double y, double z, double time, int me);
  double scalarpotx(double x, double y, double z, double time, int me);
  double scalarpoty(double x, double y, double z, double time, int me);
  double scalarpotz(double x, double y, double z, double time, int me);
  double interactionpotxy(double x, double y, double z, double time, int me);
  double imagpot(long xindex, long yindex, long zindex, double time, grid g);
  double imagpotx(long xindex, long yindex, long zindex, double time, grid g);
  double imagpoty(long xindex, long yindex, long zindex, double time, grid g);
  double field(double time, int me);
  double dftpot(grid g, double x, double y, double z, double time, int me, 
		const fluid &wf_one, const wavefunction &wf_two);

 private:
  double delta_x, delta_y, delta_z;
  double (*hamopvecpotx)(double, int);
  double (*hamopvecpoty)(double, int);
  double (*hamopvecpotz)(double, int);
  double (*hamopscalarpotx)(double, double, double, double, int);
  double (*hamopscalarpoty)(double, double, double, double, int);
  double (*hamopscalarpotz)(double, double, double, double, int);
  double (*hamopinteractionpotxy)(double, double, double, double, int);
  double (*hamopimagpotx)(long, long, long, double, grid);
  double (*hamopimagpoty)(long, long, long, double, grid);
  double (*hamopfield)(double, int);
  double (*hamopdftpot)(grid, double, double, double, double, int, const fluid &, const wavefunction &);


};



#endif // hamop_h




