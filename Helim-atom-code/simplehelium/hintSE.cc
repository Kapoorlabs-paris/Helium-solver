#include<TDSE.h>
#include<wavefunction.h>
#include<fluid.h>
#include<grid.h>
#include<hamop.h>

#define pi 3.1415926


int main(int argc, char **argv)
{
  int me=0;



  FILE *file_wfdat,*file_alphasch,*file_alphadftapproxdat,*file_alphadftapproxdat;
  FILE *file_pot, *file_wf1d, *file_alphasch,*file_wfdftexacttestdat;
  FILE *file_obser, *file_obser_imag;
 // FILE *file_reading;
  
  
  // *** create some files with appropriate appendices
  
  char string_wfdat[]=       "/home/theo/kapoor/d34/res/florian/wftest.dat";
  char string_pot[]=        "/home/theo/kapoor/d34/res/florian/pottest.dat";
  char string_wf1d[]=       "/home/theo/kapoor/d34/res/florian/wf1dtest.dat";
  char string_alphasch[]=       "/home/theo/kapoor/d34/res/florian/alphasch.dat";
   char string_wfdftexactdat[]=      "/home/theo/kapoor/d34/res/florian/wfdftexacttest.dat";
char string_alphadftapproxdat[]=      "/home/theo/kapoor/d34/res/florian/alphadftapproxtest.dat";
char string_alphadftexactdat[]=      "/home/theo/kapoor/d34/res/florian/alphadftexacttest.dat";
char string_obser[]=       "/home/theo/kapoor/d34/res/florian/observtest.dat";
  char string_obser_imag[]=  "/home/theo/kapoor/d34/res/florian/observimagtest.dat";
    
  
  file_wfdat = fopen(string_wfdat,"w");
  file_alphasch = fopen(string_alphasch,"w");
 file_wfdftexactdat = fopen(string_wfdftexactdat,"w");
   file_alphadftapproxdat = fopen(string_alphadftapproxdat,"w");
 file_alphadftexactdat = fopen(string_alphadftexactdat,"w");
file_pot = fopen(string_pot,"w");
  file_wf1d = fopen(string_wf1d,"w");
  file_obser = fopen(string_obser,"w");
  file_obser_imag = fopen(string_obser_imag,"w");
  //file_reading = fopen(string_reading,"r");
  
  
  long  yindex;
  
  double deltx=0.2;
double delty=0.2;
  
  long ngpsx=80;
  long ngpsy=80;
  // *** declare grid ***
  grid g;
  g.set_dim(16);            // propagation mode is 16 for xy cartesian
  g.set_ngps(ngpsx,ngpsy,1);    // N_x, N_y, 1 was 100,100,1
  g.set_delt(deltx,delty,1.0);  // delta_x, delta_y, 1.0 was 0.5,0.5,1.0
  g.set_offs(ngpsx/2,ngpsy/2,0);    // origin  (usually at N_x/2, N_y/2)
  
  // *** declare smaller grid for reading***
  grid g_small;
  g_small.set_dim(16);
  g_small.set_ngps(ngpsx,ngpsx,1);
  g_small.set_delt(deltx,delty,1.0);
  g_small.set_offs(ngpsx/2,ngpsy/2,0);
  
  // *** declare grid for 1D ***
  grid g1d;
  gtwo.set_dim(15);
  gtwo.set_ngps(ngpsx,1,1);
  gtwo.set_delt(deltx,1.0,1.0);
  gtwo.set_offs(ngpsx/2,0,0);
  
  
  // *** declare rest of variables ***
  double imag_timestep=0.2;
  double real_timestep=0.05;
  long   no_of_imag_timesteps=1000;
  long   no_of_real_timesteps=10000;
  int    obs_output_every=1;
  long   wf_output_every=100; 
  int    dumpingstepwidth=1;
  int    vecpotflag=1;   // don't touch this
  int    box=15;          //was 50
  double masses[]={1.0,1.0};    // don't touch this
  double charge=0.0;    // don't touch this
  
  hamop hamilton(g,vecpot_x,vecpot_y,vecpot_z,scalarpotx,scalarpoty,scalarpotz,interactionpotxy,imagpotx,imagpoty,field,dftpot);
  hamop hamiltontwo(gtwo,vecpot_x,vecpot_y,vecpot_z,scalarpotxtwo,scalarpotytwo,scalarpotztwo,interactionpotxytwo,imagpotxtwo,imagpotytwo,field,dftpot);


  wavefunction wf(g.ngps_x()*g.ngps_y()*g.ngps_z());
 
  wavefunction wf1d(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  wavefunction lowerimag(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  wavefunction upper(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  wavefunction lower(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  wavefunction potential(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  wavefunction realpot(g1d.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  wavefunction wfdftexact(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
 wavefunction pot(g1d.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());


    pot.nullify();
  
  complex<double> timestep;
  double time=0.0;
  
  long counter_i=0;
  long counter_ii=0;
  long counter_iii=0;

  long testindex1=0;
  long testindex2=0;
  long testindex3=0;

  
  // initialization
  long outputofinterest=0;
  // ground
 //  wfread.init(g_small,99,0.1,0.0,0.0,file_reading,outputofinterest);
//   wf.nullify();
//   wf.regrid(g,g_small,wfread);
//   wf*=1.0/sqrt(wf.norm(g));
//   wfground=wf;




  wf.init(g,1,7.0,0.0,0.0);
  wf*=1.0/sqrt(wfgs.norm(g));




  fclose(file_reading);
  

  cout << "norm wf    : " << wf.norm(g) << "\n";

  

  wavefunction staticpot_x(g.ngps_x());
  staticpot_x.calculate_fixed_potential_array_x(g,hamilton,0.0,me);
  
  wavefunction staticpot_y(g.ngps_y());
  staticpot_y.calculate_fixed_potential_array_y(g,hamilton,0.0,me);
  
  wavefunction staticpot_xy(g.ngps_x()*g.ngps_y()*g.ngps_z());
  staticpot_xy.calculate_fixed_potential_array_xy(g,hamilton,0.0,me);
  
  wavefunction staticpotdft(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  staticpot1d.calculate_fixed_potential_array(gtwo,hamiltontwo,0.0,me);

  wavefunction tddftpot(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  tddftpot.nullify();


  complex<double> complenerg;
  complex<double> groundstatepop, excitedstatepop;

  // ************* imag timeprop
  long ts;
  long no_of_timesteps=no_of_imag_timesteps;
  for (ts=0; ts<no_of_timesteps; ts++)
    {   
      
      cout << "Imag: " << ts << endl;
      
      counter_i++;
      counter_ii++;
      timestep=complex<double>(0.0*real_timestep,-1.0*imag_timestep);  
      
      
      // and now the actual propagation
      wf.propagate(timestep,0.0,g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);

      wf*=1.0/sqrt(wfgs.norm(g));

      
      

  alphasch=wf.alphasch(g);
 
     
      
      if (counter_ii==obs_output_every)
      	{
 	  complenerg=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge); 
	  
	  fprintf(file_obser_imag,"%li %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le\n",
		  ts,real(complenerg),imag(complenerg),wf.norm(g),wf.non_ionized(g,box),wf.sing_ionized(g,box),wf.doub_ionized(g,box),wf.expect_x(g),wf.expect_y(g));
	  counter_ii=0;
	  cout << real(complenerg) << endl;
	};
      
    };
  
  
  fclose(file_obser_imag);
  
   alphasch.dump_to_file(g_only_x,file_alphaschdat,dumpingstepwidth);
  
  counter_i=0;
  counter_ii=0;
  int coun_dip=0;
  


  // ************* real timeprop
  
  timestep=complex<double>(real_timestep,0.0);
  
  no_of_timesteps=no_of_real_timesteps;


  wf.dump_to_file(g,file_wfdat,dumpingstepwidth);

  wf1d=wf.dens1d(g);
  wf1d.dump_to_file(gtwo,file_wf1d,dumpingstepwidth);


  complenerg=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge);
	 

  wfdftexact=wf1d;

  for (ts=0; ts<no_of_timesteps; ts++)
    {
      counter_i++;
      counter_ii++;
      time=real(timestep*(complex<double>)(ts));
      
      cout << "Real: " << ts << "  (" << me << ")" << endl;
      
      wf.propagate(timestep,time+0.5*real(timestep),g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);   
      
     
	  
	  counter_ii=0;
	};

      lowerimag=wf1d;
      wf1d=wf.dens1d(g);


      upper=wf1d;
      lower=lowerimag;

      upper.propagate(-0.5*timestep,time+0.5*real(timestep),gtwo,hamilton1d,me,vecpotflag,staticpotdft,pot,charge);
      lower.propagate(0.5*timestep,time+0.5*real(timestep),gtwo,hamilton1d,me,vecpotflag,staticpotdft,pot,charge);

      
      potential=upper/lower;


      realpot=potential.arccos(gtwo,real(timestep));


      wfdftexact.propagate(timestep,time+0.5*real(timestep),gtwo,hamilton1d,me,vecpotflag,staticpotdft,realpot,charge);


      if (counter_i==wf_output_every)
       	{
	  wf.dump_to_file(g,file_wfdat,dumpingstepwidth);

	  wf1d.dump_to_file(gtwo,file_wf1d,dumpingstepwidth);
	  wfdftexact.dump_to_file(gtwo,file_wfdftexactdat,dumpingstepwidth);
	  realpot.dump_to_file(gtwo,file_pot,dumpingstepwidth);

	  counter_i=0;
 	};
      
       if (counter_ii==obs_output_every)
      	{
	  complenerg=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge);
	  

	


	  cout << real(complenerg) << endl;

	  fprintf(file_obser,"%.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le\n",
		  (time+real(timestep)),
		  real(complenerg),
		  imag(complenerg),
		  wf.norm(g),
		 
		  );
      
      
      };
 alphasch=wf.alphasch(g);
  alphadftexact=wfdftexact.alphadft(gtwo);
    wf.dump_to_file(g,file_wfdat,dumpingstepwidth);

    wf1d=wf.dens1d(g);
    wf1d.dump_to_file(gtwo,file_wf1d,dumpingstepwidth);
  alphasch.dump_to_file(gtwo,file_alphaschdat,dumpingstepwidth);
  alphadftexact.dump_to_file(gtwo,file_alphadftexactdat,dumpingstepwidth);






    fclose(file_obser);
    fclose(file_wfdat);
    fclose(file_dftexactdat);
    fclose(file_pot);
    fclose(file_wf1d);

    cout << me << ": Hasta la vista, ... " << endl;

 

 }



double vecpot_x(double time, int me)
{
 

  frequref=0.057; // <--------- put same in alpha_x !!!

  vecpotampl = 2.0;//;alphahat*frequref;


  n=3.0;


  ww=frequref/(2.0*n);
  dur=n*2*M_PI/frequref;
 
  if ((time>0.0) && (time<dur))
    {
      result=-vecpotampl*sin(ww*time)*sin(ww*time)*sin(frequref*time);
    };
  
  return result;
}  


double vecpot_xxx(double time, int me)
{

     double result;
     double vecpotampl, frequref,frequofinterest;
  double ramping, constperiod, downramp;

   frequref=0.058;

   vecpotampl =2*frequref;
    //alphahat*frequref;



   ramping=1.0*2*pi/frequref;
   constperiod=1*2*pi/frequref;
   downramp=1.0*2*pi/frequref;


   // the first pulse with the reference frequency

   if (time<ramping)
     {
       result=-vecpotampl/ramping*time*cos(frequref*time) + vecpotampl/(frequref*ramping)*sin(frequref*time);
     };
   if ( (time >= ramping) && (time < ramping+constperiod) ) 
     {
       result=-vecpotampl*cos(frequref*time);
     }; 
   if ( (time >= ramping+constperiod) && (time < ramping+constperiod+downramp) ) 
     {
       result=vecpotampl/downramp*(time-downramp-constperiod-ramping)*cos(frequref*time)-vecpotampl/(frequref*downramp)*sin(frequref*time);
     }; 
   if ( time > ramping+constperiod+downramp )
     {
       result=0.0;
     };

}  



double vecpot_xx(double time, int me)
{

  double result=0.0;
  return result;
}  





double vecpot_y(double time, int me)
{

 double result=0;
  double vecpotampl, frequref;
  double ramping, constperiod, downramp;
  double ww, dur, n;

  frequref=0.057; // <--------- put same in alpha_x !!!

  vecpotampl = 2.0;//;alphahat*frequref;


  n=3.0;


  ww=frequref/(2.0*n);
  dur=n*2*M_PI/frequref;
 
  if ((time>0.0) && (time<dur))
    {
      result=-vecpotampl*sin(ww*time)*sin(ww*time)*sin(frequref*time);
    };
  
  return result;
  

}  


double vecpot_z(double time, int me)

{
  return 0.0;
}  

double scalarpotx(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  return -2.0/sqrt(x*x+eps);
  // return 0.0;
}

double scalarpoty(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  return -2.0/sqrt(y*y+eps);
  // return 0.0;
}

double scalarpotz(double x, double y, double z, double time, int me)
{
  double result=0.0;
  return result;
}

double interactionpotxy(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  return 1.0/sqrt((x-y)*(x-y)+eps);
  //  return 0.0; // attention --- no interaction !!!!!!!!!!!!!!
}

  
double field(double time, int me)
{
  double result=0.0;


  return result;


}  

double imagpotx(long xindex, long yindex, long zindex, double time, grid g)
{
  double x,y,z;
  double ampl=0.0; // switch imaginary potential off
  //  double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
      x=((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x())
	*((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x());
  
      return ampl*x*x*x*x*x*x*x*x;
    }
  else
    {
       return 0.0;
    };
      
}  


double imagpoty(long xindex, long yindex, long zindex, double time, grid g)
{
  double y;
  double ampl=0.0; // switch imaginary potential off
  //  double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       y=((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y())
	    *((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y());
       return ampl*y*y*y*y*y*y*y*y;
    }
  else
    {
       return 0.0;
    };
      
} 



// dft is not used
double dftpot(grid g, double x, double y, double z, double time, int me, 
	      const fluid &v_null, const wavefunction &v_eins)
{
  double result;
  result=0.0;
  return result;
}



// 1d stuff


double scalarpotxtwo(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  return -2.0/sqrt(x*x+eps);
  // return 0.0;
}

double scalarpotytwo(double x, double y, double z, double time, int me)
{
//   double eps=1.00;
//   return -2.0/sqrt(y*y+eps);
  return 0.0;
}

double scalarpotztwo(double x, double y, double z, double time, int me)
{
  double result=0.0;
  return result;
}

double interactionpotxytwo(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  //  return 1.0/sqrt((x-y)*(x-y)+eps);
  return 0.0; // attention --- no interaction !!!!!!!!!!!!!!
}

  

double imagpotxtwo(long xindex, long yindex, long zindex, double time, grid g)
{
  double x,y,z;
  double ampl=0.0; // switch imaginary potential off
  //  double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
      x=((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x())
	*((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x());
  
      return ampl*x*x*x*x*x*x*x*x;
    }
  else
    {
       return 0.0;
    };
      
}  


double imagpotytwo(long xindex, long yindex, long zindex, double time, grid g)
{
  double y;
  double ampl=0.0; // switch imaginary potential off
  // double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       y=((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y())
	    *((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y());
       return ampl*y*y*y*y*y*y*y*y;
    }
  else
    {
       return 0.0;
    };
      
} 
