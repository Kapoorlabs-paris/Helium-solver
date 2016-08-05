
 
#include<TDSE.h>
#include<wavefunction.h>
#include<fluid.h>
#include<grid.h>
#include<hamop.h>

#define pi 3.1415926

int main(int argc, char **argv)
{
  int me=0;



  	FILE *file_wfdat,*file_wf1dat,*file_reading;
  	FILE *file_obser, *file_obser_imag;
	// FILE *file_reading,*file_readingexcited ;    // *** create some files with appropriate appendices
  
  	char string_wfdat[]=  "res/wf_laser_25.dat";
  	char string_obser[]=       "res/obser_laser_25.dat";
  	char string_obser_imag[]=  "res/obserimag_laser_25.dat";
       
	char string_reading[]=  "res/wf_2000_01.dat";
   	
  
  ///////////////////////////// Turn on reading once you generate ground state file to be read from, 
  //then uncomment the readin in string declaration,filename and in fopen ////////////////////
  
     	file_wfdat = fopen(string_wfdat,"w");
    	//file_wf1dat = fopen(string_wf1dat,"w");

  	file_obser = fopen(string_obser,"w");
  	file_obser_imag = fopen(string_obser_imag,"w");
  	file_reading = fopen(string_reading,"r");
	//file_readingexcited = fopen(string_readingexcited,"r");

  long index, rrindex, xindex, yindex,zindex,index2;
  complex<double> imagi(0.0,1.0);
  double deltx=0.2;
  double delty=0.2;
  double deltz=1;
  long ngpsx=1500;
  long ngpsy=1500;
  long ngpsz=1;

  // *** declare grid ***
  grid g;
  g.set_dim(16);            		// propagation mode is 3 for xy cartesian
  g.set_ngps(ngpsx,ngpsy,ngpsz);    	// N_x, N_y, 1 was 100,100,1
  g.set_delt(deltx,delty,deltz);  	// delta_x, delta_y, deltz was 0.1,0.1,0.1
  g.set_offs(ngpsx/2,ngpsy/2,0);    	// origin  (usually at N_x/2, N_y/2,N_z/2)
  
  // *** declare smaller grid for reading***
  grid g_small;
  g_small.set_dim(16);
  g_small.set_ngps(ngpsx/2,ngpsy/2,ngpsz/2);
  g_small.set_delt(deltx,delty,deltz);
  g_small.set_offs(ngpsx/4,ngpsy/4,0);
  
  // *** declare rest of variables ***
  double imag_timestep=0.25;
  double real_timestep=0.1;
  long   no_of_imag_timesteps=1000;
  long   no_of_real_timesteps=20000;
  int    obs_output_every=1;
  long   wf_output_every=20000; 
  int    dumpingstepwidth=1;
  int    vecpotflag=1;   		// don't touch this
  int    box=50;          		//was 50
  double masses[]={1.0,1.0};    	// don't touch this
  double charge=0.0;    		// don't touch this
  double epsilon=0.00001;
  hamop interactionhamil(g,vecpot_x,vecpot_y,vecpot_z,interactionpotxy,scalarpoty,scalarpotz,interactionpotxy,imagpotx,imagpoty,field,dftpot);
  hamop hamilton(g,vecpot_x,vecpot_y,vecpot_z,scalarpotx,scalarpoty,scalarpotz,interactionpotxy,imagpotx,imagpoty,field,dftpot);
 
    	wavefunction wf(g.ngps_x()*g.ngps_y()*g.ngps_z());   
   	wavefunction wf1(g.ngps_x()*g.ngps_y()*g.ngps_z()); 
   	wavefunction wfread(g.ngps_x()*g.ngps_y()*g.ngps_z());

  	wavefunction wfini(g.ngps_x()*g.ngps_y()*g.ngps_z());
  	wavefunction wfeverythingeven(g.ngps_x()*g.ngps_y()*g.ngps_z());
  	wavefunction wfeverythingodd(g.ngps_x()*g.ngps_y()*g.ngps_z());
 
  complex<double> timestep;
  double time=0.0;
  
  long counter_i=0;
  long counter_ii=0;
 

  ////////////////////////// To get first excited spin singlet state two orthognalizations are needed in imaginary propagation and so on (two more and then more) to obtain higher spin singlet states //////////////////////
  
	// initialization

  long outputofinterest=0;
  
  
	wf.init(g,3,2.0,0.0,0.0);
      		// ground
	//wfread.init(g,99,0.1,0.0,0.0,file_reading,outputofinterest);
  	//wf.nullify();
  	//wf.regrid(g,g,wfread);

	wf*=1.0/sqrt(wf.norm(g));
		//wfground=wf;
		//wf=wfground;

	//wf1.init(g,3,2.0,0.0,0.0);
	//wf1*=1.0/sqrt(wf1.norm(g));


  fclose(file_reading);
 
  cout << "norm wf    : " << wf.norm(g) <<  "\n";

  

  wavefunction staticpot_x(g.ngps_x());
  staticpot_x.calculate_fixed_potential_array_x(g,hamilton,0.0,me);
  

  wavefunction staticpot_y(g.ngps_y());
  staticpot_y.calculate_fixed_potential_array_y(g,hamilton,0.0,me);
  
    
  wavefunction staticpot(g.ngps_x()*g.ngps_y()*g.ngps_z());
  staticpot.calculate_fixed_potential_array(g,hamilton,0.0,me);
   

  wavefunction staticpot_xy(g.ngps_z()*g.ngps_y()*g.ngps_x());
  staticpot_xy.calculate_fixed_potential_array_xy(g,hamilton,0.0,me);
  


 complex<double> dftapproxenerg;
 complex<double> complenerg;
 complex<double> groundstatepop, excitedstatepop;

	// ************* imag timeprop

  long ts;
  long no_of_timesteps=no_of_imag_timesteps;
  for (ts=0; ts<no_of_timesteps; ts++)
    {   
                                                                                  
      cout << "Imag: " << ts <<"  "  <<  " energy  "  << real(complenerg) << "  "  << endl; 
          
      counter_i++;
      counter_ii++;
      timestep=complex<double>(0.0*real_timestep,-1.0*imag_timestep);  
      time=-imag(timestep*(complex<double>)(ts));
      
      	// and now the actual propagation

  	wf.propagate(timestep,0.0,g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);
	wf*=1.0/sqrt(wf.norm(g));

	//wf.propagate(timestep,0.0,g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);
	//wf=wf-((wf1*wf)*g.delt_x()*g.delt_y())*wf1;
	//wf*=1.0/sqrt(wf.norm(g));
      
      if (counter_ii==obs_output_every)
      	{
		complenerg=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge);

	  	  fprintf(file_obser_imag,"%li %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le\n",
			  ts,real(complenerg),imag(complenerg),real(complenerg),wf.norm(g),wf.expect_x(g),wf.expect_y(g),wf.doub_ionized(g,box),wf.expect_x(g)); //wf.doub_ionized(g,box),wf.expect_x(g),wf.expect_y(g)
	  counter_ii=0;
	 
	};
      
    };
  
  
  fclose(file_obser_imag);
  
	wfini=wf;
  
      counter_i=0;
      counter_ii=0;
    
	wf=wfini;
    
	// ************* real timeprop
  
  timestep=complex<double>(real_timestep,0.0);
  no_of_timesteps=no_of_real_timesteps;
  
  for (ts=0; ts<no_of_timesteps; ts++)
    {
      counter_i++;
      counter_ii++;
      time=real(timestep*(complex<double>)(ts));
        cout << "Real: " << ts << "  "  << " erwart_x wf : " << wf.expect_x(g)  << " " << " Energie " << real(complenerg) << " " << endl; //\n";
   
      wf.propagate(timestep,time,g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);   

   if (counter_ii==obs_output_every)
      	{
	  complenerg=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge);
	     
	       groundstatepop=wf*wfini*g.delt_x()*g.delt_y();


	  	  fprintf(file_obser," %.14le %.14le %.14le %.14le %.14le  %.14le  %.14le %.14le\n ",
	    (time+real(timestep)),real(complenerg), imag(complenerg), wf.norm(g), wf.expect_x(g), wf.expect_y(g), hamilton.vecpot_x(time,me), real(conj(groundstatepop)*groundstatepop));
	  
	  counter_ii=0;
	};
    
 
/*
  if (counter_i==wf_output_every)
           {
     	wf.dump_to_file(g,file_wfdat,dumpingstepwidth);
	counter_i=0;
};
*/

};

wf.dump_to_file(g,file_wfdat,dumpingstepwidth);
	
    fclose(file_obser);
   
    fclose(file_wfdat);
   
    cout << me << ": Hasta la vista, ... " << endl;

 }

	//potentials

double vvvvecpot_x(double time, int me)
{
  double result=0.0;
double alphahat=0.1;
  double frequ=0.1;
  double n=100.0;
  double ampl=alphahat*frequ;
  double phi=0.0;
  double dur=n*2.0*M_PI/frequ;
double gap=2.0*M_PI/frequ;
 double ww=0.5*frequ/n;
  

 
   if ((time>0.0)  )
    {
	result=-0.001;
	//result=ampl*sin(ww*time)*sin(ww*time)*sin(frequ*time-phi);
 };

  return result;

}  


double vecpot_x(double time, int me)
{

  double result=0.0;

  double frequ= 1.556;
 
 
 
double alphahat=0.1;
  
  double n=400.0;
  double ampl=alphahat*frequ;
  
  double dur=n*2.0*M_PI/frequ;

 double ww=0.5*frequ/n;
 

 
   if ((time>0.0))
    {
	//result=0.0;
	//result=-0.001;
	result=ampl*sin(ww*time)*sin(ww*time)*sin(frequ*time);
};
 

	
	/*

 
  if (time>0)
    {

      if (time<ramping){
	result=-ampl/ramping*time*cos(frequ*time) + ampl/(frequ*ramping)*sin(frequ*time);
	}
      
      if ( (time>=ramping)&&(time<ramping+constperiod))
      {
	result=-ampl*cos(frequ*time);//-ampl*cos(frequ2*time);
      }
     
      if ( (time >= ramping+constperiod) && (time < ramping+constperiod+downramp)) 
	{
		  result=ampl/downramp*(time-downramp-constperiod-ramping)*cos(frequ*time)-ampl/(frequ*downramp)*sin(frequ*time);
	}
      
      if ( (time >= ramping+constperiod+downramp )) //&& (time<ramping+constperiod+downramp +gap)) 
      {
      result=0.0;
 }
  };
 
*/
    
  return result;
}  



double vvvvvvecpot_x(double time, int me)
{

  double result=0.0;
  return result;
}  


double alpha_y(double time, int me)
{
  double result=0.0;
  double vecpotampl, frequref;
  double ramping, constperiod, downramp;

   // <-------- put same in vecpot_y !!!
  vecpotampl = 0;


  result=0;

  return result;

}  


double alpha_x(double time, int me)
{
 double result=0.0;
  double vecpotampl, frequref;
  double ramping, constperiod, downramp;

  // <-------- put same in vecpot_y !!!
  vecpotampl =0;


  result=0;

  return result;

 

}  

double alpha_z(double time, int me)
{
 

  return alpha_y(time,me);

}  


double vecpot_y(double time, int me)
{
 return vecpot_x(time,me);
}  


double vecpot_z(double time, int me)

{

 
 
    

  return 0.0;
 
 // return vecpot_x(time,me);
}  





double scalarpotx(double x, double y, double z, double time, int me)
{
  double eps=1.0;
 double result;
 
 //result= 0.5*x*x-10; // harmonic
 result= -2.0/sqrt(x*x+eps*eps); // coulomb  
	return result;
}

double scalarpoty(double x, double y, double z, double time, int me)
{
  double eps=1.0;
    
double result;

 //result= 0.5*y*y-10; // harmonic
 result= -2.0/sqrt(y*y+eps*eps);  // coulomb
      return result;

}

double scalarpotz(double x, double y, double z, double time, int me)
{
double eps=1.0;

 return 0.0;

}

double interactionpotxy(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  
  return 1.0/sqrt((x-y)*(x-y)+eps*eps); //coulomb
  //return -0.5*(x-y)*(x-y); // harmonic
  

}

double interactionpotyz(double x, double y, double z, double time, int me)
{
  double eps=1;
   return 0.0;  
 
 
}
double interactionpotxz(double x, double y, double z, double time, int me)
{
  double eps=1;
    return 0.0;  
 
 
}
  
double field(double time, int me)
{
  double result=0.0;


  return result;


}  







double imagpotx(long xindex, long yindex, long zindex, double time, grid g)
{
  double x,y,z;
  
 //double ampl=0.0; // switch imaginary potential off
 double ampl=50.0; // switch imaginary potential on


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
  //double ampl=0.0; // switch imaginary potential off
  
   double ampl=50.0; // switch imaginary potential on


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

double imagpotz(long xindex, long yindex, long zindex, double time, grid g)
{
  double z;
  // double ampl=0.0; // switch imaginary potential off
   double ampl=0.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       z=((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z())
	    *((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z());
       return ampl*z*z*z*z*z*z*z*z;
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




double dfthartree(grid gone, double x, double y, double z, double time, int me, 
	      const wavefunction & wf, double eps, long box)
{
  double result=0.0;
  return result;
}


