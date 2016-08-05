#ifndef wavefunction_h
#define wavefunction_h wavefunction_h
#include<assert.h>
#include<complex>
#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>

using namespace std;

class fluid;
class hamop;
class grid;

class wavefunction
{
  public:
    wavefunction(long x=0) : wf_dim(x), start(new complex<double>[x]) { }

    wavefunction(const wavefunction& v) 
      {
	wf_dim  = v.wf_dim;
	start = new complex<double>[wf_dim];
	for (long i = 0; i < wf_dim; i++)
	  start[i] = v.start[i];
      }

    virtual ~wavefunction() { delete [] start;}

    long  wf_size() const {return wf_dim;}

    complex<double>&  operator[](long index)
    {
      //      assert(index >= 0  &&  index < wf_dim);
      return start[index];
    }

    const complex<double>&  operator[](long index) const   
    {   
      //      assert(index >= 0  &&  index < wf_dim);
      return start[index];
    }

    wavefunction& operator=(const wavefunction&);
    wavefunction& operator=(const fluid&);
    
    void init(grid g, int inittype, double width, double k, double time, FILE* filename, int output_of_interest);
    void init(grid g, int inittype, double width, double k, double time);
    complex<double>* begin() { return start;}
    complex<double>* end()   { return start + wf_dim;}
    const complex<double>* begin() const { return start;}
    const complex<double>* end()   const { return start + wf_dim;}

 wavefunction Static_KS(grid g);
    wavefunction& operator *= (double z); 
    wavefunction& operator *= (complex<double> z); 
 wavefunction clip(grid gone);
 wavefunction alog(grid g, double deltat);
    wavefunction dens1d(grid g);
    wavefunction phaseoutdft(grid g);
    wavefunction phaseinsch(grid g);
    wavefunction phasekohnshamorbital(grid gone, grid g, const wavefunction &wf,long box);
     wavefunction phasekohnshamorbitalnull(grid gone, grid g, const wavefunction &wf,long box);
    wavefunction denskohnsham(grid g);
    wavefunction denskohnshamcorrector(grid gone,grid g,const wavefunction &wfheliumplus,const wavefunction &wf);
    wavefunction approxdenscorector(grid gone,grid g,const wavefunction &wfheliumplus,const wavefunction &wf);
void even(grid g);
    void regrid(grid g, grid g_small, const wavefunction &v);
    void propagate(complex<double> delta_time, 
		   double time, 
		   grid g, 
		   hamop hamil, 
		   int me, 
		   int vecpotflag,
		   const wavefunction &staticpot_x,
		   const wavefunction &staticpot_y,
		   
		   const wavefunction &staticpot_xy,
		   double charge);
void oddinz(grid g);
    void calc_hartree_electronic(grid g, grid gi, hamop interactionhamil, const wavefunction &wfI);
   void propagatedft(complex<double> delta_time, 
		   double time, 
		   grid g, 
		   hamop hamil, 
		   int me, 
		   int vecpotflag,
		   const wavefunction &staticpot,
		   const wavefunction &tddftpot,
		   double charge);
wavefunction staticdftpot(grid g);
    complex<double> energy(double time, grid g, hamop hamil, int me, const double masses[], const wavefunction &staticpot_x, const wavefunction &staticpot_y, const wavefunction &staticpot_xy, double charge);
    complex<double> energydft(double time, grid g, hamop hamil,int me, const double masses[], const wavefunction &staticpot, const wavefunction &tddftpot, double charge);
    wavefunction  setmyboundary(grid g, hamop hamil, double time, int me);
    double expect_x(grid g);
    double expect_y(grid g);
    double expect_dipole(grid g);
    complex<double> expect_transition_dipole(grid g, const wavefunction &a);
    double expect_dipole_acceleration(grid g, double eps);
    wavefunction alpha1d(grid g);
    void  setzeroboundary(grid g);
wavefunction  forcezeroboundary(grid g);
    double norm(grid g);
    complex<double> floquet(grid g, const wavefunction &wf);
wavefunction timederivalpha(grid g, const wavefunction &wf, const wavefunction &dens, double deltat);
wavefunction notimederivalpha(grid g, const wavefunction &wf, const wavefunction &dens, double deltat);
    double non_ionized(grid g, long box);
    double sing_ionized(grid g, long box);
    double doub_ionized(grid g, long box);
    double molecule_ionized(grid g, long box);
double expect_xcb(grid g);
double expect_ycb(grid g);
double expect_xsq(grid g);
double expect_ysq(grid g);
double entropy(grid g);
    void nuclenergysurface(double time, grid g, hamop hamil,int me, const double masses[], const wavefunction &staticpot_x, const wavefunction &staticpot_y, const wavefunction &staticpot_xy, double charge, wavefunction &nuclpot);
    void solve(const wavefunction &a, const wavefunction &b, 
	       const wavefunction &c, const wavefunction &psi, long dimens);
    void solve_toep(double a, double b, double b_upperleft, double b_lowerright, 
	       double c, const wavefunction &psi, long dimens);

    void apply(const wavefunction &a, const wavefunction &b, 
	       const wavefunction &c, const wavefunction &psi);

    void dump_to_file(grid g, FILE* os, int stepwidth);
    void  calculate_fixed_potential_array(grid g, hamop hamil, double time, int me);
 void  calculate_fixed_potential_arraydft(grid g, hamop hamil, double time, int me);
    void dump_dens1d_to_file(grid g, FILE* os, int stepwidth);

    void dump_wf1d_to_file(grid g, FILE* os, int stepwidth);
    void initmine(grid g, int inittype, double width, double k, double time, FILE* filename, int output_of_interest);

    void nullify();

    void calculate_fixed_potential_array_x(grid g, hamop hamil, double time, int me);
    void calculate_fixed_potential_array_y(grid g, hamop hamil, double time, int me);
     void calculate_fixed_potential_array_z(grid g, hamop hamil, double time, int me);
    void calculate_fixed_potential_array_xy(grid g, hamop hamil, double time, int me);
wavefunction fracdenskohnshamcorrector(grid gone,const wavefunction &wfheliumplus,double epsilon);
wavefunction fracphasekohnshamorbital(grid gone,double epsilon);
wavefunction justalpha(grid g, const wavefunction &wf,const wavefunction &dens, double deltat);
wavefunction justalphax(grid g, const wavefunction &wf,const wavefunction &dens, double deltat);
wavefunction simplepotential(grid gone);
  private:
    long  wf_dim;
    complex<double>   *start; 

    complex<double> int_P0_lp(const long r_point, const long l_l, grid g);
    complex<double> int_P1_lp(const long l, const long l_l, grid g);
    complex<double> int_P2_lp(const long l, const long l_l, grid g);



    void do_cn_step_x(complex<double> timestep, double time, grid g, hamop hamil, int me, int vecpotflag, long yindex, long zindex, const fluid &wf_one, const wavefunction &wf_two);

    void do_cn_step_x_muller(complex<double> timestep, double time, grid g, hamop hamil, int me, const wavefunction &staticpot, const wavefunction &tddftpot, long yindex, long zindex);

    void do_cn_step_xy_muller(complex<double> timestep, double time, grid g, hamop hamil, int me, const wavefunction &staticpot_x, const wavefunction &staticpot_y, const wavefunction &staticpot_xy, long yindex, long zindex);

    void do_muller_crapola(complex<double> timestep, double time, grid g, hamop hamil, const wavefunction &staticpot, int me, fluid &wf_one, wavefunction &wf_two, long noofgridpoints[], double charge);

    void do_muller(complex<double> timestep, double time, grid g, hamop hamil, const wavefunction &staticpot, int me, long noofgridpoints[], double charge);

    void do_cn_step_xy(complex<double> timestep, double time, grid g, hamop hamil, int me, int vecpotflag, long zindex, const fluid &wf_one, const wavefunction &wf_two);

    void do_osc_coupling_prop(complex<double> timestep, double time, grid g, hamop hamil, int me, int vecpotflag, long zindex, const fluid &wf_one, const wavefunction &wf_two, const double masses[]);



};


ostream& operator<<(ostream& os, const wavefunction& v);
istream& operator>>(istream& is, wavefunction& v);

complex<double> operator * (const wavefunction &v, const wavefunction &w );
wavefunction operator * (double z, const wavefunction &v);
wavefunction operator * (complex<double> z, const wavefunction &v);
wavefunction operator * (const wavefunction &v, double z);
wavefunction operator * (const wavefunction &v, complex<double> z);
wavefunction operator / (const wavefunction &v, double z);
wavefunction operator / (const wavefunction &v, const wavefunction &w);
wavefunction operator | (const wavefunction &v, const wavefunction &w);

wavefunction operator - (const wavefunction &v, const wavefunction &w );
wavefunction operator + (const wavefunction &v, const wavefunction &w );
wavefunction operator + (const wavefunction &v, const fluid &w );
wavefunction operator + (const fluid &v, const wavefunction &w );
wavefunction operator + (const wavefunction &v, double z );

#endif // wavefunction_h






