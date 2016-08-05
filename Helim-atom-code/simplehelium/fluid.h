#ifndef fluid_h
#define fluid_h fluid_h
#include<assert.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>

using namespace std;

class grid;
class hamop;
class propmanage;
class wavefunction;


class fluid
{
 public:
  fluid(long x=0) : fluid_dim(x), start(new double[x]) { }
  
  fluid(const fluid& v) 
    {
      fluid_dim  = v.fluid_dim;
      start = new double[fluid_dim];
      for (long i = 0; i < fluid_dim; i++)
	start[i] = v.start[i];
    }
  
  virtual ~fluid() { delete [] start;}
  
  long  wf_size() const {return fluid_dim;}
  
  double&  operator[](long index)
    {
      //      assert(index >= 0  &&  index < fluid_dim);
      return start[index];
    }
  
  const double&  operator[](long index) const   
    {   
      //      assert(index >= 0  &&  index < fluid_dim);
      return start[index];
    }
  
  fluid& operator=(const fluid&);
  
  double* begin() { return start;}
  double* end()   { return start + fluid_dim;}
  const double* begin() const { return start;}
  const double* end()   const { return start + fluid_dim;}
  
  fluid& operator *= (double z); 
  
  
  double check_total_charge(grid g) const;
  

  void dump_to_file(grid g,FILE* os, int stepwidth);
  
    
 private:
  long  fluid_dim;
  double   *start; 
  
  
};


ostream& operator<<(ostream& os, const fluid& v);
istream& operator>>(istream& is, fluid& v);
istream& operator>>(ifstream& is, fluid& v);
double operator * (const fluid &v, const fluid &w );
fluid operator * (double z, const fluid &v);
fluid operator * (const fluid &v, double z);
fluid operator + (const fluid &v, const fluid &w );
fluid operator - (const fluid &v, const fluid &w );

#endif // fluid_h






