// *** grid.h ***  for the hop_cart project -----------------  db 09/98
#ifndef grid_h
#define grid_h grid_h

#include<iostream>


using namespace std;

class grid
{
  
 public:

  grid(int a=0, long b=0, double c=0.0, long d=0)
    {
      dimension=a;
      n_gps_x=b;
      n_gps_y=b;
      n_gps_z=b;
      delta_x=c;
      delta_y=c;
      delta_z=c;
      offse_x=d;
      offse_y=d;
      offse_z=d;
    }

  inline int dimens() const { return dimension; }

  inline long ngps_x() const { return n_gps_x; }
  inline long ngps_y() const { return n_gps_y; }
  inline long ngps_z() const { return n_gps_z; }
  
  inline double delt_x() const { return delta_x; }
  inline double delt_y() const { return delta_y; }
  inline double delt_z() const { return delta_z; }
  
  inline long offs_x() const { return offse_x; }
  inline long offs_y() const { return offse_y; }
  inline long offs_z() const { return offse_z; }
  
  inline double volelem() const { return delta_x*delta_y*delta_z; }
  inline long totgps() const { return n_gps_x*n_gps_y*n_gps_z; }

  inline void set_dim(int a) { dimension=a; }

  inline void set_ngps(long a, long b, long c)
    {
      n_gps_x=a;
      n_gps_y=b;
      n_gps_z=c;
    }

  inline void set_delt(double a, double b, double c)
    {
      delta_x=a;
      delta_y=b;
      delta_z=c;
    }

  inline void set_offs(long a, long b, long c)
    {
      offse_x=a;
      offse_y=b;
      offse_z=c;
    }

  inline void showit() const
    {
      cout << dimension << endl << n_gps_x << ", " << n_gps_y << ", "
	   << n_gps_z << endl << delta_x << ", " << delta_y << ", "
	   << delta_z << endl << offse_x << ", " << offse_y << ", "
	   << offse_z << endl;
    }

  inline long index(long xindex, long yindex, long zindex)
    {
      long result;
      result=zindex*(ngps_x()*ngps_y())+yindex*ngps_x()+xindex;
      return result;
    }







  inline double r(long index)
    {
      return (index+1)*delta_x;
    }

  inline double rho(long index) // Attention! This is for TF where we have
                                      // extra points for the boundaries.
     {
       return (index-0.5)*delta_x;
    }


  inline double x(long index)
    {
      return (index-offse_x+0.5)*delta_x;
    }
  
  inline double y(long index)
    {
      return (index-offse_y+0.5)*delta_y;
    }
  
  inline double z(long index)
    {
      return (index-offse_z+0.5)*delta_z;
    }

  inline long rindex(double r)
    {
      return long( r/delta_x - 0.5);
    }

  inline long rhoindex(double rho)
    {
      return long( rho/delta_x + 0.5);
    }

  inline long yindex(double y)
    {
      return long( y/delta_y - 0.5 + offse_y);
      //return long( y/delta_y + offse_y);
    }

  inline long xindex(double x)
    {
      return long( x/delta_x - 0.5 + offse_x);
    }

  
 private:
  
  int dimension;
  long n_gps_x, n_gps_y, n_gps_z;
  double delta_x, delta_y, delta_z;
  long offse_x, offse_y, offse_z;

};



#endif 





