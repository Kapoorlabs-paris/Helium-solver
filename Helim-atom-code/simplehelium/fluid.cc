#include<fluid.h>
#include<hamop.h>
#include<wavefunction.h>
#include<grid.h>

#define pi 3.1415926

class wavefunction;



// *** now some operators follow
fluid &fluid::operator=(const fluid &v)
{
    if (this != &v)
    {
       delete[] start;
       fluid_dim  = v.fluid_dim;
       start = new double[fluid_dim];
       for (long i = 0; i < fluid_dim; i++)
           start[i] = v.start[i];
    }
    return *this;
}

fluid& fluid::operator *= (double z)
{
  for(long i=0; i<fluid_dim; i++)
    start[i]=start[i]*z;
  return *this;
}

fluid operator * (
			 double z, 
			 const fluid &v
			 )
{
  fluid temp=v;
  return temp *= z;
}

fluid operator * (
			 const fluid &v, 
			 double z
			 )
{
  fluid temp=v;
  return temp *= z;
}

double operator * (
			    const fluid &v, 
			    const fluid &w 
			    )
{
  double result=0.0;
  for(long i=0; i<v.wf_size(); i++)
    {
      result+=v[i]*w[i];
    };
  return result;
}

fluid operator + (
			    const fluid &v, 
			    const fluid &w 
			    )
{
  fluid result=v;
  for(long i=0; i<v.wf_size(); i++)
    {
      result[i]=v[i]+w[i];
    };
  return result;
}

fluid operator - (
			    const fluid &v, 
			    const fluid &w 
			    )
{
  fluid result=v;
  for(long i=0; i<v.wf_size(); i++)
    {
      result[i]=v[i]-w[i];
    };
  return result;
}


// *** and now the <<, >> -stuff
ostream& operator<<(
		    ostream& os, 
		    const fluid& v)
{  
  for(long i = 0; i < v.wf_size(); i++)
    {   os << v[i] << endl;
      //os << v[i] << '\t';
    //    if ((i+1)%8==0 || i==v.wf_size()-1) os << endl;
    }
  return os;
}


istream& operator>>(
		    istream& is, 
		    fluid& v
		    )
{  
  double tmpre;
  for(long i = 0; i < v.wf_size(); i++)
    {   
      is >> tmpre;
      v[i]=tmpre;
    }
  return is;
}

istream& operator>>(
		    ifstream& is, 
		    fluid& v
		    )
{  
  double tmpre;
  for(long i = 0; i < v.wf_size(); i++)
    {   
      is >> tmpre;
      v[i]=tmpre;
    }
  return is;
}



double fluid::check_total_charge(grid g) const
{
  long index, xindex, yindex;
  double result=0.0;
  double cosy, cosr, sinr, beta;

  beta=g.delt_z();
  
  switch (g.dimens())
    {

    case 89 :
      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{
	  for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      result+=start[index];
	    };
	};
      
      return result*6.283185307*g.delt_x()*g.delt_y();
      
      break;


    case 88 :

      for (xindex=1; xindex<g.ngps_x()-1; xindex++)
	{
	  for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	    {
	      index=g.index(xindex,yindex,0);
	      result+=g.rho(xindex)*start[index];
	    };
	};
      
      return result*6.283185307*g.delt_x()*g.delt_y();
      
      break;

    case 87 :
      
      xindex=0;
      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	{
	  index=g.index(xindex,yindex,0);
	  result+=start[index];
	};
      
      return result*g.delt_y();
      
      break;

    case 86 :

      xindex=0;
      for (yindex=1; yindex<g.ngps_y()-1; yindex++)
	{
	  index=g.index(xindex,yindex,0);
	  cosy=cos(g.y(yindex));
	  result+=start[index]/(cosy*cosy);
	};
      
      return result*g.delt_y()*beta;
      
      break;


    };
}


void fluid::dump_to_file(grid g,FILE* os, int stepwidth)
{

  long i;



      for (i=0; i<g.ngps_x(); i+=stepwidth)
	{ 
	  fprintf(os,"%.14le\n",start[i]);
	};

}


