//----------------------------------------------------------------------------//
// “очный метод численного интегрировани€
//----------------------------------------------------------------------------//
struct SatOrbitPoint
{
	double posx;
	double posy;
	double posz;
	double posr;

	double vx;
	double vy;
	double vz;

	int flag;
};

SatOrbitPoint RunGRK( double x, double y, double z, double vx, double vy, double vz, double Tstart, double Tend, double inh, double EEE );

void rk8cc( int n, double t0, double *x0, double t, double h, double *p, int *np, int kp, double *e, int (*fun) (double,double *,double *) );

//----------------------------------------------------------------------------//