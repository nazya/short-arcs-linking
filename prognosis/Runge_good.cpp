//-----------------------------------------------------------------------//
// ��������� �.�.
// �������������� ��������� ��������
//-----------------------------------------------------------------------//
#include <math.h>
#include <stdio.h>
#include "GRK8.h"
#include "Runge_good.h"

#include "Influence/InfluenceForce.h"

Force::InfluenceForce *IForce;

//-----------------------------------------------------------------------//
// ���������� ����� �������
//-----------------------------------------------------------------------//
double GetDistTwoPoints( double x1, double y1, double z1, double x2, double y2, double z2 )
{
	double dist = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) );
	return dist;
}
//-----------------------------------------------------------------------//
// ���������� ������ �����
//-----------------------------------------------------------------------//
int RightFxyzv( double t, double *xt, double *fx )
{
	// ���� ������ ����������, ��������, �������� ������
	// ������� ���������� ���� ��� GPU
	//IF->rh_fun_grav( t, xt, fx );

	// ���� ���� �����������
	IForce->rh_fun( t, xt, fx );

	// ������������ ��������
	//IF->rh_fun_kepler( t, xt, fx );

	return 0;
};

int RightFxyzvAp( double t, double *xt, double *fx)
{
	// ���� ������ ����������, ��������, �������� ������
	// ������� ���������� ���� ��� GPU
	//IForce->rh_fun_grav( t, xt, fx );

	// ���� ���� �����������
	//IForce->rh_fun( t, xt, fx );

	// ������������ ��������
	IForce->rh_fun_kepler( t, xt, fx );

	return 0;
}

void  GetTELPos( double date1, double time1, double *Tpos, double *Telicrf )
{
	//-------------------------------------------------------//
	// ��������� �������
	double int1, ajd1, delt1;
			
	IForce->set_time(date1, time1, &ajd1, &delt1, &int1 );

	// ������� �������� � ������ �������
	double Arot[9];
	IForce->iers_update_matrix( int1, Arot, ajd1, delt1 );
	// ������� �������� �� ������ � ���������� �������
	double invArot[9];
	IForce->transpose( Arot, invArot ); 

	// ������� � ICRF ��������� ���������
	//double Telicrf[3];
	IForce->matVecMul(  invArot, Tpos, Telicrf );
	
	/* Output omitted
	printf( "Station Position icrf\n" );
	printf( "%f\t %f\t %f\n", Telicrf[0], Telicrf[1], Telicrf[2] );
	*/
}

//==============================================================================//
// ������� ������������������ ��������� � ���������� RA DEC
//==============================================================================//
void ConvertXYZtoRADEC( double *resultPosition, double *inTelescopePosition, double *Ra, double *Dec )
{
	// ������� ������ �������� � ICRF
	// ������ ����������� � ������� ICRF
	double x = resultPosition[0] - inTelescopePosition[0];
	double y = resultPosition[1] - inTelescopePosition[1];
	double z = resultPosition[2] - inTelescopePosition[2];

	double r = atan2( y, x );
	double d = atan2( z, sqrt( x*x + y*y ) );

	double pi = 3.1415926535;
	
	if( r < 0 )
		r = 2.0*pi + r;

	*Ra = r;
	*Dec = d;
}
//-----------------------------------------------------------------------//
// ������������� ����������� �����������
//-----------------------------------------------------------------------//
double OrbitIntegration_goodInit( )
{
	// init
	IForce = new Force::InfluenceForce;
	IForce->Init_CPU();

	// for test rotation only
	//IForce->TestRotation();

	return 0;
}
//-----------------------------------------------------------------------//
// ������� ��������������
//-----------------------------------------------------------------------//
double OrbitIntegration_good( SatParamToPredict &sptr,double H0, double E0)
{
	// ����� ����������
	int n = 6;
	double date_end = sptr.data_end;
	double time_end = sptr.time_end;
	double date = sptr.data_start;
	double time = sptr.time_start;

	// ��������� ������ ���������
	double XS[6];
	for( int it = 0; it < 6; it++ )
		XS[it] = sptr.inX[it];

	double t, ajd0, delt0;
	IForce->set_time(date, time, &ajd0, &delt0, &t );
	double t_e = IForce->get_time( date_end, time_end, ajd0, delt0 );
	IForce->S_ajd0 = ajd0;
	IForce->S_delt0 = delt0;

	const double Satm = sptr.atm;
	const double Ssun = sptr.sun;
	IForce->SetSigmaAtm( Satm );
	IForce->SetSigmaSun( Ssun );

	double h = H0;
	double *p = new double[ 8*n+9];
	int *np = new int[13];

	int kp = 1;
	double *e = new double[6];
	for( int it = 0; it < 6; it++ )
		e[it] = E0;

	rk8cc( n, t, XS, t_e, h, p, np, kp, e, RightFxyzv );

	for( int it = 0; it < 6; it++ )
	{
		sptr.outX[it] = XS[it];
	};
	delete[] e;
	delete[]np;
	delete[] p;
	return 0;
}
//-----------------------------------------------------------------------//
double OrbitIntegration_approx( SatParamToPredict &sptr,double H0, double E0 )
{
	// ����� ����������
	int n = 6;
	double date_end = sptr.data_end;
	double time_end = sptr.time_end;
	double date = sptr.data_start;
	double time = sptr.time_start;

	// ��������� ������ ���������
	double XS[6];
	for( int it = 0; it < 6; it++ )
		XS[it] = sptr.inX[it];

	double t, ajd0, delt0;
	IForce->set_time(date, time, &ajd0, &delt0, &t );
	double t_e = IForce->get_time( date_end, time_end, ajd0, delt0 );
	IForce->S_ajd0 = ajd0;
	IForce->S_delt0 = delt0;

	const double Satm = sptr.atm;
	const double Ssun = sptr.sun;
	IForce->SetSigmaAtm( Satm );
	IForce->SetSigmaSun( Ssun );

	double h = H0;
	double *p = new double[ 8*n+9];
	int *np = new int[13];

	int kp = 1;
	double *e = new double[6];
	for( int it = 0; it < 6; it++ )
		e[it] = E0;

	rk8cc( n, t, XS, t_e, h, p, np, kp, e, RightFxyzvAp);

	for( int it = 0; it < 6; it++ )
	{
		sptr.outX[it] = XS[it];
	};
	delete[] e;
	delete[]np;
	delete[] p;
	return 0;
}
