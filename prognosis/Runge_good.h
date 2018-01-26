//-----------------------------------------------------------------------//
// Андрианов Н.Г.
// Интегрирование уравнения движения
//-----------------------------------------------------------------------//
#ifndef SPARAM
#define SPARAM

struct SatParamToPredict
{
	double inX[6];
	double sun;
	double atm;

	double data_start;
	double time_start;

	double outX[6];
	double data_end;
	double time_end;
};


double OrbitIntegration_goodInit( );
double OrbitIntegration_good( SatParamToPredict &sptr,double H0 = 5.0E-5, double E0 = 1.0E-11 );
double OrbitIntegration_approx( SatParamToPredict &sptr,double H0 = 5.0E-5, double E0 = 1.0E-11);

void  GetTELPos( double date1, double time1, double *Tpos, double *Telicrf );
void ConvertXYZtoRADEC( double *resultPosition, double *inTelescopePosition, double *Ra, double *Dec );

#endif