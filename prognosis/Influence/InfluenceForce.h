//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force

// files ftp://maia.usno.navy.mil/ser7
//==============================================================================//
#ifndef _InfluenceForce_H_
#define _InfluenceForce_H_

#include "DefineParam.h"

// name
namespace Force
{
	class InfluenceForce
	{
	private:
	
		//================================================//
		// InfluencePlanet.cpp
		// ������ ������
		int DE403_clu;
		int DE403_step;
		int DE403_size;
		double *FileMemDE403;
		// �������� ��������
		void DE403LoadFileToMemory();
		// ������� ������
		void DE403DeleteMemory();
		// Planet Position
		double PLCOORD[11*3];

		// ���������� �������
		void pln_coords( int i, double *pos );
		// ����������� ��� �������
		double mu_plan( int i );
		// ����� ������� ���������
		bool pln_flag( int i );
		// ������� ��������
		void cheb2( int mode, double x, int degree, double *coeff, double *w, double *dw, double *ddw );
		// ���������� ������� �� �������� ��������
		void DE403( double t, int n_pl, int mode, double* x_pl, double ajd0, double delt0 );
		// ���������� ��������� ������ ������������ �����
		void planets_update_geo( double t, double ajd0, double delt0 );
		// ���������� ��������� ������ ������������ ������
		void planets_update_sol( double t, double ajd0, double delt0 );
		// ���������� �� ������
		void planets_grav( double *x, double *f_gr );
		//================================================//

		//================================================//
		//InfluenceSun.cpp
		// ���������� ������ ������
		double shadow_geo( double *r_ka, double *r_sun );
		// ����
		void sp_cannonball( double *x, double *f_sp, double *pln_coords11, double mu_plan11, double sp_q );
		// �������� ������
		void sp_cannonballForce( double *x, double sp_q, double *f_sp );
		//================================================//

		//================================================//
		//InfluenceEGM96.cpp
		double *h_egm96;

		// ��������� �����
		void GetF_Harm_egm96( double *x_b, int n_harm, double *f_harm );
		// ���� ��-�� ��������
		void GetHarmForce( double *x, double *Fharm );
		// init egm
		void InitHarmInCPU();
		//================================================//

		//================================================//
		// InfluenceAtmosphere.cpp
		// ��������� �����
		void Atm_drag( double *x, double t, double *f, double sigma_up, double ajd0, double delt0 );
		// ���� ���������
		//double Roa2004_2( double time, double *x, double ajd0, double delt0 );
		//================================================//

		//================================================//
		// InfluenceTime.cpp
		// �������� �������, ������� ������� � ������ �������
		// ������ � ���������� �� �����
		int TAUUTC_size;
		int TAUUTC_kmax;
		double *TAIUTCCorrect;
		
		// ������� ������� � ���������� �� �����
		int InitTAU_UTCcorrect();
		// �������� �������
		void DeInitTAU_UTCcorrect();
		// ��������� �������� ������� - ����� �����
		int GetTAU_UTCkmax();
		// ������ �������
		int GetTAU_UTC_size();
		//��������� ����� TDB � UTC
		double tdb_utc( double ajd, double delt, double t );
		// YYYYMMDD.0d0 � HHMMSS.SSSS... � ������� �������������
		void jddeltt( double dt, double tm, double *ajd, double *delt, double *t );
		//================================================//

		//================================================//
		//InfluenceNutationEarth.cpp
		// ��������� � ������� �����
		//������ �������� ������� ���������
		double E2000( double E1, double E2);
		void FA2000( double AED, double *FA );
		//���������� �������� �������
		void N2000( int N, double AJD, double *HYT );
		//����� ������� �������
		void NM2000( double E, double *HUT );
		//================================================//

		//================================================//
		// InfluiencePrecessionEarth.cpp
		//������ ������� ��������� 
		void PM2000( double E1, double E2, double *P );
		//================================================//

		//================================================//
		// InfluencePoleEarth.cpp
		// ������� ��������
		int finals_n;
		double *finals_tab;
	
		// ������������� ������
		void iers_init();
		// �������� ��������
		void iers_delete();
		// ��������� ��������� �������� ��������� ������ ����� � ������
		void get_xyt ( double t, double *xyt, double ajd0, double delt0 );
		// ������� �������� ��������������� �������� �������
		void PM_mat ( double x, double y, double *A_pole );
		//================================================//

		//================================================//
		// InfluenceEarthRotation.cpp
		// ��������� �������� ������� ��������� �� EME2000 � ����������� ��� �������� ��
		void state_to_itrf( double t, double *x_2000, double *x_g, double *A_rot );
		// ��������� ��������� ������� ��������������� ��������� �������� �����
		void ER_mat ( double t, double d_UT1, double *A_rotat, double ajd0, double delt0 );
		// ��������� ��������� ������� �������� �� ITRF � ICRF
		void iers_mat( double t, double *A, double ajd0, double delt0 );
		//================================================//

		//================================================//
		// AdditionalFunction.cpp
		// �������������� �������
		void matMul( double *in1, double *in2, double *out );
		void matVecMul_V6( double *inMat, double *inVec, double *outVec );
		double DMOD( double X, double Y );
		double DDIM( double X, double Y  );
		double DSIGN( double X, double Y );
		double Dmax( double X, double Y );
		double Dmin( double X, double Y );
		//================================================//

	public:
		InfluenceForce();
		~InfluenceForce();

		// ������������
		int ST;
		double SIGMA_ATM;
		double SIGMA_SUN;
		// �������� �������
		double S_ajd0;
		double S_delt0;

		// ������������� ������ �� CPU
		void Init_CPU();
		void DeInit_CPU();
		// Init koeff
		void SetSigmaAtm( double sg ) { SIGMA_ATM = sg; };
		void SetSigmaSun( double sg ) { SIGMA_SUN = sg; };

		//================================================//
		// InfluenceTime.cpp
		void set_time( double date, double time, double *ajd0, double *delt0, double *t );
		double get_time( double date, double time, double ajd0, double delt0 );
		double Ajd_dt ( double ajd);
		double dt_ajd ( double dt );
		// �� �������� ������������� � YYYYMMDD.0d0 � HHMMSS.SSSS
		void dttm( double ajd, double delt, double t, double *dt, double *tm );
		//================================================//

		//================================================//
		// InfluenceEarthRotation.cpp
		void GetTemeMatrix( double t, double *A, double ajd0, double delt0 );
		void iers_update_matrix( double t, double *A_rot, double ajd0, double delt0 );
		//================================================//

		//================================================//
		// AdditionalFunction.cpp
		void matVecMul( double *inMat, double *inVec, double *outVec );
		void transpose( double *A, double *B );
		//================================================//

		//================================================//
		// InfluenceForce.cpp
		// ���������� ������ ������
		int rh_fun( double t, double *x, double *f );
		int rh_fun_kepler(double t, double *xt, double *fx );
		int rh_fun_grav( double t, double *x, double *f );
		//================================================//

		//================================================//
		// Test
		int TestRotation();
		//================================================//
	};
};
//==============================================================================//
#endif