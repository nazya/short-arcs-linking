//==============================================================================//
// Andrianov N.G.
// opbit predict 
// module find Influence Force
// main
//==============================================================================//
#include <math.h>
#include <stdio.h>

#include "InfluenceForce.h"

namespace Force
{
	//==============================================================================//
	// create and delete
	//==============================================================================//
	InfluenceForce::InfluenceForce()
	{
	};
	InfluenceForce::~InfluenceForce()
	{
		
	}
	//==============================================================================//
	// задание массивов с константами
	// загрузка эфемерид планет
	// загрузка таблис с веременем
	// копирование константных коэффициентов на GPU
	// !!! массивы с коэффициентами повторяются!! важно учесть при изменениях
	//==============================================================================//
	void InfluenceForce::Init_CPU()
	{
		printf("Init Influence Force module CPU.... \n");
		// установка поправок для TAI UTC
		InitTAU_UTCcorrect();
		// загрузка эфемерид
		DE403LoadFileToMemory();
		// иницилизация модуля движения земли
		iers_init();
		// гармоники
		InitHarmInCPU();
		ST = 0;

		printf("OK\n");
	}
	//==============================================================================//
	// Удаление ресурсов
	//==============================================================================//
	void InfluenceForce::DeInit_CPU()
	{
		printf( "Influence delete CPU\n" );
		
		DE403DeleteMemory();
		DeInitTAU_UTCcorrect();
		iers_delete();
		delete h_egm96;

		printf("OK\n");
	}
	//==============================================================================//
	// Процедура расчета правых частей уравнений движения + вариационных уравнений.
	//==============================================================================//
	int InfluenceForce::rh_fun( double t, double *x, double *f )
	{
		ST++;
		double x_g[6];
		double f_hrm[3];
		double f_atm[3];
		double f_sp[3];
		double f_ah[3];
		double f_gr[3];
		double  jd_t, r_abs;
		double A_rot[9];

		// матрица поворота
		iers_update_matrix( t, A_rot, S_ajd0, S_delt0 );
		// положение планет
		planets_update_geo( t, S_ajd0, S_delt0 );

		for( int it = 0; it < 6; it++ )
			f[it] = 0.0;

		f[0] = x[3];
		f[1] = x[4];
		f[2] = x[5];

		// юлианская дата
		jd_t = S_ajd0 + ( S_delt0 + t )/86.40;
		// радиус вектор
		r_abs = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
		// поворот вектора в систему координат чемли
		state_to_itrf( t, x, x_g, A_rot );

		// влияние гармоник
		//GetF_Harm_egm96( x_g, 36, f_hrm );
		GetF_Harm_egm96( x_g, 75, f_hrm );

		// влияние атмосферы
		double sigma_up = SIGMA_ATM;// 0.3E-2;
		Atm_drag( x_g, t, f_atm, sigma_up, S_ajd0, S_delt0 );
		//f_atm[0] = 0.0;
		//f_atm[1] = 0.0;
		//f_atm[2] = 0.0;
		
		//FILE *fre;
		//fre = fopen( "fatm.txt", "at" );
		//fprintf( fre, "%.12e %.12e %.12e %e\n", f_atm[0], f_atm[1], f_atm[2], t );
		//// 1 1.12523154311973767E-010 -5.93712291122998987E-010 -5.89755768296451942E-010
		//fclose( fre);
		
		// солнечное давление
		double sp_q = SIGMA_SUN; //0.5E-05;
		sp_cannonballForce( x, sp_q, f_sp );
		//f_sp[0] = 0.0;
		//f_sp[1] = 0.0;
		//f_sp[2] = 0.0;

		// гравитация планет
		planets_grav( x, f_gr );

		// перенос вектора воздействия гармоник и атмосферы
		double invA[9];
		transpose( A_rot, invA ); 
		f_hrm[0] = f_atm[0] + f_hrm[0];
		f_hrm[1] = f_atm[1] + f_hrm[1];
		f_hrm[2] = f_atm[2] + f_hrm[2];
		matVecMul( invA, f_hrm, f_ah );

		// суммирование воздействий
		f[3] = f[3] + f_ah[0] + f_gr[0] + f_sp[0];
		f[4] = f[4] + f_ah[1] + f_gr[1] + f_sp[1];
		f[5] = f[5] + f_ah[2] + f_gr[2] + f_sp[2];

		return 0;
	}
	//==============================================================================//
	// функция в парой части для кеплеровского движения
	//==============================================================================//
	int InfluenceForce::rh_fun_kepler(double t, double *xt, double *fx )
	{
		double EARTHR = 6370000.0;
		double ACCELERATION = 9.822;
		double mm = ACCELERATION*EARTHR*EARTHR;

		double x = xt[0];
		double y = xt[1];
		double z = xt[2];

		double vx = xt[3];
		double vy = xt[4];
		double vz = xt[5];

		double r = 1000000*sqrt( x*x + y*y + z*z );
		double R3 = r*r*r;

		double ffx = vx;
		double ffy = vy;
		double ffz = vz;
		double ffvx = -x*mm/R3;
		double ffvy = -y*mm/R3;
		double ffvz = -z*mm/R3;

		fx[0] = ffx;
		fx[1] = ffy;
		fx[2] = ffz;

		fx[3] = ffvx;
		fx[4] = ffvy;
		fx[5] = ffvz;

		return 0;
	};
	//==============================================================================//
	// учет только гравитации, гармоник, давления солнца
	//==============================================================================//
	int InfluenceForce::rh_fun_grav( double t, double *x, double *f )
	{
		double f_gr[3];
		double A_rot[9];

		double x_g[6];
		double f_hrm[3];
		double f_ah[3];
		double f_sp[3];

		// матрица поворота
		iers_update_matrix( t, A_rot, S_ajd0, S_delt0 );
	
		//if( Nam < 6 ){
		//	for( int j = 0; j < 9; j++ )
		//		printf( "%.10f ", A_rot[j] );
		//	printf( "%.10e %.10e %.10e\n", t, S_ajd0, S_delt0);
		//}
	
		// положение планет
		planets_update_geo( t, S_ajd0, S_delt0 );

		for( int it = 0; it < 6; it++ )
			f[it] = 0.0;

		f[0] = x[3];
		f[1] = x[4];
		f[2] = x[5];

		// гравитация планет
		planets_grav( x, f_gr );
		
		// поворот вектора в систему координат чемли
		state_to_itrf( t, x, x_g, A_rot );

		// влияние гармоник
		GetF_Harm_egm96( x_g, 75, f_hrm );

		// солнечное давление
		double sp_q = 0.5E-05;
		sp_cannonballForce( x, sp_q, f_sp );

		// перенос вектора воздействия гармоник и атмосферы
		double invA[9];
		transpose( A_rot, invA ); 
		matVecMul( invA, f_hrm, f_ah );
		
		//if( Nam < 6 )
		//printf( "%.10e %.10e %.10e\n", f_ah[0], f_ah[1], f_ah[2] );
		
		// суммирование воздействий
		f[3] = f[3] + f_ah[0] + f_gr[0] + f_sp[0];
		f[4] = f[4] + f_ah[1] + f_gr[1] + f_sp[1];
		f[5] = f[5] + f_ah[2] + f_gr[2] + f_sp[2];

		// суммирование воздействий
		//f[3] = f_gr[0];
		//f[4] = f_gr[1];
		//f[5] = f_gr[2];

		//Nam++;
		return 0;
	};
	//==============================================================================//
};