static SatParamToPredict satParamToPredict;

void prognosisInit() {
	OrbitIntegration_goodInit();
}


Vec keplerExt(Vec vecToPredict, Doub endDate) {
	//SatParamToPredict satParamToPredict;
	Date dateStart(vecToPredict.jd);
	Date dateEnd(endDate);
	satParamToPredict.data_start = dateStart.day;
	satParamToPredict.time_start = dateStart.time;
	satParamToPredict.data_end = dateEnd.day;
	satParamToPredict.time_end =  dateEnd.time;
	for (int i = 0; i < 6; i++) {
		satParamToPredict.inX[i] = vecToPredict.vec[i];
	}
	satParamToPredict.atm = 0.3E-2;
    satParamToPredict.sun = 0.5E-05;
	OrbitIntegration_approx(satParamToPredict);
	vecToPredict.jd = endDate;
	for (int i = 0; i < 6; i++) {
		vecToPredict.vec[i] = satParamToPredict.outX[i];
	}
	return vecToPredict;
}

Vec preciseExt(Vec vecToPredict, Doub endDate) {
	//SatParamToPredict satParamToPredict;
	Date dateStart(vecToPredict.jd);
	Date dateEnd(endDate);
	satParamToPredict.data_start = dateStart.day;
	satParamToPredict.time_start = dateStart.time;
	satParamToPredict.data_end = dateEnd.day;
	satParamToPredict.time_end =  dateEnd.time;
	for (int i = 0; i < 6; i++) {
		satParamToPredict.inX[i] = vecToPredict.vec[i];
	}
	satParamToPredict.atm = 0.3E-2;
    satParamToPredict.sun = 0.5E-05;
	OrbitIntegration_good(satParamToPredict);
	vecToPredict.jd = endDate;
	for (int i = 0; i < 6; i++) {
		vecToPredict.vec[i] = satParamToPredict.outX[i];
	}
	return vecToPredict;
}

Vec rflPreciseExt(Vec vecToPredict, Doub endDate) {
	//SatParamToPredict satParamToPredict;
	Date dateStart(vecToPredict.jd);
	Date dateEnd(endDate);
	satParamToPredict.data_start = dateStart.day;
	satParamToPredict.time_start = dateStart.time;
	satParamToPredict.data_end = dateEnd.day;
	satParamToPredict.time_end =  dateEnd.time;

	Doub r = vecToPredict.vec[0];
	VecDoub xyz(3);
	Doub ra = vecToPredict.vec[1];
	Doub dec = vecToPredict.vec[2];
	
	xyz[0] = r * cos(ra) * cos(dec);
	xyz[1] = r * sin(ra) * cos(dec);
	xyz[2] = r * sin(dec);

	

	for (int i = 0; i < 3; i++) {
		satParamToPredict.inX[i] = xyz[i];
		satParamToPredict.inX[i+3] = vecToPredict.vec[i+3];
	}
	satParamToPredict.atm = 0.3E-2;
    satParamToPredict.sun = 0.5E-05;
	OrbitIntegration_good(satParamToPredict);
	vecToPredict.jd = endDate;

	r = 0;
	for (int i = 0; i < 3; i++) {
		r += satParamToPredict.outX[i] * satParamToPredict.outX[i];
		vecToPredict.vec[i+3] = satParamToPredict.outX[i+3];
	}
	r = sqrt(r);
	vecToPredict.vec[0] = r;
	ra = atan2(satParamToPredict.outX[1], satParamToPredict.outX[0]);
	if( ra < 0 )
		ra = 2.0*3.1415926535 + ra;


	vecToPredict.vec[1] = ra;
	vecToPredict.vec[2] = atan2(satParamToPredict.outX[2], sqrt(satParamToPredict.outX[0]*satParamToPredict.outX[0] + satParamToPredict.outX[1]*satParamToPredict.outX[1]));


	return vecToPredict;
}