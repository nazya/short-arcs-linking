typedef struct Date {
	Doub day;
	Doub time;
	Date(){
	}
	Date(Doub day, Doub time) : day(day), time(time) {
	}
	Date(Doub jDate) {
		jDateToDate(jDate, this);
	}


	void jDateToDate(Doub &jDate, Date *date) {
		date->day = JDtoYYYYMMDD(jDate);
		date->time = SECtoHHMMSS(date->day, jDate);
	}

	void dateToJDate(Date &date, Doub *jDate) {
		*jDate = YYYYMMDDtoJD(*jDate);
		*jDate += date.time;
	}


	int NINTC( double A ) {
		int ia = (int)A;
		double da = ia;

		if( A > 0 )
		{
			int delta = (int)(10.0*(A-da));
			if( delta >= 5 )
				ia = ia+1;
		}
		else
		{
			int delta = (int)(10.0*(A-da));
			if( delta <= -5 )
				ia = ia-1;
		}
		return ia;
	}
	//==============================================================================//
	// ������� ���������� ����, ���������� � ������� YYYYMMDD.0000d0, 
	// �� ��������� ����.
	//==============================================================================//
	double JDtoYYYYMMDD( double ajd )
	{
		int y, d, j, m;

		j = NINTC(ajd-1721119.0);
		y = (4*j-1)/146097;
		d = (4*j-1-146097*y)/4;
		j = (4*d+3)/1461;
		d = (4*d+7-1461*j)/4;
		m = (5*d-3)/153;
		d = (5*d+2-153*m)/5;
		y = 100*y+j;
		if (m < 10)
		{
			m=m+3;
		}
		else
		{
			m=m-9;
			y=y+1;
		}
		double ajd_dt = y*10000+m*100+d;
		return ajd_dt;
	}
	//==============================================================================//
	// ������� ��������� ��������� ���� �� ����, ���������� � �������
	// YYYYMMDD.0000d0, �.�. ������� ���� �� �������� ����� ������.
	// ���������� ��������� ���� �������� ������� ������.
	//==============================================================================//
	double YYYYMMDDtoJD( double dt )
	{
		double dt_ajd;

		int year, month, day, century, year_c;
		year = floor(dt/10000.0);
		month = floor((dt-year*10000.0)/100.0);
		day = NINTC(dt) - year*10000 - month*100;
		if (month > 2)
		{
			month = month - 3;
		}
		else
		{
			month = month + 9;
			year = year - 1;
		}

		century = year/100;
		year_c = year - 100*century;
		dt_ajd = (146097*century)/4+(1461*year_c)/4+ (153*month+2)/5+day+1721119;

		return  dt_ajd;
	} 
	//==============================================================================//
	// ��������� �� ������� ��������� ���� HHMMSS.SSSS
	//==============================================================================//
	double SECtoHHMMSS( double data, double jd )
	{
		// ��������� ���� ������ �����
		double jdstart = YYYYMMDDtoJD( data ) - 0.5;	

		// ����� ������ �� ������ �����
		double Day = jd - jdstart;
		// ����� ����
		int hour = Day*86400.0/3600.0;	
		// ����� �����
		int minute = ( Day*86400.0 - ((double)hour)*3600.0 )/60.0;	
		// ������� ����� ������, �������
		double sec = ( Day*86400.0 - ((double)hour)*3600.0  - (double)minute*60.0);
		// ����� � �������  HHMMSS.SSSS 
		double d2 = (double)hour*10000.0 + (double)minute*100.0 + sec;

		return d2;
	}
};