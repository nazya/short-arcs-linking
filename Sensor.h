struct Sensor {
	VecDoub pos;

	Sensor () {
	};

	Sensor(VecDoub pos) : pos(pos) {
	}

	VecDoub icrfPos(Doub jd) {
		Date date(jd);
		Doub telICRF[3];
		VecDoub ICRFPos(3);
		GetTELPos(date.day, date.time, &pos[0], telICRF);
		for (int i = 0; i < 3; i++) {
			ICRFPos[i] = telICRF[i];
		}
		return ICRFPos;
	}

	Vec toObs(Vec& orb) {
		Date date(orb.jd);
		Doub telICRF[3];
		Vec obs(2);
		GetTELPos(date.day, date.time, &pos[0], telICRF);
		ConvertXYZtoRADEC(&orb.vec[0], telICRF, &obs.vec[0], &obs.vec[1]);
		obs.jd = orb.jd;
		return obs;
	}

	Vec toObsRFLVec(Vec& rfl) {
		Vec orb = rfl2xyzVec(rfl);
		Date date(orb.jd);
		Doub telICRF[3];
		Vec obs(2);
		GetTELPos(date.day, date.time, &pos[0], telICRF);
		ConvertXYZtoRADEC(&orb.vec[0], telICRF, &obs.vec[0], &obs.vec[1]);
		obs.jd = orb.jd;
		return obs;
	}
};