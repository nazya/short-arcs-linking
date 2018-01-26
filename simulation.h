Doub residual(Vec vec, Observs track) {
	Vec b;
	Doub resRa, resDec;
	Doub res = 0.;
	for (Int i = 0; i < track.obs.size(); i++) {
		vec = preciseExt(vec, track.obs.at(i).jd);
		b = track.sen.toObs(vec);
		resRa = b.vec[0] - track.obs.at(i).vec[0];
		resDec = b.vec[1] - track.obs.at(i).vec[1];
		res += resRa * resRa / track.obs.at(i).cov[0][0] + resDec * resDec / track.obs.at(i).cov[1][1];
	}
	return sqrt(res / 2. / track.obs.size());
}

Observs modelTrack(Vec orbit, Sensor tel, Int isAddNoise = 1) {
	Int nBlips = 10;
	Doub samplePeriod = 1.     /60./60./24.; 
	Observs observs(tel);
	
	//noise
	VecDoub ev(2, 0.);
	MatDoub cov(2, 2, 0.);
	cov[0][0] = 2.35e-13;
	cov[1][1] = 2.35e-13;
	Ullong rand = 14898;
	//Ullong rand = 14898;
	Multinormaldev blipNoise(rand, ev, cov);
	VecDoub noise;
	
	Est trackBlip(2);
	for (Int i = 0; i < nBlips; i++) {
		trackBlip = tel.toObs(orbit);

		if ( i==0 || i==1) {
			p2d.point(trackBlip.vec[0], trackBlip.vec[1], 3);
			
		}

		if (isAddNoise) {
			noise = blipNoise.dev();
			trackBlip.vec[0] += noise[0];
			trackBlip.vec[1] += noise[1];
		}

		trackBlip.cov = cov;

		if ( i==0 || i==1) {
			p2d.point(trackBlip.vec[0], trackBlip.vec[1], 5);
		}

		observs.obs.push_back(trackBlip);
		orbit = preciseExt(orbit, orbit.jd + samplePeriod);
	}
	return observs;
}

/*
Vec findSecondRun(Telescope& tel, Vec& orbit) {
	Est refBlip, curBlip;
	//Vec vec2 = predict(orbit, orbit.jDate + 1.);
	refBlip = toBlip(orbit);

	Vec vec2 = preciseExt(orbit, orbit.jDate + 2.);

	do { 
		vec2 = preciseExt(vec2, vec2.jDate + 1. / 24.);
		curBlip = toBlip(vec2);
	} while (0);//blipsRes(refBlip, curBlip) > 0.05);
	return vec2;
}

*/