struct STrack {
	Observs refTrack;
	VecDoub raabc, decabc;
	Doub tscale;


	STrack(Observs track) : refTrack(track) {
		tscale = 10000. / 1.;//60. * 60. * 24.;
		raabc = VecDoub(3);
		decabc = VecDoub(3);
		smooth(track);
	}


	Vec blip(Doub jd) {
		Doub _t = ti(jd);
		Vec b(2);
		b.jd = jd;
		b.vec[0] = raabc[0]*_t*_t + raabc[1]*_t + raabc[2];
		b.vec[1] = decabc[0]*_t*_t + decabc[1]*_t + decabc[2];
		return b;
	}

	Vec dBlip(Doub jd) {
		Doub _t = ti(jd);
		Vec b(2);
		b.jd = jd;
		b.vec[0] = 2.*raabc[0]*_t * tscale + raabc[1] * tscale;
		b.vec[0] *= 1. / 60. / 60. / 24.;

		b.vec[1] = 2.*decabc[0]*_t*tscale + decabc[1] * tscale;
		b.vec[1] *= 1. / 60. / 60. / 24.;

		return b;
	}

	Vec ddBlip(Doub jd) {
		Vec b(2);
		b.jd = jd;
		b.vec[0] = 2.*raabc[0] * tscale * tscale;
		b.vec[0] *= 1. / 60. / 60. / 24.;
		b.vec[0] *= 1. / 60. / 60. / 24.;
		b.vec[1] = 2.*decabc[0] * tscale * tscale;
		b.vec[1] *= 1. / 60. / 60. / 24.;
		b.vec[1] *= 1. / 60. / 60. / 24.;

		return b;
	}
	
	Doub ti(Doub jd) {
		Doub t = jd - refTrack.obs[0].jd;// - (track.obs.back().jd - track.obs.front().jd) / 2.;
		return (t * tscale);
	}


	void smooth(Observs track) {
		Doub t, _t;
		MatDoub a(3, 3, 0.);
		VecDoub raright(3, 0.), decright(3, 0.);
		for (int i = 0; i < track.obs.size(); i++) {
			raright[2] += track.obs[i].vec[0];
			decright[2] += track.obs[i].vec[1];

			t = _t = (track.obs[i].jd - track.obs[0].jd) * tscale;
		
			raright[1] += track.obs[i].vec[0] * _t;
			decright[1] += track.obs[i].vec[1] * _t;
			a[2][1] += _t;
		
			_t *= t; //t2
			raright[0] += track.obs[i].vec[0] * _t;
			decright[0] += track.obs[i].vec[1] * _t;
			a[2][0] += _t;

			_t *= t; //t3
			a[1][0] += _t;
			
			_t *= t; //t4
			a[0][0] += _t;
			
		}
		a[1][2] = a[2][1];

		a[0][2] = a[2][0];
		a[1][1] = a[2][0];

		a[0][1] = a[1][0];

		a[2][2] = Doub(track.obs.size());

		SVD svd(a);
		svd.solve(raright, raabc);
		svd.solve(decright, decabc);
		
		for (int i = 0; i < track.obs.size(); i++) {
			//p2d.point(track.obs[i].jd, track.obs[i].vec[0], 3);
			//p2d.point(track.obs[i].jd, blip(track.obs[i].jd).vec[0]);
		}
		//p2d.show();

	}
};
 

void estim(STrack track) {
	int m = 0;
	
	Vec b = track.blip(track.refTrack.obs[m].jd);
	Doub ra = b.vec[0], dec = b.vec[1];
	b = track.dBlip(track.refTrack.obs[m].jd);
	print(b.vec[0]);
	Doub dradt = b.vec[0], ddecdt = b.vec[1];
	b = track.ddBlip(track.refTrack.obs[m].jd);
	Doub d2radt2 = b.vec[1], d2decdt2 = b.vec[1];



	Doub cosra = cos(ra), sinra = sin(ra), cosdec = cos(dec), sindec = sin(dec);


	VecDoub e(3), dedt(3), d2edt2(3);
	e[0] = cosra*cosdec; e[1] = sinra*cosdec; e[2] = sindec;
	dedt[0] = -sinra*cosdec*dradt - cosra*sindec*ddecdt; dedt[1] = cosra*cosdec*dradt - sinra*sindec*ddecdt; dedt[2] = cosdec*ddecdt;
	d2edt2[0] = -cosra*cosdec*dradt*dradt + sinra*sindec*dradt*ddecdt - sinra*cosdec*d2radt2 + sinra*sindec*ddecdt*dradt - cosra*cosdec*ddecdt*ddecdt - cosra*sindec*d2decdt2;
	d2edt2[1] = -sinra*cosdec*dradt*dradt - cosra*sindec*dradt*ddecdt + cosra*cosdec*d2radt2 - cosra*sindec*dradt*ddecdt - sinra*cosdec*ddecdt*ddecdt - sinra*sindec*d2decdt2;
	d2edt2[2] = -sindec*ddecdt + cosdec*d2decdt2;

	Doub dt = 1.     /60./60./24.;
	VecDoub s = track.refTrack.sen.icrfPos(track.refTrack.obs[m].jd), dsdt(3), d2sdt2(3);
	VecDoub _s = track.refTrack.sen.icrfPos(track.refTrack.obs[m].jd - dt);
	VecDoub s_ = track.refTrack.sen.icrfPos(track.refTrack.obs[m].jd + dt);
	for (int i = 0; i < 3; i++) {
		dsdt[i] = (s_[i] - _s[i]) / 2.;// / (2. * dt);
	}
	//_s = track.refTrack.sen.icrfPos(track.refTrack.obs[m].jd - 2*dt);
	//s_ = track.refTrack.sen.icrfPos(track.refTrack.obs[m].jd + 2*dt);
	//dt = dt * 100.;
	for (int i = 0; i < 3; i++) {
		d2sdt2[i] = (s_[i] - 2 * s[i] + _s[i]);// / (1. * dt * dt);
	}

	print(d2sdt2);

	VecDoub out(3);
	VecDoub right(3);
	Doub d;
	Doub k = 398.6004415;

	out[0] = dOrb(vec1, track.refTrack.sen);
	d =  out[0];

	while(1) {
	d =  out[0];


	Doub r = sqrt(innerProduct(s, s) + d*d*innerProduct(e, e) + 2*d*innerProduct(s, e));

	print("R");
	print(r);

	for (int i = 0; i < 3; i++) {
		right[i] = -(s[i] + d * e[i]) * k  /r/r/r - d2sdt2[i];
	}

	MatDoub a(3, 3);
	for (int i = 0; i < 3; i++) {
		a[i][0] = d2edt2[i];
		a[i][1] = 2*dedt[i];
		a[i][2] = e[i];
	}

	SVD svd(a);

	

	svd.solve(right, out);

	VecDoub v(3);
	for (int i = 0; i < 3; i++) {
		v[i] = out[1] * e[i] + out[0] * dedt[i] + dsdt[i];
	}

	print("v");
	print(v);
	print(out[0]);

	
	}
}


Vec prelim(Doub d, STrack strack) {
	VecDoub icrfTel = strack.refTrack.sen.icrfPos(strack.refTrack.obs.front().jd);
	Doub dt = 0.;
	
	Vec blip = strack.blip(strack.refTrack.obs.front().jd);

	//VecDoub x0 = radec2xyz(blip.vec[0], blip.vec[1]);
	VecDoub x0 = radec2xyz(strack.refTrack.obs.front().vec[0], strack.refTrack.obs.front().vec[1]);
	for (int i = 0; i < 3; i++) {
		x0[i] = x0[i] * d + icrfTel[i];
	}

	VecDoub x(3);
	VecDoub e(3);
	VecDoub xs(3);
	VecDoub y(3, 0.);
	MatDoub p;
	VecDoub pXx;
	MatDoub a(3, 3, 0.);
	VecDoub g = gravityVector(x0);

	Int it = 0;
	while ( it < strack.refTrack.obs.size() && dt < 31.) {
		dt = 60. * 60. * 24. * (strack.refTrack.obs.at(it).jd - strack.refTrack.obs.front().jd);
		xs = strack.refTrack.sen.icrfPos(strack.refTrack.obs.at(it).jd);
		for (Int i = 0; i < 3; i++) {
			x[i] = -dt * (x0[i] - xs[i] + 0.5 * dt * dt * g[i]);
		}
		e = radec2xyz(strack.refTrack.obs.at(it).vec[0], strack.refTrack.obs.at(it).vec[1]);
		p = projectionMatrix(e);
		pXx = matXVec(p, x);
		for (Int i = 0; i < 3; i++) {
			y[i] += pXx[i];
		}
		for (Int i = 0; i < a.nrows(); i++) {
			for (Int j = 0; j < a.ncols(); j++) {
				a[i][j] += dt * dt * p[i][j];
			}
		}
		it++;
	}

	SVD svd(a);
	VecDoub v(3);
	svd.solve(y, v);
	
	Vec orb;
	for (Int i = 0; i < 3; i++) {
		orb.vec[i] = x0[i];
	}
	for (Int i = 3; i < 6; i++) {
 		orb.vec[i] = v[i - 3] * 1000;
	}
	orb.jd = strack.refTrack.obs.front().jd;
	return orb;
}

/*
Est estimate(STrack st) {
	for (Doub d = 0.1; d < 80.; d+=1) {
		Est est = prelim(d, st);

	}
}
*/
