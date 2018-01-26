Doub orbitEnergy(Vec orbit) {
	Doub v2 = 0.;		//velocity
	Doub r = 0.;		//distance
	Doub k = 398600.4415;	//gravity parameter km3 / s2
	//k *= 60. * 60.;


	for (Int i = 0; i < 3; i++) {
		r += orbit.vec[i] * orbit.vec[i];
		v2 += orbit.vec[i + 3] * orbit.vec[i + 3];
	}
	r = sqrt(r);
	return v2 - 2. * k / r / 1000.;

}

VecDoub orbitSigma(Vec orbit) {
	VecDoub pos(3);
	VecDoub vel(3);
	for (Int i = 0; i < 3; i++) {
		pos[i] = orbit.vec[i] * 1000.;
		vel[i] = orbit.vec[i + 3];
	}
	VecDoub sigma = crossProduct(pos, vel);
	return sigma;
}

Doub isOrbit2(Vec orbit) {
	Doub h = orbitEnergy(orbit);				//energy
	if (true) {
		if (h >= -398600.4415 / 70000.) {
			printf("not eleptic orbit\n");
			return 0;
		}
	}
	if (true) {
		VecDoub pos(3);
		VecDoub vel(3);
		for (Int i = 0; i < 3; i++) {
			pos[i] = orbit.vec[i] * 1000.;
			vel[i] = orbit.vec[i + 3];
		}
		VecDoub sigma = crossProduct(pos, vel);
		VecDoub lambda = crossProduct(vel, sigma);
		VecDoub re = normilize(pos);
		for (Int i = 0; i < 3; i++) {
			lambda[i] -= 398600.4415 * re[i];
		}
		Doub pericentre = (-398600.4415 + scalar(lambda)) / h;
		printf("per: %lf\n", pericentre);
		printf("apo: %lf\n", (-398600.4415 - scalar(lambda)) / h);
		if (pericentre <= 6370. + 50.) {
			printf("pericentre is too low\n");
			return 0;
		}
	}
	Doub b = -398600.4415 / h;
	Doub q = b / 6371.11;
	return 84.3/60./24. * q * sqrt(q);
	return 2. * 3.1415926535 / sqrt(398600.4415) * b * sqrt(b) /60./60./24.;
}


Doub isOrbit(Vec orbit) {
	Doub h = orbitEnergy(orbit);				//energy
	if (false) {
		if (h >= 0.) {
			printf("not eleptic orbit\n");
			return 0;
		}
	}
	if (true) {
		VecDoub pos(3);
		VecDoub vel(3);
		for (Int i = 0; i < 3; i++) {
			pos[i] = orbit.vec[i] * 1000.;
			vel[i] = orbit.vec[i + 3];
		}
		VecDoub sigma = crossProduct(pos, vel);

		Doub sigma2 = innerProduct(sigma, sigma);
		Doub eps = sqrt(1. + h * sigma2 / 398600.4415 / 398600.4415);
		Doub p = sigma2  / 398600.4415;
		Doub pericentre = p / (1. + eps);
			//printf("per: %lf\n", pericentre);
			//printf("apo: %lf\n", p / (1. - eps));

		if (pericentre < 6370. + 99.) {
			//printf("per: %lf\n", pericentre);
			//printf("pericentre is too low\n");
			return 0;
		}
		if (p / (1. - eps) > 100001.) {
			printf("per: %lf\n", p / (1. - eps));
			printf("apocentre is too high\n");
			return 0.;

		}
	}
	Doub b = -398600.4415 / h;
	Doub q = b / 6371.11;
	return 84.3/60./24. * q * sqrt(q);
	return 2. * 3.1415926535 / sqrt(398600.4415) * b * sqrt(b) /60./60./24.;
}



VecDoub gravityVector(VecDoub atPoint) {
	Doub g = 0.009819900504 / 1000.;
	Doub earthRadius = 6.37111;
	Doub r2 = 0.;
	for (Int i = 0; i < atPoint.size(); i++) {
		r2 += atPoint[i] * atPoint[i];
	}
	for (Int i = 0; i < atPoint.size(); i++) {
		atPoint[i] *= -g * earthRadius * earthRadius / r2 / sqrt(r2);
	}
	return atPoint;
}

VecDoub radec2xyz(Doub ra, Doub dec) {
	VecDoub xyz(3);
	xyz[0] = cos(ra) * cos(dec);
	xyz[1] = sin(ra) * cos(dec);
	xyz[2] = sin(dec);
	return xyz;
}



MatDoub radec2xyzJac(Vec& b) {
	MatDoub jac(3, 2);
	jac[0][0] = -sin(b.vec[0]) * cos(b.vec[1]);
	jac[1][0] = cos(b.vec[0]) * cos(b.vec[1]);
	jac[2][0] = sin(b.vec[1]);
	jac[0][1] = -cos(b.vec[0]) * sin(b.vec[1]);
	jac[1][1] = -sin(b.vec[0]) * sin(b.vec[1]);
	jac[2][1] = cos(b.vec[1]);
	return jac;
}

MatDoub covWithFixedD(Observs track, Doub d, Int startBlip = 0, Doub areaSize = 1.) {
	VecDoub icrfTel = track.sen.icrfPos(track.obs.at(startBlip).jd);
	Doub ra = track.obs.at(startBlip).vec[0], dec = track.obs.at(startBlip).vec[1];
	Doub cosra = cos(ra), sinra = sin(ra), cosdec = cos(dec), sindec = sin(dec);

	VecDoub x0(3), nra(3), ndec(3);
	x0[0] = cosra*cosdec;
	x0[1] = sinra*cosdec;
	x0[2] = sindec;
	for (int i = 0; i < 3; i++) {
		x0[i] = x0[i] * d + icrfTel[i];
	}
	x0 = normilize(x0);

	MatDoub cov(3, 3);
	cov[0][0] = x0[0];
	cov[1][0] = x0[1];
	cov[2][0] = x0[2];
	cov[0][1] = -sinra*cosdec;
	cov[1][1] = cosra*cosdec;
	cov[2][1] = 0;
	cov[0][2] = -cosra*sindec;
	cov[1][2] = -sinra*sindec;
	cov[2][2] = cosdec;

	MatDoub sig(3, 3, 0.);
	Doub zn = 1e-4;
	Doub nzn = d * sqrt(2*zn);
	sig[0][0] = zn;
	sig[1][1] = nzn;
	sig[2][2] = nzn;
	sig = mXm(cov, sig);
	cov = mXm(sig, trans(cov));
	return cov;
}


Vec xyz2rfl(Vec& xyz, Sensor tel) {
	Date date(xyz.jd);
	Vec rfl;
	for (int i = 0; i < 3; i++) {
		rfl.vec[0] += xyz.vec[i] * xyz.vec[i];
		rfl.vec[i + 3] = xyz.vec[i + 3];
	}
	rfl.vec[0] = sqrt(rfl.vec[0]);
	rfl.jd = xyz.jd;
	Doub telICRF[3];
	GetTELPos(date.day, date.time, &tel.pos[0], telICRF);
	ConvertXYZtoRADEC(&xyz.vec[0], telICRF, &rfl.vec[1], &rfl.vec[2]);
	return rfl;
}

Vec rfl2xyz(Vec& rfl, Sensor tel) {
	VecDoub icrfTel = tel.icrfPos(rfl.jd);
	Doub dt = 0.;
	VecDoub x0 = radec2xyz(rfl.vec[1], rfl.vec[2]);
	for (int i = 0; i < 3; i++) {
		x0[i] = x0[i] * rfl.vec[0] + icrfTel[i];
	}
	return Vec();

}

Vec prelimWithFixedD(Observs track, Doub d, Int startBlip = 0) {
	VecDoub icrfTel = track.sen.icrfPos(track.obs.at(startBlip).jd);
	Doub dt = 0.;
	VecDoub x0 = radec2xyz(track.obs.at(startBlip).vec[0], track.obs.at(startBlip).vec[1]);
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
	while ( it < track.obs.size() && abs(dt) < 31.) {
		dt = 60. * 60. * 24. * (track.obs.at(it).jd - track.obs.at(startBlip).jd);
		xs = track.sen.icrfPos(track.obs.at(it).jd);
		for (Int i = 0; i < 3; i++) {
			x[i] = -dt * (x0[i] - xs[i] + 0.5 * dt * dt * g[i]);
		}
		e = radec2xyz(track.obs.at(it).vec[0], track.obs.at(it).vec[1]);
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
	/*
	MatDoub jac = radec2xyzJac(track.obs.first());
	MatDoub blipCov(2, 2, 0.);
	blipCov[0][0] = track.obs.t.cov.vec[0];
	blipCov[1][1] = track.obs.t.cov.vec[1];
	MatDoub tmp = mult(jac, blipCov);
	MatDoub xyzCov = mult(tmp, trans(jac));
	*/
	Vec orb;
	for (Int i = 0; i < 3; i++) {
		orb.vec[i] = x0[i];
	}
	for (Int i = 3; i < 6; i++) {
 		orb.vec[i] = v[i - 3] * 1000;
	}
	orb.jd = track.obs.at(startBlip).jd;
	return orb;//preciseExt(orb, track.obs[0].jd);
}

Vec prelimWithFixedDLast(Observs track, Doub d) {
	VecDoub icrfTel = track.sen.icrfPos(track.obs.back().jd);
	Doub dt = 0.;
	VecDoub x0 = radec2xyz(track.obs.back().vec[0], track.obs.back().vec[1]);
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

	Int it = track.obs.size() - 1;
	while (it >= 0 && abs(dt) < 31.) {
		dt = 60. * 60. * 24. * (track.obs.at(it).jd - track.obs.back().jd);
		xs = track.sen.icrfPos(track.obs.at(it).jd);
		for (Int i = 0; i < 3; i++) {
			x[i] = -dt * (x0[i] - xs[i] + 0.5 * dt * dt * g[i]);
		}
		e = radec2xyz(track.obs.at(it).vec[0], track.obs.at(it).vec[1]);
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
		it--;
	}

	SVD svd(a);
	VecDoub v(3);
	svd.solve(y, v);
	/*
	MatDoub jac = radec2xyzJac(track.obs.first());
	MatDoub blipCov(2, 2, 0.);
	blipCov[0][0] = track.obs.t.cov.vec[0];
	blipCov[1][1] = track.obs.t.cov.vec[1];
	MatDoub tmp = mult(jac, blipCov);
	MatDoub xyzCov = mult(tmp, trans(jac));
	*/
	Vec orb;
	for (Int i = 0; i < 3; i++) {
		orb.vec[i] = x0[i];
	}
	for (Int i = 3; i < 6; i++) {
 		orb.vec[i] = v[i - 3] * 1000;
	}
	orb.jd = track.obs.back().jd;
	return orb;
}

Est estWithFixedD(Observs track, Doub d, Int startBlip = 0) {
	VecDoub icrfTel = track.sen.icrfPos(track.obs.at(startBlip).jd);
	Doub dt = 0.;
	VecDoub x0 = radec2xyz(track.obs.at(startBlip).vec[0], track.obs.at(startBlip).vec[1]);
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

	Int it = startBlip;
	while ( it < 2 && track.obs.size() && dt < 31.) {
		dt = 60. * 60. * 24. * (track.obs.at(it).jd - track.obs.at(startBlip).jd);
		xs = track.sen.icrfPos(track.obs.at(it).jd);
		for (Int i = 0; i < 3; i++) {
			x[i] = -dt * (x0[i] - xs[i] + 0.5 * dt * dt * g[i]);
		}
		e = radec2xyz(track.obs.at(it).vec[0], track.obs.at(it).vec[1]);
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
	/*
	MatDoub jac = radec2xyzJac(track.obs.first());
	MatDoub blipCov(2, 2, 0.);
	blipCov[0][0] = track.obs.t.cov.vec[0];
	blipCov[1][1] = track.obs.t.cov.vec[1];
	MatDoub tmp = mult(jac, blipCov);
	MatDoub xyzCov = mult(tmp, trans(jac));
	*/
	Est orb;
	for (Int i = 0; i < 3; i++) {
		orb.vec[i] = x0[i];
	}
	for (Int i = 3; i < 6; i++) {
 		orb.vec[i] = v[i - 3] * 1000;
	}
	orb.jd = track.obs.at(startBlip).jd;

	MatDoub rcov, vcov;
	//dcov
	Doub ra = track.obs.at(startBlip).vec[0];
	Doub dec = track.obs.at(startBlip).vec[1];
	Doub cosra = cos(ra), sinra = sin(ra), cosdec = cos(dec), sindec = sin(dec);
	x0 = radec2xyz(ra, dec);;
	MatDoub rot(3, 3);
	rot[0][0] = x0[0];
	rot[1][0] = x0[1];
	rot[2][0] = x0[2];
	rot[0][1] = -sinra*cosdec;
	rot[1][1] = cosra*cosdec;
	rot[2][1] = 0;
	rot[0][2] = -cosra*sindec;
	rot[1][2] = -sinra*sindec;
	rot[2][2] = cosdec;
	MatDoub sig(3, 3, 0.);
	Doub dRa = track.obs.at(startBlip).cov[0][0], dDec = track.obs.at(startBlip).cov[1][1];
	sig[0][0] = 1e-15 * 1. * 3.;//3 * rStep * rStep;
	sig[1][1] = 3 * sqrt(dRa) * d;
	sig[2][2] = 3 * sqrt(dDec) * d;
	sig = mXm(rot, sig);
	rot = mXm(sig, trans(rot));
	rcov = rot;
	
	//vcov
	rot[0][0] = x0[0];
	rot[1][0] = x0[1];
	rot[2][0] = x0[2];
	v = normilize(v);
	VecDoub n = crossProduct(x0, v);
	//n = normilize(n);
	rot[0][1] = n[0];
	rot[1][1] = n[1];
	rot[2][1] = n[2];
	v = crossProduct(x0, n);
	//v = normilize(v);
	rot[0][2] = v[0];
	rot[1][2] = v[1];
	rot[2][2] = v[2];
	sig = MatDoub(3, 3, 0.);
	sig[0][0] = 12;
	sig[1][1] = 1e-1;
	sig[2][2] = 12;
	sig = mXm(rot, sig);
	rot = mXm(sig, trans(rot));
	vcov = rot;

	return orb;//preciseExt(orb, track.obs[0].jd);
}

void bestBlip(Observs track, Doub d) {
	Int startBlip = 0;
	VecDoub icrfTel = track.sen.icrfPos(track.obs.at(startBlip).jd);
	Doub dt = 0.;

	Doub dRa = track.obs.at(startBlip).cov[0][0], dDec = track.obs.at(startBlip).cov[1][1];
	Plot3D p3d;
	p3d.point(track.sen.toObs(vec1).vec[0], track.sen.toObs(vec1).vec[1], 0, 3);
	p3d.point(track.obs.at(startBlip).vec[0] + 0 / 100. * 3 * sqrt(dRa), track.obs.at(startBlip).vec[1] + 0 / 100. * 3 * sqrt(dDec), 1, 3);
	for (int i = -99; i < 100; i+=19) {
		for (int j = -99; j < 100; j+=19) {
			VecDoub x0 = radec2xyz(track.obs.at(startBlip).vec[0] + i / 100. * 3 * sqrt(dRa), track.obs.at(startBlip).vec[1] + j / 100. * 3 * sqrt(dDec));
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

			Int it = startBlip;
			while ( it < track.obs.size() && dt < 16.) {
				dt = 60. * 60. * 24. * (track.obs.at(it).jd - track.obs.at(startBlip).jd);
				xs = track.sen.icrfPos(track.obs.at(it).jd);
				for (Int i = 0; i < 3; i++) {
					x[i] = -dt * (x0[i] - xs[i] + 0.5 * dt * dt * g[i]);
				}
				e = radec2xyz(track.obs.at(it).vec[0], track.obs.at(it).vec[1]);
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
			/*
			MatDoub jac = radec2xyzJac(track.obs.first());
			MatDoub blipCov(2, 2, 0.);
			blipCov[0][0] = track.obs.t.cov.vec[0];
			blipCov[1][1] = track.obs.t.cov.vec[1];
			MatDoub tmp = mult(jac, blipCov);
			MatDoub xyzCov = mult(tmp, trans(jac));
			*/
			Vec orb;
			for (Int i = 0; i < 3; i++) {
				orb.vec[i] = x0[i];
			}
			for (Int i = 3; i < 6; i++) {
 				orb.vec[i] = v[i - 3] * 1000;

			}
			orb.jd = track.obs.front().jd;

			Doub h = orbitEnergy(orb);
			if (0) {
				VecDoub v(3), r(3);
				for (Int i = 0; i < 3; i++) {
					v[i] = orb.vec[i+3];
					r[i] = orb.vec[i];
				}
				v = normilize(v);
				Doub k = 398.6004415;	//gravity parameter km3 / s2
				Doub vMod = sqrt(2. * k / scalar(r));

				for (Int i = 0; i < 3; i++) {
					orb.vec[i+3] = v[i] * vMod;
				}
			}


			Doub res;
			//print(orbitEnergy(orb));
			if (isOrbit(orb) > 0.) {
				orb.print();
				//print("kodekdok");
				res = residual(orb, track);
				p3d.point(track.obs.at(startBlip).vec[0] + i / 100. * 3 * sqrt(dRa), track.obs.at(startBlip).vec[1] + j / 100. * 3 * sqrt(dDec), res);
			}
			
		}
			p3d.show();

	}

}
