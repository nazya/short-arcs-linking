Doub dOrb(Vec ov, Sensor sen) {

	VecDoub icrf = sen.icrfPos(ov.jd);
	Doub d;
	for (int k = 0; k < 3; k++) {
		icrf[k] = ov.vec[k] - icrf[k];
	}
	return sqrt(innerProduct(icrf));
}

Doub rOrb(Vec& ov) {
	Doub r = 0.;		//distance
	for (Int i = 0; i < 3; i++) {
		r += ov.vec[i] * ov.vec[i];
	}
	return sqrt(r);
}


vector<Vec> fakeJac(Vec vec1, Observs track1, Observs track2, Est* lp) {
	Est paramsEst;// = estimate(track1);
	paramsEst.print();

	MatDoub varcov(6, 6, 0.);

	for (int i = 0; i < 3; i++) {
		varcov[i][i] = 1e-11;
		varcov[i+3][i+3] = 1e-11;
	}
	VecDoub paramsDev(6);
	Multinormaldev dev(124, vec1.vec, varcov);
	vector<Vec> prevars;
	vector<Doub> reses;
	paramsEst.jd = vec1.jd;
	while(prevars.size() < 7) {
		paramsEst.vec = dev.dev();
		Doub res = residual(paramsEst, track2);
		printf("%lf\n", res);
		if(res < 3.) {
			prevars.push_back(paramsEst);
			reses.push_back(res);
	
		}
	}

	Doub bres = 9.;
	Int b;
	for (int i = 0; i < 6; i++) {
		if(reses[i] < bres) {
			b = i;
			bres = reses[i];
		}
	}

	paramsEst.vec = prevars[b].vec;
	vector<Vec> vars;
	for (int i = 0; i < 7; i++) {
		if (i != b) {
			vars.push_back(prevars[i]);
		}
	}
	*lp = paramsEst;
	return vars;

}

Doub rCov(Est est) {
	
	Doub r = rOrb(est);
	MatDoub jac(1, 3);
	MatDoub xyzcov(3, 3);
	for (int i = 0; i < 3; i++) {
		jac[0][i] = est.vec[i] / r;
		for (int j = 0; j < 3; j++) {
			xyzcov[i][j] = est.cov[i][j];
		}
	}

	print(mXm(mXm(jac, xyzcov), trans(jac)));

	return sqrt(mXm(mXm(jac, xyzcov), trans(jac))[0][0]);
}


void dv(Vec vec1, Est paramsEst, Observs track1) {
	Plot3D dplt, vplt;
	VecDoub vel(3), dis(3), dtrue(3), vtrue(3);
	MatDoub vcov(3, 3), dcov(3, 3);
	for (int i = 0; i < 3; i++) {
		dtrue[i] = vec1.vec[i];
		vtrue[i] = vec1.vec[i+3];
		dis[i] = paramsEst.vec[i];
		vel[i] = paramsEst.vec[i+3];
		for (int j = 0; j < 3; j++) {
			dcov[i][j] = paramsEst.cov[i][j];
			vcov[i][j] = paramsEst.cov[i+3][j+3];
		}
	}
	dplt.point(dtrue);
	dplt.ellipsoid(dis, dcov);
	vplt.point(vtrue, 3);
	vplt.ellipsoid(vel, vcov);
	dplt.show();
	vplt.show();

	vector<VecDoub> points;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	
	

	Doub d = .0;
	Observs track = track1;
	Int startBlip = 0;
	int n = 0;
	Int nR = 0, nRa = n, nDec = n;
	Doub ra0 = track.obs.at(startBlip).vec[0], dec0 = track.obs.at(startBlip).vec[1];
	Doub r = rOrb(vec1), ra, dec;
	VecDoub icrf = track.sen.icrfPos(track.obs.at(startBlip).jd);
	for (int k = 0; k < 3; k++) {
		icrf[k] = vec1.vec[k] - icrf[k];
	}
	for (int k = 0; k < 3; k++) {
		d += icrf[k] * icrf[k];
	}
	d = sqrt(d);
	r = d;
	
	Doub dRa = track.obs.at(startBlip).cov[0][0], dDec = track.obs.at(startBlip).cov[1][1];
	dRa = 3 * sqrt(dRa) / (2 * nRa + 1), dDec = 3 * sqrt(dDec) / (2 * nDec + 1);
	for (int iR = -nR; iR < nR + 1; iR++) {
		for (int iRa = -nRa; iRa < nRa + 1; iRa++) {
			for (int iDec = -nDec; iDec < nDec + 1; iDec++) {
				//form fl grid
				//r = d + iR * 1e-1;
				ra = ra0 + iRa * dRa;
				dec = dec0 + iDec * dDec;
				//form xyz
				VecDoub icrfTel = icrf = track.sen.icrfPos(track.obs.at(startBlip).jd);
				Doub dt = 0.;
				VecDoub x0 = radec2xyz(ra, dec);
				for (int i = 0; i < 3; i++) {
					x0[i] = x0[i] * r + icrfTel[i];
				}
				VecDoub x(3), e(3), xs(3), y(3, 0.), pXx, g = gravityVector(x0), v(3);
				MatDoub p, a(3, 3, 0.);
				Int it = startBlip;
				while ( it < track.obs.size() && dt < 30.) {
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
				svd.solve(y, v);
				Est orb;
				for (Int i = 0; i < 3; i++) {
					orb.vec[i] = x0[i];
				}
				for (Int i = 3; i < 6; i++) {
 					orb.vec[i] = v[i - 3] * 1000;
				}
				orb.jd = track.obs.at(startBlip).jd;
				vplt.point(vel);
				//dcov
				Doub cosra = cos(ra), sinra = sin(ra), cosdec = cos(dec), sindec = sin(dec);
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
				sig[0][0] = 1e-15;
				sig[1][1] = dRa * dRa * 9. * 9. ;
				sig[2][2] = dDec * dDec * 9. * 9.;
				sig = mXm(cov, sig);
				cov = mXm(sig, trans(cov));
				for (Int i = 0; i < 3; i++) {
					for (Int j = 0; j < 3; j++) {
						orb.cov[i][j] = cov[i][j];
					}
				}
				for (Int i = 3; i < 6; i++) {
					for (Int j = 3; j < 6; j++) {
						orb.cov[i][i] = paramsEst.cov[i][j];
					}
				}
				for (int i = 0; i < 3; i++) {
					dis[i] = orb.vec[i];
					vel[i] = orb.vec[i+3];
					for (int j = 0; j < 3; j++) {
						dcov[i][j] = orb.cov[i][j];
						vcov[i][j] = orb.cov[i+3][j+3];
					}
				}
				vplt.ellipsoid(vel, vcov);
				vplt.show();
			}
		}
	}

}


int diff(Vec v1, Vec v2) {
	Vec out;
	for (Int i = 0; i < 6; i++) {
		out.vec[i] = abs(v1.vec[i] - v2.vec[i]);
	}

	Doub v = 0.;		//velocity
	Doub r = 0.;		//distance
	for (Int i = 0; i < 3; i++) {
		r += out.vec[i];
		v += out.vec[i + 3];
	}
	v = sqrt(v);
	r = sqrt(r);
	if (v < 0.01 && r < 0.000001)
		return 1;
	return 0;
}

vector<VecDoub> fixeDVelRavine4(Vec orb, VecDoub evrav) {
	evrav = normilize(evrav);
	Doub v02 = 0.;		//velocity
	Doub r = 0.;		//distance
	Doub k = 398600.4415;	//gravity parameter km3 / s2
	//Doub b = 0.;
	VecDoub pos(3), v0(3);
	for (Int i = 0; i < 3; i++) {
		pos[i] = orb.vec[i] * 1000.;
		v0[i] = orb.vec[i+3];
		r += orb.vec[i] * orb.vec[i];
		v02 += orb.vec[i + 3] * orb.vec[i + 3];
	//	b += evrav[i] * v0[i];
	}
	//b *= 2;
	r = sqrt(r) * 1000.;
	VecDoub nxr, v0xr;
	Doub nv0, nxr2, v0xr2, v0xrnxr, k2, k3;
	k2 = k * k;
	k3 = k2 * k;
	nxr = crossProduct(evrav, pos);
	nxr2 = innerProduct(nxr);
	v0xr = crossProduct(v0, pos);
	v0xr2 = innerProduct(v0xr);
	nv0 = innerProduct(evrav, v0);
	v0xrnxr = innerProduct(v0xr, nxr);
	Doub down, up;
	//energy gate: max energy => max v2
	/*
	Doub c = v02 - 2. * k / r;
	up = (-b - sqrt(b * b - 4 * c)) / 2.;
	down = (-b + sqrt(b * b - 4 * c)) / 2.;
	*/
	Doub abc[3];
	Doub per = 6370. + 150., ap = 25000.;
	Doub per2 = per * per, ap2 = ap * ap;
	vector<Doub> gate;
	abc[0] = nxr2 + ap2;
	abc[1] = 2. * (v0xrnxr + ap2 * nv0 );
	abc[2] = v0xr2 + ap2 * v02 - 2. * ap2 * k / r;
	gate.push_back((-abc[1] + sqrt(abc[1] * abc[1] - 4 * abc[0] * abc[2])) / 2. / abc[0]);
	print(gate.back());
	gate.push_back((-abc[1] - sqrt(abc[1] * abc[1] - 4 * abc[0] * abc[2])) / 2. / abc[0]);
	print(gate.back());

	abc[0] = nxr2 - per2;
	abc[1] = 2. * (v0xrnxr - nv0 * per2);
	abc[2] = v0xr2 - v02 * per2 + 2. * per2 * k / r - 2. * per * k;
	gate.push_back((-abc[1] + sqrt(abc[1] * abc[1] - 4. * abc[0] * abc[2])) / 2. / abc[0]);
	print(gate.back());
	gate.push_back((-abc[1] - sqrt(abc[1] * abc[1] - 4. * abc[0] * abc[2])) / 2. / abc[0]);
	print(gate.back());
	for (int i = 0; i < 3; i++) {
		for (int j = i + 1; j < gate.size(); j++) {
			if (gate[j] > gate[i]) {
				gate[i] += gate[j];
				gate[j] = gate[i] - gate[j];
				gate[i] -= gate[j];
			}
		}
	}
	vector<VecDoub> ret;


	//for (int j = 0; j < gate.size(); j++) {
	//	for (Int i = 0; i < 3; i++) {
	//		vel[i] = orb.vec[i + 3] + gate[j] * evrav[i];
	//	}
	//	ret.push_back(vel);
	//}
	//return ret;

	up = gate[1];
	down = gate[2];
	for (Int i = 0; i < 3; i++) {
		vel[i] = orb.vec[i + 3] + down * evrav[i];
	}
	ret.push_back(vel);
	for (Int i = 0; i < 3; i++) {
		vel[i] = orb.vec[i + 3] + up * evrav[i];
	}
	ret.push_back(vel);




	return ret;
}



vector<VecDoub> fixeDVelRavine2(Vec orb, VecDoub evrav) {
	evrav = normilize(evrav);
	Doub v02 = 0.;		//velocity
	Doub r = 0.;		//distance
	Doub k = 398600.4415;	//gravity parameter km3 / s2
	//Doub b = 0.;
	VecDoub pos(3), v0(3);
	for (Int i = 0; i < 3; i++) {
		pos[i] = orb.vec[i] * 1000.;
		v0[i] = orb.vec[i+3];
		r += orb.vec[i] * orb.vec[i];
		v02 += orb.vec[i + 3] * orb.vec[i + 3];
	//	b += evrav[i] * v0[i];
	}
	r = sqrt(r) * 1000.;

	Doub abc[3], D;
	vector<Doub> gate;
	
	Doub n2, v2, vn;
	n2 = innerProduct(evrav);
	vn = innerProduct(evrav, v0);
	v2 = innerProduct(v0);

	//print(k / (6370. + 100.));
	//print(k / (50000.));

	abc[0] = n2;
	abc[1] = 2. * vn;
	abc[2] = v2 - 2 * k / r + k / (50000.);
	D = abc[1] * abc[1] - 4 * abc[0] * abc[2];
	if (D > 0) {
		gate.push_back((-abc[1] + sqrt(abc[1] * abc[1] - 4 * abc[0] * abc[2])) / 2. / abc[0]);
		//print(gate.back());
		gate.push_back((-abc[1] - sqrt(abc[1] * abc[1] - 4 * abc[0] * abc[2])) / 2. / abc[0]);
		//print(gate.back());
	}
	vector<VecDoub> ret;
	if(gate.size() < 2)
		return ret;
	for (Int i = 0; i < 3; i++) {
		vel[i] = orb.vec[i + 3] + gate[0] * evrav[i];
	}
	ret.push_back(vel);
	for (Int i = 0; i < 3; i++) {
		vel[i] = orb.vec[i + 3] + gate[1] * evrav[i];
	}
	ret.push_back(vel);
	return ret;

}
vector<Doub> crossGate(vector<vector<Doub>> in) {
	vector<Doub> c;
	for (int i = in.size() - 1; i > 0 ; i--) {
		Int isCChanged = 0;
		vector<Doub> a = in[i], b = in[i-1];
		c.clear();
		if (a.size() == 0 || b.size() == 0) {
			//print("cross");
			//print(c);
			//print("");
			return c;
		}
		for (int k = 0, j = 0; k < b.size() &&  j < a.size();) {
			//print(j);
			//print(k);
			isCChanged = 0;
			if (a[j] >= b[k]) {
				if (b[k+1] >= a[j]) {
					c.push_back(a[j]);
					isCChanged = 1;
					if (b[k+1] <= a[j+1]) {
						c.push_back(b[k+1]);
						isCChanged = 1;
					} else {
						c.push_back(a[j+1]);
						isCChanged = 1;
					}
				}
			} else {
				if (a[j+1] >= b[k]) {
					c.push_back(b[k]);
					isCChanged = 1;
					if (a[j+1] <= b[k+1]) {
						c.push_back(a[j+1]);
						isCChanged = 1;
					} else {
						c.push_back(b[k+1]);
						isCChanged = 1;
					}
				}
			}

			if(c.size() == 0 || !isCChanged) {
				if (a[j+1] <= b[k])
					j+=2;
				else if ( b[k+1] <= a[j])
					k+=2;
			} else {
				if (a[j+1] == c.back())
					j+=2;
				else if (b[k+1] == c.back())
					k+=2;
			}
			isCChanged = 0;

		}
		in.pop_back();
		in.pop_back();
		in.push_back(c);
		
	}
	//print("cross");
	//print(in[0]);
	//print("");
	return in[0];
}

vector<Doub> quadraticInequality(Doub* abc, Doub sign) {
	vector<Doub> tmp;
	Doub D = abc[1] * abc[1] - 4. * abc[0] * abc[2];
	if (abc[0] * sign < 0) {
		if (D > 0) {
			tmp.push_back((-abc[1] + sign * sqrt(abc[1] * abc[1] - 4. * abc[0] * abc[2])) / 2. / abc[0]);
			//print(tmp.back());
			tmp.push_back((-abc[1] - sign * sqrt(abc[1] * abc[1] - 4. * abc[0] * abc[2])) / 2. / abc[0]);
			//print(tmp.back());
		} else {
			//print("0");
		}
	} else if (abc[0] * sign > 0) {
		if (D > 0) {
			tmp.push_back(-DBL_MAX);
			//print(tmp.back());
			tmp.push_back((-abc[1] - sign * sqrt(abc[1] * abc[1] - 4. * abc[0] * abc[2])) / 2. / abc[0]);
			//print(tmp.back());
			tmp.push_back((-abc[1] + sign * sqrt(abc[1] * abc[1] - 4. * abc[0] * abc[2])) / 2. / abc[0]);
			//print(tmp.back());
			tmp.push_back(DBL_MAX);
			//print(tmp.back());
		} else {
			tmp.push_back(-DBL_MAX);
			//print(tmp.back());
			tmp.push_back(DBL_MAX);
			//print(tmp.back());
		}
	} else {
		print("quadraticInequality a=0 exception");
		exit(0);
	}
	return tmp;
}

vector<VecDoub> fixeDVelRavine3saved(Vec orb, VecDoub evrav) {
	evrav = normilize(evrav);
	Doub v02 = 0.;		//velocity
	Doub r = 0.;		//distance
	Doub k = 398600.4415;	//gravity parameter km3 / s2
	//Doub b = 0.;
	VecDoub pos(3), v0(3);
	for (Int i = 0; i < 3; i++) {
		pos[i] = orb.vec[i] * 1000.;
		v0[i] = orb.vec[i+3];
		r += orb.vec[i] * orb.vec[i];
		v02 += orb.vec[i + 3] * orb.vec[i + 3];
	//	b += evrav[i] * v0[i];
	}
	//b *= 2;
	r = sqrt(r) * 1000.;
	VecDoub nxr, v0xr;
	Doub nv0, nxr2, v0xr2, v0xrnxr;
	
	nxr = crossProduct(pos, evrav);
	nxr2 = innerProduct(nxr);
	v0xr = crossProduct(pos, v0);
	v0xr2 = innerProduct(v0xr);
	nv0 = innerProduct(evrav, v0);
	v0xrnxr = innerProduct(v0xr, nxr);
	
	Doub abc[3];
	Doub per = 6370. + 100.;
	Doub ap = 50000;
	Doub per2 = per * per, ap2 = ap * ap;
	vector<vector<Doub>> gate;
	abc[0] = nxr2 - ap2;
	abc[1] = 2. * (v0xrnxr - ap2 * nv0 );
	abc[2] = v0xr2 - ap2 * v02 + 2. * ap2 * k / r - 2. * k * ap;
	gate.push_back(quadraticInequality(abc, +1));
	abc[0] = nxr2;
	abc[1] = 2. * (v0xrnxr);
	abc[2] = v0xr2 - ap * k;
	//gate.push_back(quadraticInequality(abc, -1));
	
	abc[0] = 1.;
	abc[1] = 2. * (nv0);
	abc[2] = v02 - 2 * k / r + k / per;
	gate.push_back(quadraticInequality(abc, +1));
	abc[0] = 1;
	abc[1] = 2. * (nv0);
	abc[2] = v02 - 2 * k / r + k / 50000.;
	gate.push_back(quadraticInequality(abc, -1));


	abc[0] = nxr2 - per2;
	abc[1] = 2. * (v0xrnxr - nv0 * per2);
	abc[2] = v0xr2 - v02 * per2 + 2. * per2 * k / r - 2. * per * k;
	gate.push_back(quadraticInequality(abc, 1.));
	abc[0] = nxr2;
	abc[1] = 2. * (v0xrnxr);
	abc[2] = v0xr2 - k * per;
	//gate.push_back(quadraticInequality(abc, 1.));
	
	vector<Doub> tmp = crossGate(gate);
	VecDoub v(3);
	vector<VecDoub> ret;
	for (int j = 0; j < tmp.size(); j++) {
		for (int i = 0; i < 3; i++) {
			v[i] = v0[i] + tmp[j] * evrav[i];
		}
		ret.push_back(v);
	}
	if (ret.size() > 2) {
		print("exception fixdravine3");
		exit(0);
	}
	return ret;
}

vector<VecDoub> fixeDVelRavine3(Vec orb, VecDoub evrav) {
	evrav = normilize(evrav);
	Doub v02 = 0.;		//velocity
	Doub r = 0.;		//distance
	Doub k = 398600.4415;	//gravity parameter km3 / s2
	//Doub b = 0.;
	VecDoub pos(3), v0(3);
	for (Int i = 0; i < 3; i++) {
		pos[i] = orb.vec[i] * 1000.;
		v0[i] = orb.vec[i+3];
		r += orb.vec[i] * orb.vec[i];
		v02 += orb.vec[i + 3] * orb.vec[i + 3];
	//	b += evrav[i] * v0[i];
	}
	//b *= 2;
	r = sqrt(r) * 1000.;
	VecDoub nxr, v0xr;
	Doub nv0, nxr2, v0xr2, v0xrnxr;
	
	nxr = crossProduct(pos, evrav);
	nxr2 = innerProduct(nxr);
	v0xr = crossProduct(pos, v0);
	v0xr2 = innerProduct(v0xr);
	nv0 = innerProduct(evrav, v0);
	v0xrnxr = innerProduct(v0xr, nxr);
	
	Doub abc[3];
	Doub per = 6370. + 100.;
	Doub ap = 100000;
	Doub per2 = per * per, ap2 = ap * ap;
	vector<vector<Doub>> gate;
	abc[0] = nxr2 - ap2;
	abc[1] = 2. * (v0xrnxr - ap2 * nv0 );
	abc[2] = v0xr2 - ap2 * v02 + 2. * ap2 * k / r - 2. * k * ap;
	gate.push_back(quadraticInequality(abc, +1));
	
	abc[0] = 1.;
	abc[1] = 2. * (nv0);
	abc[2] = v02 - 2 * k / r + k / per;
	gate.push_back(quadraticInequality(abc, +1));
	abc[0] = 1;
	abc[1] = 2. * (nv0);
	abc[2] = v02 - 2 * k / r + k / ap;
	gate.push_back(quadraticInequality(abc, -1));


	abc[0] = nxr2 - per2;
	abc[1] = 2. * (v0xrnxr - nv0 * per2);
	abc[2] = v0xr2 - v02 * per2 + 2. * per2 * k / r - 2. * per * k;
	gate.push_back(quadraticInequality(abc, 1.));
	
	vector<Doub> tmp = crossGate(gate);
	VecDoub v(3);
	vector<VecDoub> ret;
	for (int j = 0; j < tmp.size(); j++) {
		for (int i = 0; i < 3; i++) {
			v[i] = v0[i] + tmp[j] * evrav[i];
		}
		ret.push_back(v);
	}
	if (ret.size() > 2) {
		print("exception fixdravine3");
		exit(0);
	}
	return ret;
}

vector<VecDoub> fixeDVelRavine(Vec orb, VecDoub evrav) {
	evrav = normilize(evrav);
	Doub v02 = 0.;		//velocity
	Doub r = 0.;		//distance
	Doub k = 398600.4415;	//gravity parameter km3 / s2
	//Doub b = 0.;
	VecDoub pos(3), v0(3);
	for (Int i = 0; i < 3; i++) {
		pos[i] = orb.vec[i] * 1000.;
		v0[i] = orb.vec[i+3];
		r += orb.vec[i] * orb.vec[i];
		v02 += orb.vec[i + 3] * orb.vec[i + 3];
	//	b += evrav[i] * v0[i];
	}
	//b *= 2;
	r = sqrt(r) * 1000.;
	VecDoub nxr, v0xr;
	Doub nv0, nxr2, v0xr2, v0xrnxr;
	
	nxr = crossProduct(pos, evrav);
	nxr2 = innerProduct(nxr);
	v0xr = crossProduct(pos, v0);
	//v0xr2 = innerProduct(v0xr);
	nv0 = innerProduct(evrav, v0);
	//v0xrnxr = innerProduct(v0xr, nxr);
	
	Doub abc[3];
	Doub per = 6370. + 100.;
	Doub ap = 50000;
	Doub per2 = per * per, ap2 = ap * ap;
	
	abc[0] = 1;
	abc[1] = 2. * (nv0);
	abc[2] = v02 - 2 * k / r;
	vector<Doub> tmp = quadraticInequality(abc, -1);

	
	VecDoub v(3);
	vector<VecDoub> ret;
	for (int j = 0; j < tmp.size(); j++) {
		for (int i = 0; i < 3; i++) {
			v[i] = v0[i] + tmp[j] * evrav[i];
		}
		ret.push_back(v);
	}
	return ret;
}
