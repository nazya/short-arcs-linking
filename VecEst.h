struct Vec {
	VecDoub vec;
	Doub jd;

	Vec(Int dim = 6) : vec(dim, 0.) {
	}

	void print() {
		print(stdout);
	}

	void print(FILE* p) {
		fprintf(p, "vec: ");
		for (Int i = 0; i < vec.size(); i++) {
			fprintf(p, "%+6.4lf ", vec[i]);
		}
		Date date(jd);
		fprintf(p, "# %.0lf %.0lf", date.day, date.time);
		fprintf(p, "\n");
	}
	
};

struct Est : Vec {
	MatDoub cov;

	Est() : cov(vec.size(), vec.size(), 0.) {
	}

	Est(Vec v) : cov(v.vec.size(), v.vec.size(), 0.) {
		vec = v.vec;
		jd = v.jd;
	}


	void print() {
		print(stdout);
	}

	void print(FILE* p) {
		fprintf(p, "date:\n%lf\nvec:\n", jd);
		for (Int i = 0; i < vec.size(); i++) {
			fprintf(p, "%+6.5lf ", vec[i]);
		}
		fprintf(p, "\ncov:\n");
		for (Int i = 0; i < vec.size(); i++) {
			for (int j = 0; j < vec.size(); j++) {
				fprintf(p, "%+6.4lf ", cov[i][j]);
			}
			fprintf(p, "\n");
		}
	}
};

typedef Vec Vec_I;
typedef Vec Vec_O;
typedef Vec Vec_IO;

typedef Est Est_I;
typedef Est Est_O;
typedef Est Est_IO;

Vec xyz2rflVec (Vec xyz) {
	Vec out;
	Doub r = 0;
	for (int i = 0; i < 3; i++) {
		r += xyz.vec[i] * xyz.vec[i];
		out.vec[i+3] = xyz.vec[i+3];
	}
	r = sqrt(r);
	out.vec[0] = r;
	Doub ra = atan2(xyz.vec[1], xyz.vec[0]);
	if( ra < 0 )
		ra = 2.0*3.1415926535 + ra;
	out.vec[1] = atan2(xyz.vec[1], xyz.vec[0]);
	out.vec[2] = atan2(xyz.vec[2], sqrt(xyz.vec[0]*xyz.vec[0] + xyz.vec[1]*xyz.vec[1]));
	out.jd = xyz.jd;
	return out;
}

Vec rfl2xyzVec(Vec rfl) {
	Vec out;
	Doub r = rfl.vec[0];
	VecDoub xyz(3);
	Doub ra = rfl.vec[1];
	Doub dec = rfl.vec[2];
	xyz[0] = r * cos(ra) * cos(dec);
	xyz[1] = r * sin(ra) * cos(dec);
	xyz[2] = r * sin(dec);


	for (int i = 0; i < 3; i++) {
		out.vec[i] = xyz[i];
		out.vec[i+3] = rfl.vec[i+3];
	}
	out.jd = rfl.jd;
	return out;
}