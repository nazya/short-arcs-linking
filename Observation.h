struct Observation {
	VecDoub obs;
	MatDoub cov;
	Doub jd;

	Observation(Int dim = 2) : obs(dim, 0.), cov(dim, dim, 0.) {
		cov[0][0] = 2.35e-12;
		cov[1][1] = 2.35e-12;
	}

	void print(FILE* p = stdout) {
		fprintf(p, "blip: ");
		fprintf(p, "ra: %+6.5lf dec: %+6.5lf jdate: %+6.5lf", obs[0], obs[1], jd);
		fprintf(p, "\n");
	}

};
