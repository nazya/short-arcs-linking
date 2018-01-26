void print(Vec orb, FILE* p = stdout) {
	fprintf(p, "vec: ");
	for (Int i = 0; i < orb.vec.size(); i++) {
		fprintf(p, "%+6.4lf ", orb.vec[i]);
	}
	Date date(orb.jd);
	fprintf(p, "# %.0lf %.0lf", date.day, date.time);
	fprintf(p, "\n");
	
}

void print(vector<Doub> rav) {
	if (rav.size() == 0)
		printf("null");
	for (int i = 0; i < rav.size(); i++) {
		printf("%lf ", rav[i]);
	}
	printf("\n");
}
Sensor initTelescope (Char* filePath) { //init telescope
	VecDoub pos(3);
	FILE* in = fopen(filePath, "r");
	for (Int i = 0; i < 3; i++) {
		fscanf(in, "%lf ", &pos[i]);
	}
	Sensor tel(pos);
	fclose(in);
	return tel;
}

Vec initVec (Char* filePath) {
	Vec vec;
	FILE* in = fopen(filePath, "r");
    fscanf(in, "%lf", &vec.jd);

	for (Int i = 0; i < 6; i++) {
		fscanf(in, "%lf", &(vec.vec[i]));
	}
	fclose(in);
	return vec;
}

Est initEst (Char* filePath) {
	Est est;
	FILE* in = fopen(filePath, "r");
    fscanf(in, "%lf", &est.jd);

	for (Int i = 0; i < 6; i++) {
		fscanf(in, "%lf", &(est.vec[i]));
	}
	
	for (Int i = 0; i < 6; i++) {
        for (Int j = 0; j < 6; j++) {
            fscanf(in, "%lf", &est.cov[i][j]);
		}
	}
	fclose(in);
	return est;
}

void println() {
	printf("\n");
}

void print(Doub m) {
	printf("%+.25lf\n", m);
}


void print(MatDoub m) {
	for (Int i = 0; i < m.nrows(); i++) {
		for (Int j = 0; j < m.ncols(); j++) {
			printf("%+6.6lf ", m[i][j]);
		}
		printf("\n");
	}
}

void print(VecDoub& v) {
	for (Int i = 0; i < v.size(); i++) {
		printf("%+6.16lf ", v[i]);
	}
	printf("\n");
}

void print(char* s) {
	printf("%s\n", s);
}

void saveOrbvec (Char* filePath, Vec& vec) {
	FILE* out = fopen(filePath, "w");
    fprintf(out, "%.15lf\n", vec.jd);
	for (Int i = 0; i < 6; i++) {
		fprintf(out, "%.15lf ", vec.vec[i]);
	}
	fclose(out);
}

void saveEstvec (Char* filePath, Est& est) {
	FILE* out = fopen(filePath, "w");
    fprintf(out, "%.15lf\n", est.jd);
	for (Int i = 0; i < 6; i++) {
		fprintf(out, "%.15lf ", est.vec[i]);
	}
	fprintf(out, "\n");
	for (Int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			fprintf(out, "%.15lf ", est.cov[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
}