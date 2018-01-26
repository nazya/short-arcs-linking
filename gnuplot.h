static Int gnuplotfile = 0;

struct Graph {
	string filename;
	Int style;
	Char type;
	Int dim;

	Graph(string filename, Char type, Int style, Int dim) : filename(filename), style(style), type(type), dim(dim) {
	}

	Graph add(vector<VecDoub> points) {
		FILE* temp = fopen(filename.c_str(), "w");
		for (int i = 0; i < points.size(); i++) {
			for (int j = 0; j < dim; j++) {
				fprintf(temp, "%.24lf ", points[i][j]); //Write the data to a temporary file
			}
			fprintf(temp, "\n");
		}
		fclose(temp);
		return *this;
	}

	string tocmd() {
		if (type == 'l')
			return " '" + filename + "'"  + " with lines ls " + to_string(style);
			//return " '" + filename + "'";//  + " with lines ls " + to_string(style);
		else if (type == 'p')
			return " '" + filename + "'" + " with points ls " + to_string(style);
			//return " '" + filename + "'";// + " with points ls " + to_string(style);
        else
            return " '" + filename + "'";// + "with points ls " + to_string(4);
	}
};

struct Plot {
	string path;
	FILE* gnuplotPipe;
	vector<Graph> graphs;

	Plot() {
		gnuplotPipe = popen ("gnuplot", "w");//_popen ("gnuplot -persistent", "w");
		path = "";
		gnuplotcmd("unset key");
		//gnuplotcmd("set view equal xyz");
		gnuplotcmd("set style line 1 lc rgb 'red' ps 0 lt 1 lw 0 pt 1");
		gnuplotcmd("set style line 2 lc rgb 'red' lt 1 lw 1 pt 1");
		gnuplotcmd("set style line 3 lc rgb 'blue' lt 1 lw 1 pt 2");
		gnuplotcmd("set style line 4 lc rgb 'black' lt 1 lw 1 pt 6 ps 3");
		gnuplotcmd("set style line 5 lc rgb 'black' lt 1 lw 1 ps 3");
		gnuplotcmd("set style line 6 lc rgb 'red' lt 1 lw 1 pt 0 ps 1");
		//gnuplotcmd("set dgrid3d");
		//gnuplotcmd("set pm3d corners2color c2");
	}

	void gnuplotcmd(const char* cmd) {
		fprintf(gnuplotPipe, "%s\n", cmd); //Send commands to gnuplot one by one.
		fflush(gnuplotPipe);
	}

	void pop_back() {
		remove((graphs.back()).filename.c_str());
		graphs.pop_back();
	}

	void del() {
		for (Int i = 0; i < graphs.size(); i++) {
			remove(graphs[i].filename.c_str());
		}
		for (Int i = 0; i < graphs.size(); i++) {
			graphs.pop_back();
		}
		graphs.clear();
	}

	void clear() {
		gnuplotcmd("clear\n");
		del();
	}

	~Plot(){
		clear();
		gnuplotcmd("q");
		fclose(gnuplotPipe);
	}
};

struct Plot3D : Plot {
	void point(Doub x, Doub y, Doub z, Int style = 2) {
		VecDoub v (3);
		v[0] = x;
		v[1] = y;
		v[2] = z;
		vector<VecDoub> points;
		points.push_back(v);
		graphs.push_back(Graph(path + to_string(++gnuplotfile),  'p', style, 3).add(points) );
	}

	void point(VecDoub v,  Int style = 2) {
		vector<VecDoub> points;
		points.push_back(v);
		graphs.push_back(Graph(path + to_string(++gnuplotfile),  'p', style, 3).add(points) );
	}

	void points(vector<VecDoub> points,  Int style = 1) {
		graphs.push_back(Graph(path + to_string(++gnuplotfile),  'p', style, 3).add(points) );
	}

	void line(vector<VecDoub> points,  Int style = 3) {
		graphs.push_back(Graph(path + to_string(++gnuplotfile),  'l', style, 3).add(points) );
	}

	void pellipsoid(VecDoub centre, MatDoub cov) {
		Multinormaldev testtdev(123, centre, cov);
		VecDoub vec(3);
		vector<VecDoub> p;
		for (Int i = 0; i < 1000; i++) {
			centre = testtdev.dev();
			for (Int i = 0; i < 3; i++) {
				vec[i] = centre[i];
			}
			p.push_back(vec);
		}
		points(p);
	}

	void ellipsoid(VecDoub ccentre, MatDoub ccov) {
		VecDoub centre(3);
		MatDoub cov(3, 3);
		for (int i = 0; i < 3; i++) {
			centre[i] = ccentre[i];
			for (int j = 0; j < 3; j++) {
				cov[i][j] = ccov[i][j] * 9.2103;
			}
		}
		Jacobi jac(cov);
		MatDoub evals(3, 3, 0.);
		for (int i = 0; i < 3; i++) {
			evals[i][i] = sqrt(jac.d[i]);
		}
		vector<VecDoub> points;
		VecDoub v(3);
		Int nPoints = 40;
		Doub angleStep = 2. * 3.14159265358979 / (Doub)(nPoints - 1);
		for (int k = 0; k < 3; k++) {
			for (int i = 0; i < nPoints; i++) {
				v[(k + 0) % 3] = cos( angleStep * (Doub)i);
				v[(k + 1) % 3] = sin( angleStep * (Doub)i);
				v[(k + 2) % 3] = 0;
				v = mXv(evals, v);
				v = mXv(jac.v, v);
				for (int i = 0; i < 3; i++) {
					v[i] += centre[i];
				}
				points.push_back(v);
			}
			line(points);
			points.clear();
		}
	}

    void surface(vector<VecDoub> points) {
        graphs.push_back(Graph(path + to_string(++gnuplotfile),  's', 1, 3).add(points));
    }
           

	void show() {
		string tocmd = "splot ";
		for (Int i = 0; i < graphs.size(); i++) {
			if (i) {
				tocmd += ", ";
			}
			tocmd += graphs.at(i).tocmd();
		}
		//tocmd += " w pm3d";
		gnuplotcmd(tocmd.c_str());
	}

   	void show2() {
        //gnuplotcmd("set style data lines");
        gnuplotcmd("set dgrid3d");
		gnuplotcmd("set pm3d corners2color c2");
		string tocmd = "splot ";
		for (Int i = 0; i < graphs.size(); i++) {
			if (i) {
				tocmd += ", ";
			}
			tocmd += graphs.at(i).tocmd();
		}
		tocmd += " w pm3d";
		gnuplotcmd(tocmd.c_str());
	}
};

struct Plot2D : Plot {
	void point(Doub x, Doub y, Int style = 2) {
		VecDoub v (2);
		v[0] = x;
		v[1] = y;
		vector<VecDoub> points;
		points.push_back(v);
		graphs.push_back(Graph(path + to_string(++gnuplotfile),  'p', style, 2).add(points) );
	}

	void point(VecDoub v,  Int style = 2) {
		vector<VecDoub> points;
		points.push_back(v);
		graphs.push_back(Graph(path + to_string(++gnuplotfile),  'p', style, 2).add(points) );
	}

	void points(vector<VecDoub> points,  Int style = 1) {
		graphs.push_back(Graph(path + to_string(++gnuplotfile),  'p', style, 2).add(points) );
	}

	void line(vector<VecDoub> points,  Int style = 3) {
		graphs.push_back(Graph(path + to_string(++gnuplotfile),  'l', style, 2).add(points) );
	}
	void ellipsoid(VecDoub centre, MatDoub cov) {
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				cov[i][j] *= 9.2103;
			}
		}
		Jacobi jac(cov);
		MatDoub evals(2, 2, 0.);
		for (int i = 0; i < 2; i++) {
			evals[i][i] = sqrt(jac.d[i]);
		}
		vector<VecDoub> points;
		VecDoub v(2);
		Int nPoints = 40;
		Doub angleStep = 2. * 3.14159265358979 / (Doub)(nPoints - 1);
		for (int i = 0; i < nPoints; i++) {
			v[0] = cos( angleStep * (Doub)i);
			v[1] = sin( angleStep * (Doub)i);
			v = mXv(evals, v);
			v = mXv(jac.v, v);
			for (int j = 0; j < 2; j++) {
				v[j] += centre[j];
			}
			points.push_back(v);		
		}
		point(centre, 3);
		line(points);
	}

	void show() {
		string tocmd = "plot ";
		for (Int i = 0; i < graphs.size(); i++) {
			if (i) {
				tocmd += ", ";
			}
			tocmd += graphs.at(i).tocmd();
		}
		gnuplotcmd(tocmd.c_str());
	}
};

vector<VecDoub> ps;
Plot2D p2d;
Plot2D eplt;

Plot3D pplt, vplt;
VecDoub vel(3), pos(3), ptrue(3), vtrue(3);
MatDoub vcov(3, 3), pcov(3, 3);


/*
void ellipsoid(VecDoub centre, MatDoub cov) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				cov[i][j] *= 1.;
			}
		}
		Jacobi jac(cov);

		MatDoub evals(3, 3, 0.);
		for (int i = 0; i < 3; i++) {
			evals[i][i] = jac.d[i];
		}

		vector<VecDoub> points;
		VecDoub v(3);
		
		Int nPoints = 30;
		Doub angleStep = 2. * 3.14159265358979 / (Doub)(nPoints - 1);
		for (int i = 0; i < nPoints; i++) {
			v[0] = cos( angleStep * (Doub)i);
			v[1] = sin( angleStep * (Doub)i);
			v[2] = 0;
			v = mXv(evals, v);
			v = mXv(jac.v, v);
			for (int i = 0; i < 3; i++) {
				v[i] += centre[i];
			}
			points.push_back(v);
		}
		line(points);
		
		points.clear();
		for (int i = 0; i < nPoints; i++) {
			v[0] = cos( angleStep * (Doub)i);
			v[2] = sin( angleStep * (Doub)i);
			v[1] = 0;
			v = mXv(evals, v);
			v = mXv(jac.v, v);
			for (int i = 0; i < 3; i++) {
				v[i] += centre[i];
			}
			points.push_back(v);
		}
		line(points);
		
		points.clear();
		for (int i = 0; i < nPoints; i++) {
			v[1] = cos( angleStep * (Doub)i);
			v[2] = sin( angleStep * (Doub)i);
			v[0] = 0;
			v = mXv(evals, v);
			v = mXv(jac.v, v);
			for (int i = 0; i < 3; i++) {
				v[i] += centre[i];
			}
			points.push_back(v);
		}
		line(points);
	}
	*/
