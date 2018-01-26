
MatDoub mXm(MatDoub a, MatDoub b){
	MatDoub c(a.nrows(), b.ncols(), 0.);
    Int i, j, k;
    for (i = 0; i < a.nrows(); i++) {
        for (j = 0; j < b.ncols(); j++) {
            for (k = 0; k < a.ncols(); k++) {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }
    return c;
}

VecDoub mXv(MatDoub& m, VecDoub& v) {
	VecDoub out(v.size(), 0.);
    Int i, j;
    for (i = 0; i < m.nrows(); i++) {
        for (j = 0; j < m.ncols(); j++) {
			out[i] += m[i][j] * v[j];
        }
    }
    return out;
}

MatDoub trans(MatDoub& a) {
    int i, j;
    MatDoub r (a.ncols(), a.nrows());
    for (i = 0; i < a.nrows(); i++) {
        for (j = 0; j < a.ncols(); j++) {
            r[j][i] = a[i][j];
        }
    }
    return r;
}

VecDoub crossProduct(VecDoub& a, VecDoub& b) {
	VecDoub c(3);
	c[0] = a[2]*b[1] - a[1]*b[2];
	c[1] = a[0]*b[2] - a[2]*b[0];
	c[2] = a[1]*b[0] - a[0]*b[1];
	return c;
}

Doub innerProduct(VecDoub& a, VecDoub& b) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


Doub innerProduct(VecDoub& a) {
	return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
}

Doub module(VecDoub v) {
	Doub module = 0.;
	for (Int i = 0; i < v.size(); i++) {
		module += v[i] * v[i];
	}
	return sqrt(module);
}

VecDoub normilize(VecDoub v) {
	Doub len = 0.;
	for (Int i = 0; i < v.size(); i++) {
		len += v[i] * v[i];
	}
	len = sqrt(len);
	for (Int i = 0; i < v.size(); i++) {
		v[i] /= len;
	}
	return v;
}

Doub scalar(VecDoub& v) {
	Doub r = 0.;
	for (Int i = 0; i < v.size(); i++) {
		r += v[i] * v[i];
	}
	r = sqrt(r);
	return r;
}


MatDoub iMat(Int n) {
	MatDoub e(n, n, 0.);
	for (Int i = 0; i < n; i++) {
		e[i][i] = 1.;
	}
	return e;
}

VecDoub matXVec(MatDoub& m, VecDoub& v) {
	VecDoub r(v.size(), 0.);
	for (Int i = 0; i < m.nrows(); i++) {
		for (Int j = 0; j < m.ncols(); j++) {
			r[i] += m[i][j] * v[j];
		}
	}
	return r;
}

MatDoub projectionMatrix(VecDoub& v) {
	MatDoub p(v.size(), v.size(), 0.);
	for (Int i = 0; i < v.size(); i++) {
		for (Int j = 0; j < v.size(); j++) {
			p[i][j] = (Doub)(Int)(i==j) - v[i] * v[j];
		}
	}
	return p;
}

//MatDoub AXBXAt(MatDoub& a, MatDoub& b) {
//}
