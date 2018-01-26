
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>

#include "nr3/nr3.h"
#include "nr3/cholesky.h"
#include "nr3/ran.h"
#include "nr3/svd.h"
#include "nr3/multinormaldev.h"
#include "nr3/ludcmp.h"
#include "nr3/eigen_sym.h"
#include "nr3/mins.h"
#include "nradd.h"

#include "prognosis/Runge_good.h"


#include "Date.h"
#include "VecEst.h"
#include "Sensor.h"
#include "Observs.h"
#include "prognosis.h"

#include "gnuplot.h"
vector <VecDoub> pov;

Vec vec1;
Vec vec2;
Est paramsEst;

#include "io.h"
#include "simulation.h"
#include "prelim.h"
#include "common.h"


int main(int argc, char* argv[]) {
	printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INIT\n");
	prognosisInit();
	Sensor tel = initTelescope("init/telescope.dat");
	vec1 = initVec("init/high.dat");
	vec2 = preciseExt(vec1, vec1.jd + 3. * isOrbit(vec1));
	printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END INIT\n\n");
	printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SIMULATION\n");
	Observs track1 = modelTrack(vec1, tel);
	vec1.print();
	printf("First track: %lf\n", residual(vec1, track1));
	Observs track2 = modelTrack(vec2, tel);
	vec2.print();
	printf("Second track: %lf\n", residual(vec2, track2));
	printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END SIMULATION\n\n");

    
    for (int i = 0; i < 3; i++) {
        ptrue[i] = vec1.vec[i];
		vtrue[i] = vec1.vec[i+3];
	}
	vplt.point(vtrue, 5);
	pplt.point(ptrue, 4);

    VecDoub xyzprecise(3);
    VecDoub xyzsmooth(3);
    vector<VecDoub> precise;
    vector<VecDoub> smooth;


	VecDoub evrav = VecDoub(3);
    for (Doub d = 0.1; d < 80; d += 3.) {
        Doub ra = track1.obs.front().vec[0], dec = track1.obs.front().vec[1];
	    VecDoub icrfTel = track1.sen.icrfPos(track1.obs.front().jd);
	    VecDoub x0 = radec2xyz(ra, dec);
	    for (int i = 0; i < 3; i++) {
		    x0[i] = x0[i] * d + icrfTel[i];
	    }
	    VecDoub xn, g = gravityVector(x0);

	    Doub dt = track1.obs[1].jd - track1.obs.front().jd, dn;
	    icrfTel = track1.sen.icrfPos(track1.obs[1].jd);
	    xn = radec2xyz(track1.obs[1].vec[0], track1.obs[1].vec[1]);
	    dn = 0.1;
	    for (int i = 0; i < 3; i++) {
		    xn[i] = xn[i] * dn + icrfTel[i];
		    evrav[i] = (xn[i] - x0[i]) / dt - g[i] * dt / 2.;
	    }
	    xn = radec2xyz(track1.obs[1].vec[0], track1.obs[1].vec[1]);
	    dn = 100;
	    for (int i = 0; i < 3; i++) {
		    xn[i] = xn[i] * dn + icrfTel[i];
		    evrav[i] -= (xn[i] - x0[i]) / dt - g[i] * dt / 2.;	//evrav[i] /= 1000.;
	    }

        Est orb = prelimWithFixedD(track1, d);
	    vector<VecDoub> vlim = fixeDVelRavine(orb, evrav);
       	if (vlim.size() != 2) {
    		continue;
    	}

        for (int i = 0; i < 3; i++) {
            pos[i] = orb.vec[i];
		    vel[i] = orb.vec[i+3];
	    }
        pplt.point(pos);
        pplt.show();
        vplt.point(vel);
		vplt.line(vlim);
		vplt.show();
  
        VecDoub vmin(3), vmax(3);

        for (int i = 0; i < 3; i++) {
			evrav[i] = vlim[1][i] - vlim[0][i];
			vel[i] = orb.vec[i+3];// = vlim[0][i] + (vlim[1][i] - vlim[0][i]) / 2.;
			vmin[i] = vlim[0][i] - orb.vec[i+3];
			vmax[i] = vlim[1][i] - orb.vec[i+3];
		}
		//vplt.point(vel, 3);
		evrav = normilize(evrav);
		Doub scalarvmin = innerProduct(vmin, evrav);
		Doub scalarvmax = innerProduct(vmax, evrav);

        VecDoub de = radec2xyz(track2.obs.front().vec[0], track2.obs.front().vec[1]);
		VecDoub telICRF = track2.sen.icrfPos(track2.obs.front().jd);
        Doub k = 398.6004415;
        for(Doub v=scalarvmin; v < scalarvmax; v += 1. ) {
            Est tmp(orb);
            for (int i = 0; i < 3; i++) { 
			    tmp.vec[i+3] = orb.vec[i+3] + v * evrav[i];
		    }

             if (!(isOrbit(tmp) > 0.)) {
                continue;
            }
		    tmp = preciseExt(tmp, track2.obs.front().jd);
            xyzprecise[0] = d;
            xyzprecise[1] = v;
            xyzprecise[2] = residual(tmp, track2);
            precise.push_back(xyzprecise);
		
		    VecDoub ve(3), re(3), lamda(3);
		    for (int i = 0; i < 3; i++) {
			    re[i] = tmp.vec[i];
			    ve[i] = tmp.vec[i+3];
		    }
		    Doub reMod = scalar(re);
		    VecDoub sigma = crossProduct(re, ve);
		    Doub sigmaMod = scalar(sigma);
		    VecDoub sigmaXve = crossProduct(sigma, ve);
		    for (int i = 0; i < 3; i++) {
			    lamda[i] = -sigmaXve[i] - k / reMod * re[i];
		    }

		    VecDoub r(3, 0.), vel(3);
		    VecDoub sXl = crossProduct(telICRF, lamda);
		    VecDoub deXl = crossProduct(de, lamda);
		
		    MatDoub abc(2, 2);
		    VecDoub	in(2), out(2);
		    for (int i = 0; i < 2; i++) {
			    abc[i][0] = sigma[i+1];
			    abc[i][1] = deXl[i+1];
			    in[i] = -sXl[i+1];
		    }
		    SVD svd(abc);
		    svd.solve(in, out);
		    Doub d2 = out[1];
		    for (int i = 0; i < 3; i++) {
			    r[i] = telICRF[i] + de[i] * d2;
		    }
		    VecDoub p = normilize(r);
		    VecDoub n = normilize(crossProduct(sigma, p));

		    Doub vr = -innerProduct(lamda, n) / sigmaMod;//  scalar(crossProduct(lamda, p)) / sigmaMod;// out[0] / rMod;//-lamdaMod / sigmaMod;
		    Doub vn = -innerProduct(lamda, n) / out[0];
		
		    //vn = sqrt(h + 2 * k / rMod - vr*vr);
		    //Doub rMod = 2 * k / (vr*vr +vn*vn - h);
		    //Doub rMod = out[0] / vr;
		    //vr = out[0] / rMod;
		    //vn = sigmaMod / rMod;
		    //Doub rMod = (k + innerProduct(lamda, p)) / vn / vn;
		    //Doub rMod = (-k - sqrt(k*k + (h-vn*vn) * out[0]))/(h-vn*vn);
		    //vr = sqrt(h + 2 * k / rMod - vn*vn);
		    //rMod = out[0] / vr;
		    //vr = out[0] / rMod;
		    //rMod = out[0] / vr;
		    for (int i = 0; i < 3; i++) {
			    //r[i] = p[i] * rMod;
			    vel[i] = p[i] * vr + n[i] * vn;
			    tmp.vec[i] = r[i];
			    tmp.vec[i+3] = vel[i];
		    }
		    //Observs track3(track2.sen);
		    //track3.obs.push_back(track2.obs.back());
            xyzsmooth[0] = d;
            xyzsmooth[1] = v;
            xyzsmooth[2] = residual(tmp, track2);
            smooth.push_back(xyzsmooth);
        }

    }


	Doub trued = dOrb(vec1, tel);
    VecDoub truepoint(3);
    truepoint[0] = trued;
    Doub ra = track1.obs.front().vec[0], dec = track1.obs.front().vec[1];
	VecDoub icrfTel = track1.sen.icrfPos(track1.obs.front().jd);
    VecDoub x0 = radec2xyz(ra, dec);
    for (int i = 0; i < 3; i++) {
	    x0[i] = x0[i] * trued + icrfTel[i];
    }
    VecDoub xn, g = gravityVector(x0);
	
    Doub dt = track1.obs[1].jd - track1.obs.front().jd, dn;
    icrfTel = track1.sen.icrfPos(track1.obs[1].jd);
    xn = radec2xyz(track1.obs[1].vec[0], track1.obs[1].vec[1]);
    dn = 0.1;
    for (int i = 0; i < 3; i++) {
	    xn[i] = xn[i] * dn + icrfTel[i];
	    evrav[i] = (xn[i] - x0[i]) / dt - g[i] * dt / 2.;
    }
    xn = radec2xyz(track1.obs[1].vec[0], track1.obs[1].vec[1]);
    dn = 100;
    for (int i = 0; i < 3; i++) {
	    xn[i] = xn[i] * dn + icrfTel[i];
	    evrav[i] -= (xn[i] - x0[i]) / dt - g[i] * dt / 2.;	//evrav[i] /= 1000.;
    }
    Est orb = prelimWithFixedD(track1, trued);
    vector<VecDoub> vlim = fixeDVelRavine(orb, evrav);

    for (int i = 0; i < 3; i++) {
		evrav[i] = vlim[1][i] - vlim[0][i];
		vel[i] = orb.vec[i+3];// = vlim[0][i] + (vlim[1][i] - vlim[0][i]) / 2.;
	}

    evrav = normilize(evrav);


    Doub truev =  innerProduct(vtrue, evrav) - innerProduct(vel, evrav);
    Vec tmp(orb);
    for (int i = 0; i < 3; i++) { 
	    tmp.vec[i+3] = orb.vec[i+3] + truev * evrav[i];
	}
    print(tmp);
    truepoint[1] = truev;
    truepoint[2] = 0.;
    Plot3D sm;
    Plot3D pr;
    sm.point(truepoint);
    sm.surface(smooth);

    pr.point(truepoint);    
    pr.surface(precise);

    pr.show2();
    sm.show2();
    getchar();


	
}
