#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

/*
 *    This file is part of TEMPO2. 
 * 
 *    TEMPO2 is free software: you can redistribute it and/or modify 
 *    it under the terms of the GNU General Public License as published by 
 *    the Free Software Foundation, either version 3 of the License, or 
 *    (at your option) any later version. 
 *    TEMPO2 is distributed in the hope that it will be useful, 
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *    GNU General Public License for more details. 
 *    You should have received a copy of the GNU General Public License 
 *    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
 */

/*
 *    If you use TEMPO2 then please acknowledge it by citing 
 *    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
 *    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
 *    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
 *    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
 *    timing model.
 */

#include <stdio.h>
#include <math.h>
#include "tempo2.h"
#include "ifunc.h"

/* ------------------------------------------------------------------------- */
/*  Timing model for small-eccentricity binary pulsars, e<<1 (Wex 1998)      */
/*                                                                           */
/*  Instead of e and omega the Laplace parameters                            */
/*     epsilon1 = e*sin(omega)	                                             */
/*     epsilon2 = e*cos(omega)                                               */
/* are used as new parameters. T0 is related to the ascending node (not to   */
/* periastron as in BT, DD, ...)                                             */
/*                                                                           */
/*  Time derivatives:                                                        */
/*     nell1=0 -> fit for eps1dot,eps2dot                                    */
/*     nell1=1 -> fit for omdot,edot                                         */
/*                                                                           */
/*  Computes pulsar orbit time, torb, at time of observation t=ct(n)-pepoch. */
/*  Pulsar proper time is then TP=T+TORB.                                    */
/*  Units are such that c=G=1. Thus masses have units of seconds, with       */
/*  one solar mass = 4.925490947 usec.                                       */
/*                                                                           */
/*  Also computes the binary orbit-related values of fctn: partial           */
/*  derivatives of each arrival time residual with respect to the model      */
/*  parameters.                                                              */
/*                                                                           */
/*  Based on bnryell1.f                                                      */
/* ------------------------------------------------------------------------- */

/*  This is a modified ELL1 model.                                           */ 
/*  If the EDOT parameter is present in the par file, the Romer delay is     */
/*  calculated without the low-eccentricity approximation; otherwise, the    */
/*  original ELL1 model is used.                                             */


double ELL1model(pulsar *psr,int p,int ipos,int param,int k)
{
    double an,x0,m2,tt0,orbits,phase,e1,e2,dre,drep,drepp,brace,dlogbr,ds,da,pb;
    double eps1,eps2,eps1dot,eps2dot,si,a0,b0;
    double d2bar,torb,Csigma,Cx,Ceps1,Ceps2,Cm2,Csi,ct,t0asc,pbdot,xpbdot,x,xdot,am2;
    int norbits;
    double SUNMASS = 4.925490947e-6;
    const char *CVS_verNum = "$Id: 853536c7e64efde12efd9a9b936741f9eb9bec4a $";

    if (displayCVSversion == 1) CVSdisplayVersion("ELL1model.C","ELL1model()",CVS_verNum);

    a0 = 0.0; /* WHAT SHOULD THESE BE? */
    b0 = 0.0;

    if (psr[p].param[param_fb].paramSet[0]==1) 
        pb    = 1.0/psr[p].param[param_fb].val[0];
    else
        pb    = psr[p].param[param_pb].val[0]*SECDAY;

    if (psr[p].param[param_pbdot].paramSet[0] == 1) pbdot = psr[p].param[param_pbdot].val[0];
    else pbdot=0.0;

    an    = 2.0*M_PI/pb;
    if (psr[p].param[param_sini].paramSet[0]==1) si = getParameterValue(&psr[p],param_sini,0);
    else si = 0.0;

    if (si > 1.0)
    {
        displayMsg(1,"BIN1","SIN I > 1.0, setting to 1: should probably use DDS model","",psr[p].noWarnings);
        si = 1.0;
        psr[p].param[param_sini].val[0] = longdouble(1.0);
    }

    x0    = psr[p].param[param_a1].val[0];
    if (psr[p].param[param_a1dot].paramSet[0] == 1) 
        xdot  = psr[p].param[param_a1dot].val[0];
    else 
        xdot = 0.0;
    t0asc = psr[p].param[param_tasc].val[0];

    if (psr[p].param[param_m2].paramSet[0]==1) am2 = psr[p].param[param_m2].val[0];
    else am2 = 0.0;
    m2 = am2*SUNMASS;

    xpbdot = 0.0;
    eps1  = psr[p].param[param_eps1].val[0];
    eps2  = psr[p].param[param_eps2].val[0];
    if (psr[p].param[param_eps1dot].paramSet[0]==1) eps1dot = psr[p].param[param_eps1dot].val[0];
    else eps1dot=0;
    if (psr[p].param[param_eps2dot].paramSet[0]==1) eps2dot = psr[p].param[param_eps2dot].val[0];
    else eps2dot=0;

    ct = psr[p].obsn[ipos].bbat;      
    tt0 = (ct-t0asc)*SECDAY;
    // --- Changes to handle higher orbital-frequency derivatives (FB1, FB2, ...) ---                                               
    // 04/2015, H. J. Pletsch                                                                                                       
    orbits = tt0/pb;
    if (psr[p].param[param_fb].paramSet[1]==1) {
        int j;
        double fac = 1.0;
        for (j=1;j<psr[p].param[param_fb].aSize;j++) {
          double fbx;
            fac = fac/((double)(j+1));
            if (psr[p].param[param_fb].paramSet[j]==1) {
                fbx = psr[p].param[param_fb].val[j];
                orbits += fac * fbx * pow(tt0,j+1);
            }
        }
    } else {
        orbits -= 0.5*(pbdot+xpbdot)*pow(tt0/pb,2);
    }
      // --- End of changes to handle higher orbital-frequency derivatives ---                                                        

    if (psr[p].param[param_orbifunc].paramSet[0] == 1)
    {
        longdouble wi,t1;
        longdouble dt,speriod,tt;

        if (psr[p].param[param_ifunc].val[0] == 1) // Sinc interpolation
        {
            t1=longdouble(0.0);
            speriod = (longdouble)(psr[p].orbifuncT[1]-psr[p].orbifuncT[0]);
            //	       printf("ifuncN = %d\n",psr[p].ifuncN);
            for (k=0;k<psr[p].orbifuncN;k++)
            //	       for (k=3;k<4;k++)
            {
            //		   printf("Have %g %g\n",psr[p].ifuncT[k],psr[p].ifuncV[k]);
                dt = psr[p].obsn[ipos].bbat - (longdouble)psr[p].orbifuncT[k];
                wi=1;
                if (dt==0)
                {
                    t1 += wi*(longdouble)psr[p].orbifuncV[k];
                }
                else
                {
                    tt = M_PI/speriod*(dt);
                    t1 += wi*(longdouble)psr[p].orbifuncV[k]*sinl(tt)/(tt);
                    //		       t2 += wi*sinl(tt)/(tt);
                }
            }
            orbits += t1;
        }
        else if (psr[p].param[param_orbifunc].val[0] == 2) // Linear interpolation
        {
            double ival = ifunc(psr[p].orbifuncT,psr[p].orbifuncV,(double)psr[p].obsn[ipos].sat,psr[p].orbifuncN);
            orbits += ival;
        }
        else if (psr[p].param[param_orbifunc].val[0] == 0) // No interpolation
        {
            int k;
            double ival;
            for (k=0;k<psr[p].ifuncN-1;k++)
            {
                if ((double)psr[p].obsn[ipos].sat >= psr[p].orbifuncT[k])
                {
                    ival = psr[p].orbifuncV[k];
                    break;
                }
            }
            orbits += ival;
        }
    }

    norbits = (int)orbits;
    if (orbits<0.0) norbits = norbits-1;
    phase = 2.0*M_PI*(orbits-norbits);

    x = x0+xdot*tt0;

    e1 = eps1+eps1dot*tt0;
    e2 = eps2+eps2dot*tt0;
    
    if(psr[p].param[param_edot].paramSet[0] == 0){
    	dre  = x*(sin(phase)-0.5*(e1*cos(2.0*phase)-e2*sin(2.0*phase)));
    	drep = x*cos(phase);
    	drepp=-x*sin(phase);

        brace=1-si*sin(phase);
        dlogbr=log(brace);
    	ds=-2*m2*dlogbr;

    	/* NOTE: a0 and b0 are always zero -- they are not set in the original TEMPO!!!!! */
    	da=a0*sin(phase)+b0*cos(phase);  

    	/*  Now compute d2bar (cf. DD 52) */

    	d2bar=dre*(1-an*drep+pow(an*drep,2)+0.5*pow(an,2)*dre*drepp)+ds+da;
    	torb=-d2bar;
        
    	Csigma   = x*(cos(phase)+e1*sin(2.0*phase)+e2*cos(2*phase));
    	Cx       = sin(phase)-0.5*(e1*cos(2.0*phase)-e2*sin(2.0*phase));
    	Ceps1    = -0.5*x*cos(2*phase);
    	Ceps2    =  0.5*x*sin(2*phase);
    	Cm2      = -2*dlogbr;
    	Csi      = 2*m2*sin(phase)/brace;
    }


    if(psr[p].param[param_edot].paramSet[0] == 1){
        long double ecc,omega,u,du,su,cu,onemecu,cae,sae,ae,sw,cw,alpha,beta,bg,anhat,sqr1me2,cume;
        longdouble Comega,Ce;
        
        ecc=sqrt(e1*e1+e2*e2);
        if(ecc!=0){
            omega=atan2(e1,e2);
            if(omega<0) omega+=2*M_PI;
            phase=phase-omega; 
            if(phase<0) phase+=2*M_PI;
        }
        else{
            omega=0;
        }

        // eccz=sqrt(eps1*eps1+eps2*eps2);
        // if(eccz!=0) omz=atan2(eps1,eps2);
        // else omz=0;

        /*  Compute eccentric anomaly u by iterating Kepler's equation. */
        u=phase+ecc*sin(phase)*(1.0+ecc*cos(phase));
        do {
            du=(phase-(u-ecc*sin(u)))/(1.0-ecc*cos(u));
            u=u+du;
        } while (fabs(du)>1.0e-12);

        /*  DD equations 17b, 17c, 29, and 46 through 52 */
        
        su=sin(u);
        cu=cos(u);
        onemecu=1.0-ecc*cu;
        sw=sin(omega);
        cw=cos(omega);
        alpha=x*sw;
        beta=x*sqrt(1-pow(ecc,2))*cw;
        bg=beta;
        dre=alpha*(cu-ecc) + bg*su;
        drep=-alpha*su + bg*cu;
        drepp=-alpha*cu - bg*su;
        anhat=an/onemecu;

        /* DD equations 26, 27, 57: */
        sqr1me2=sqrt(1-pow(ecc,2));
        cume=cu-ecc;
        brace=onemecu-si*(sw*cume+sqr1me2*cw*su);
        // printf("GEORGE: si = %g, brace = %g\n",(double)si,(double)brace);
        dlogbr=log(brace);
        ds=-2*m2*dlogbr;
        // da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw);

        /*  Now compute d2bar, the orbital time correction in DD equation 42. */
        d2bar=dre*(1-anhat*drep+(pow(anhat,2))*(pow(drep,2) + 0.5*dre*drepp -
                    0.5*ecc*su*dre*drep/onemecu)) + ds;
        torb=-d2bar;
        phase=phase+omega;
        
        Comega=x*(cw*cume-sqr1me2*sw*su);
        Ce=su*Csigma-x*sw-ecc*x*cw*su/sqr1me2;
        Csigma   = x*(-sw*su+sqr1me2*cw*cu)/onemecu;
    	Cx       = sw*(cu-ecc)+sqr1me2*cw*su;
        if(ecc==0){
        	Ceps1    = -0.5*x*cos(2*phase);
    		Ceps2    =  0.5*x*sin(2*phase);
        }
    	else{
            Ceps1    = Ce*sin(omega)+(Comega-Csigma)*cos(omega)/ecc;
    		Ceps2    = Ce*cos(omega)-(Comega-Csigma)*sin(omega)/ecc;
        }
        Cm2=-2*dlogbr;
        Csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace;
    }

    if (param==-1) return torb;

    /* Now we need the partial derivatives. */
    // Csigma   = x*(cos(phase)+e1*sin(2.0*phase)+e2*cos(2*phase));
    // Cx       = sin(phase)-0.5*(e1*cos(2.0*phase)-e2*sin(2.0*phase));
    // Ceps1    = -0.5*x*cos(2*phase);
    // Ceps2    =  0.5*x*sin(2*phase);
    // Cm2      = -2*dlogbr;
    // Csi      = 2*m2*sin(phase)/brace; 

    if (param==param_pb)
        return -Csigma*an*SECDAY*tt0/(pb*SECDAY); /* Pb    */
    else if (param==param_fb)
    {
        tt0 /= 1.0e7; 
        double dfb = Csigma*an*tt0*pb;
        if (k==0) return dfb;
        for (int ik=1; ik<=k; ik++) {
            dfb *= tt0/(double)(ik+1);
        }
        return dfb;
    }
    else if (param==param_a1)
        return Cx;
    else if (param==param_eps1)
        return Ceps1;
    else if (param==param_tasc)
        return -Csigma*an*SECDAY;
    else if (param==param_eps2)
        return Ceps2;
    else if (param==param_eps1dot)
        return Ceps1*tt0;
    else if (param==param_eps2dot)
        return Ceps2*tt0;
    else if (param==param_pbdot)
        return 0.5*tt0*(-Csigma*an*SECDAY*tt0/(pb*SECDAY));
    else if (param==param_a1dot)
        return Cx*tt0;  
    else if (param==param_sini)
        return Csi;
    else if (param==param_m2)
        return Cm2*SUNMASS;

    return 0.0;
}

void updateELL1(pulsar *psr,double val,double err,int pos,int k)
{
    if (pos==param_pb)
    {
        psr->param[param_pb].val[0] += val/SECDAY;
        psr->param[param_pb].err[0]  = err/SECDAY;
    }
    else if (pos==param_fb)
    {
        psr->param[param_fb].val[k] += (val/powl(1.0e7,k+1));
        psr->param[param_fb].err[k]  = err/powl(1.0e7,k+1);
    }
    else if (pos==param_a1 || pos==param_eps1 || pos==param_eps2 || pos==param_tasc
            || pos==param_sini || pos == param_m2 
            || pos==param_eps1dot || pos==param_eps2dot
            || pos==param_a1dot)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    else if (pos==param_pbdot)
    {
        psr->param[pos].val[0] += val;
        psr->param[pos].err[0]  = err;
    }
    else if (pos==param_omdot)
    {
        psr->param[pos].val[0] += val*(SECDAY*365.25)*180.0/M_PI;
        psr->param[pos].err[0]  = err*(SECDAY*365.25)*180.0/M_PI;
    }
}
