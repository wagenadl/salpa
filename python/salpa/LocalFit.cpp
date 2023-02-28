/* artifilt/LocalFit.C: part of meabench, an MEA recording and analysis tool
** Copyright (C) 2000-2003  Daniel Wagenaar (wagenaar@caltech.edu)
**
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

// LocalFit.cpp

//#define TEST 1
#include "LocalFit.h"
#include <iostream>

#define ASYM_NOT_CHI2 1
#ifndef THIRDORDER
#define THIRDORDER 1
#endif
#define PREMATURE 1

#include <math.h>

// static float LocalFit_X=0;
// static float LocalFit_Y=0;
// raw_t LocalFit::ZERO       = LocalFit_X/LocalFit_Y; // NAN;

Error::Error(char const *issuer0, char const *cause0) {
  issuer=issuer0;
  cause=cause0;
}

Error::~Error() {
}

void Error::report(char const *src) const {
  if (src)
    fprintf(stderr,"%s: %s: %s\n",src,issuer?issuer:"",cause?cause:"");
  else
    fprintf(stderr,"%s: %s\n",issuer?issuer:"",cause?cause:"");
}

const timeref_t INFTY = (timeref_t)(-1);

bool _sanity() {
  timeref_t ttest = INFTY;
  ttest++;
  if (ttest!=0) {
    printf("Warning! INFTY does not work on your platform. Please contact\nthe author, Daniel Wagenaar.\nFor contact information please visit http://www.danielwagenaar.net.\n");
    return false;
  }
  if (sizeof(timeref_t)!=8) {
    printf("Warning! I could not create an 8-bit integer type. Please contact\nthe author, Daniel Wagenaar.\nFor contact information please visit http://www.danielwagenaar.net.\n");
    return false;
  }
  return true;
}


//--------------------------------------------------------------------
// inline functions
//
inline void LocalFit::update_X012() {
  double y_new = source[t_stream+tau];
  double y_old = source[t_stream-tau-1];
  X0 += y_new - y_old;
  X1 += tau_plus_1*y_new - minus_tau*y_old - X0;
  X2 += tau_plus_1_squared*y_new - minus_tau_squared*y_old - X0 - 2*X1;
}

inline void LocalFit::calc_alpha0() {
  alpha0 = double(T4*X0 - T2*X2)/double(T0*T4-T2*T2);
}

//--------------------------------------------------------------------
// Other LocalFit methods
//
LocalFit::LocalFit(float const *source0, float *dest0,
                   timeref_t t_end0,
		   float threshold0, int tau0):
  source(source0), dest(dest0), y_threshold(threshold0), tau(tau0),
  t_blankdepeg(BLANKDEP), t_ahead(AHEAD),
  t_chi2(TCHI2), 
  rail1(RAIL1), rail2(RAIL2),
  t_end(t_end0) {
  init_T();
  reset();
}

void LocalFit::set_t_blankdepeg(int t_blankdepeg1) {
  t_blankdepeg = t_blankdepeg1;
}

void LocalFit::set_t_ahead(int t_ahead1) {
  t_ahead = t_ahead1;
}

void LocalFit::set_t_chi2(int t_chi21) {
  t_chi2 = t_chi21;
}


bool LocalFit::isValid() {
  return _sanity();
}

void LocalFit::reset(timeref_t t_start) {
  t_peg = t_stream = t_start;
  state=PEGGED;
}

void LocalFit::setrail(raw_t rail11, raw_t rail21) {
  rail1 = rail11;
  rail2 = rail21;
}

void LocalFit::init_T() {
#if ASYM_NOT_CHI2
  my_thresh = 3.92 * t_chi2 * y_threshold*y_threshold; // 95% conf limit
#else
  my_thresh = (t_chi2-4) * y_threshold*y_threshold;
#endif
  
  tau_plus_1 = tau+1;
  tau_plus_1_squared = tau_plus_1 * tau_plus_1;
  tau_plus_1_cubed = tau_plus_1_squared * tau_plus_1;
  minus_tau = -tau;
  minus_tau_squared = minus_tau * minus_tau;
  minus_tau_cubed = minus_tau_squared * minus_tau;
    
  T0=T2=T4=T6=0;
  for (int t=-tau; t<=tau; t++) {
    int_t t2=t*t;
    int_t t4=t2*t2;
    int_t t6=t4*t2;
    T0+=1;
    T2+=t2;
    T4+=t4;
    T6+=t6;
  }
}

timeref_t LocalFit::process(timeref_t t_limit) {
  state=statemachine(t_limit, state);
  return t_stream;
}

timeref_t LocalFit::forcepeg(timeref_t t_from, timeref_t t_to) {
  state=statemachine(t_from-tau, state);
  if (state==OK) {
    // goto state PEGGING
      t0=t_stream-1;
      calc_X3();
      calc_alpha0123();
      statemachine(t_from, PEGGING);
  }
  t0=t_to;
  state=statemachine(t_to, FORCEPEG);
  return t_stream;
}

LocalFit::State LocalFit::statemachine(timeref_t t_limit, State s) {

  /* This is a straightforward implementation of the statemachine I
     described on 9/9/01.
   * //// mark boundaries program flow does not pass through.
   * On exit, t_stream == t_limit.
  */
  //  timeref_t t_check=0; //DBG
  //  fprintf(stderr, "statemachine %Li %Li %i\n", t0, t_limit, s);
  switch (s) {
  case OK: goto l_OK;
  case PEGGED: goto l_PEGGED;
  case PEGGING: goto l_PEGGING;
  case TOOPOOR: goto l_TOOPOOR;
  case DEPEGGING: goto l_DEPEGGING;
  case FORCEPEG: goto l_FORCEPEG;
  case BLANKDEPEG: goto l_BLANKDEPEG;
  default: throw Error("LocalFit","BUG! Bad State");
  }

//////////////////////////////////////////////////
 l_PEGGED: {
    if (t_stream>=t_limit)
      return PEGGED;
    // fprintf(stderr, "PEGGED %Li %g %i\n", t_stream, source[t_stream],
    //    ispegged(source[t_stream]));
    if (ispegged(source[t_stream])) {
      dest[t_stream] = ZERO;
      t_stream++;
      goto l_PEGGED;
    }
    for (int dt=1; dt<=2*tau; dt++)
      if (t_stream+dt>=t_end || ispegged(source[t_stream+dt])) {
	t0 = t_stream+dt;
	goto l_FORCEPEG;
      }
    t0 = t_stream + tau;
    calc_X012(); calc_X3();
    calc_alpha0123();
    toopoorcnt=TOOPOORCNT;
#if TEST
    t_depeg = t_stream;
#endif
    goto l_TOOPOOR;
  }

  throw Error("LocalFit","Code breach");
  
//////////////////////////////////////////////////
 l_TOOPOOR: {
    if (t_stream>=t_limit) 
      return TOOPOOR;
#if ASYM_NOT_CHI2
    double asym=0;
    double sig=0;
    for (int i=0; i<t_chi2; i++) {
      timeref_t t_i = t_stream+i;
      if (t_i > t_end) {
	t0 = t_i;
	goto l_FORCEPEG;
      }
      int dt = int(t_i - t0);
      int dt2 = dt*dt;
      int dt3 = dt*dt2;
      double dy = alpha0 + alpha1*dt + alpha2*dt2 + alpha3*dt3 - source[t_i];
      asym += dy;
      sig += dy*dy;
    }
    //    fprintf(stderr,"TOOPOOR: t=%.2f t0=%.2f asym/sqrt(t)=%g [t_chi2=%i sqrt(sig/t)=%g]\n",
    //	    t_stream/25.,t0/25.,asym/sqrt(t_chi2+0.),t_chi2,sqrt(sig/t_chi2));
    //    fprintf(stderr,"    alpha = %g  %g  %g  %g\n",alpha0,alpha1,alpha2,alpha3);
    asym*=asym;
    if (asym<my_thresh)
      toopoorcnt--;
    else
      toopoorcnt=TOOPOORCNT;
    if (toopoorcnt<=0 && asym < my_thresh/3.92) {
#if PREMATURE
      int dt = int(t_stream - t0);
      int dt2 = dt*dt;
      int dt3 = dt*dt2;
      negv =  source[t_stream]
	< float(alpha0 + alpha1*dt + alpha2*dt2 + alpha3*dt3);
#endif

#ifdef TEST
      double x0=X0,x1=X1,x2=X2,x3=X3;
#endif
      calc_X012(); calc_X3(); // for numerical stability problem!
#ifdef TEST
      if (x0!=X0 || x1!=X1 || x2!=X2 || x3!=X3)
	fprintf(stderr,"CUMULANT ERROR: X=[%g %g %g %g] dx=[%g %g %g %g] (t=%.2f, dt=%Li, C=%Li)\n",
		X0,X1,X2,X3,X0-x0,X1-x1,X2-x2,X3-x3,
                t_stream/25.0,t_stream-t_depeg,t_stream-t_check);
      t_check=t_stream;
#endif
      goto l_BLANKDEPEG;
    }
#else
    double chi2=0;
    for (int i=0; i<t_chi2; i++) {
      int t_i = t_stream+t_blankdepeg+i;
      if (t_i > t_end) {
	t0 = t_i;
	goto l_FORCEPEG;
      }
      int dt = int(t_i - t0);
      int dt2 = dt*dt;
      int dt3 = dt*dt2;
      double dy = alpha0 + alpha1*dt + alpha2*dt2 + alpha3*dt3 - source[t_i];
      chi2+=dy*dy;
    }
    ////    fprintf(stderr,"TOOPOOR: chi2=%g [%2f-%2f]\n",chi2,
    ////	    (t_stream+t_blankdepeg)/25.,(t_stream+t_blankdepeg+TOOPOORCNT)/25.);
    if (chi2 < my_thresh) {
#if PREMATURE
      int dt = int(t_stream - t0);
      int dt2 = dt*dt;
      int dt3 = dt*dt2;
      negv = source[t_stream]
	< float(alpha0 + alpha1*dt + alpha2*dt2 + alpha3*dt3);
#endif
      goto l_BLANKDEPEG;
    }
#endif

//    int dt=t_stream-t0; int dt2=dt*dt; int dt3=dt*dt2;
//    dest[t_stream]=source[t_stream]-int(alpha0 + alpha1*dt + alpha2*dt2 + alpha3*dt3);
    dest[t_stream] = ZERO;
    t_stream++; t0++;
    if (t0+tau>=t_end || ispegged(source[t0+tau])) {
      t0=t0+tau;
      goto l_FORCEPEG;
    }
    update_X0123();
#ifdef TEST
    double x0=X0,x1=X1,x2=X2,x3=X3;
#endif
    calc_X012(); calc_X3(); // for numerical stability problem!
#ifdef TEST
    if (x0!=X0 || x1!=X1 || x2!=X2 || x3!=X3)
      fprintf(stderr,"Cumulant error: X=[%g %g %g %g] dx=[%g %g %g %g] (t=%.2f, dt=%Li, C=%Li)\n",
	      X0,X1,X2,X3,X0-x0,X1-x1,X2-x2,X3-x3,
              t_stream/25.0,t_stream-t_depeg,t_stream-t_check);
    t_check=t_stream;
#endif
    calc_alpha0123();
    goto l_TOOPOOR;
  }

  throw Error("LocalFit","Code breach");

//////////////////////////////////////////////////
 l_FORCEPEG: {
    if (t_stream>=t_limit)
      return FORCEPEG;
    if (t_stream>=t0)
      goto l_PEGGED;
    dest[t_stream] = ZERO;
    t_stream++;
    goto l_FORCEPEG;
  }

  throw Error("LocalFit","Code breach");

//////////////////////////////////////////////////
 l_BLANKDEPEG: {
   // fprintf(stderr, "BLANKDEPEG %Li %Li %i %i %Li\n", t_stream, t0, tau, t_blankdepeg, t_limit);
    if (t_stream>=t_limit)
      return BLANKDEPEG;
    if (t_stream >= t0-tau+t_blankdepeg)
      goto l_DEPEGGING;
#if PREMATURE
    int dt = int(t_stream-t0);
    int dt2 = dt*dt;
    int dt3 = dt*dt2;
    float y = source[t_stream]
      - float(alpha0 + alpha1*dt + alpha2*dt2 + alpha3*dt3);
    if ((y<0) != negv) {
      dest[t_stream] = y;
      t_stream++;
      goto l_DEPEGGING;
    }
#endif
    dest[t_stream] = ZERO;
    t_stream++;
    goto l_BLANKDEPEG;
  }

  throw Error("LocalFit","Code breach");

//////////////////////////////////////////////////
 l_DEPEGGING: {
   // fprintf(stderr, "DEPEGGING %Li %Li\n", t_stream, t0);
    if (t_stream>=t_limit)
      return DEPEGGING;
#if TEST
    if (t_depeg) {
      fprintf(stderr,"%.5f %.2f\n",
	      t_stream/(1000.*25),
	      (t_stream-t_depeg)/(1.*25));
      // fprintf(stderr,"Lost time: %.2f ms on %i\n",(t_stream-t_depeg)/25.,source.hw);
      t_depeg=0;
    }
#endif
    if (t_stream==t0) {
      goto l_OK;
    }
    int dt = int(t_stream-t0);
    int dt2 = dt*dt;
    int dt3 = dt*dt2;
    dest[t_stream] = source[t_stream]
      - float(alpha0 + alpha1*dt + alpha2*dt2 + alpha3*dt3);
    t_stream++;
    goto l_DEPEGGING;
  }

  throw Error("LocalFit","Code breach");

//////////////////////////////////////////////////
 l_PEGGING: {
    if (t_stream>=t_limit)
      return PEGGING;
    if (t_stream >= t0+tau) {
      t_peg = t_stream;
      goto l_PEGGED;
    }
    int dt = int(t_stream-t0);
    int dt2 = dt*dt;
    int dt3 = dt*dt2;
    dest[t_stream] = source[t_stream]
      - float(alpha0 + alpha1*dt + alpha2*dt2 + alpha3*dt3);
    t_stream++;
    goto l_PEGGING;
  }

  throw Error("LocalFit","Code breach");

//////////////////////////////////////////////////
 l_OK: {
    if (t_stream>=t_limit)
      return OK;
    calc_alpha0();
    //    fprintf(stderr, "OK %Li %g %g\n", t_stream, source[t_stream], alpha0);
    dest[t_stream] = source[t_stream] - float(alpha0);
    t_stream++;
    if (t_stream+tau+t_ahead>=t_end ||
	ispegged(source[t_stream+tau+t_ahead])) {
      t0=t_stream-1;
      calc_X3();
      calc_alpha0123();
      goto l_PEGGING;
    } 
    update_X012();
    goto l_OK;
  }

  throw Error("LocalFit","Code breach");

//////////////////////////////////////////////////
}

void LocalFit::calc_X012() {
  X0=X1=X2=0;
  for (int t=-tau; t<=tau; t++) {
    int_t t2=t*t;
    double y=source[t0+t];
    X0+=y;
    X1+=t*y;
    X2+=t2*y;
  }
}

void LocalFit::calc_X3() {
  X3=0;
  for (int t=-tau; t<=tau; t++) {
    int_t t3=t*t*t;
    double y=source[t0+t];
    X3+=t3*y;
  }
}

void LocalFit::update_X0123() {
  // based on calc dd 9/7/01 - do recheck!
  double y_new = source[t0+tau];
  double y_old = source[t0-tau-1];
  X0 += y_new - y_old;
  X1 += tau_plus_1*y_new - minus_tau*y_old - X0;
  X2 += tau_plus_1_squared*y_new - minus_tau_squared*y_old - X0 - 2*X1;
  X3 += tau_plus_1_cubed*y_new - minus_tau_cubed*y_old - X0 - 3*X1 - 3*X2;
}

void LocalFit::calc_alpha0123() {
  double fact02 = 1./(T0*T4-T2*T2);
  alpha0 = fact02*(T4*X0 - T2*X2);
  alpha2 = fact02*(T0*X2 - T2*X0);
#if THIRDORDER
  double fact13 = 1./(T2*T6-T4*T4);
  alpha1 = fact13*(T6*X1 - T4*X3);
  alpha3 = fact13*(T2*X3 - T4*X1);
#else
  alpha1 = double(X1)/T2;
  alpha3 = 0;
#endif

  ////  report();
}

//--------------------------------------------------------------------
// debug
//
#include <stdio.h>

void LocalFit::report() {
  std::cerr << "state="
	    << (state==OK?"OK":
		state==PEGGING?"PEGGING":
		state==PEGGED?"PEGGED":
		state==TOOPOOR?"TOOPOOR":
		state==DEPEGGING?"DEPGGING":
		state==FORCEPEG?"FORCEPEG":
		state==BLANKDEPEG?"BLANKDEP":
		"???")
	    << " ";
  std::cerr << "t_stream=" << t_stream/25.0
	    << " t0=" << t0/25.0
	    << " y[t]=" <<source[t_stream]
	    << " alpha=" << alpha0 << " " << alpha1 << " " << alpha2
	    << " " << alpha3
	    << " X=" << X0 << " " << X1 << " " << X2 << " " << X3
	    << "\n";
}

void LocalFit::inirep() {
  std::cerr << "tau="<<tau<<"\n"; 
  std::cerr << "T0="<<T0<<"\n"; 
  std::cerr << "T2="<<T2<<"\n"; 
  std::cerr << "T4="<<T4<<"\n"; 
  std::cerr << "T6="<<T6<<"\n"; 
}
