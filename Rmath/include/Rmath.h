/* -*- C -*-
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998-2022  The R Core Team
 *  Copyright (C) 2004       The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *

 * Rmath.h  should contain ALL headers from R's C code in `src/nmath'
   -------  such that ``the Math library'' can be used by simply

   ``#include <Rmath.h> ''

   and nothing else.

   It is part of the API and supports 'standalone Rmath'.
   Some entries possibly are not yet documented in 'Writing R Extensions'.

*/
#ifndef RMATH_H
#define RMATH_H

/* needed for cospi etc */
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
# define __STDC_WANT_IEC_60559_FUNCS_EXT__ 1
#endif

#if defined(__cplusplus) && !defined(DO_NOT_USE_CXX_HEADERS)
# include <cmath>
#else
# include <math.h>
#endif

#ifdef NO_C_HEADERS
# warning "use of NO_C_HEADERS is defunct and will be ignored"
#endif

/*-- Mathlib as part of R --  define this for standalone : */
/* #undef MATHLIB_STANDALONE */

#define R_VERSION_STRING "4.2.2"

// Legacy defines -- C99 functions which R >= 3.5.0 requires
#ifndef HAVE_EXPM1
# define HAVE_EXPM1 1
#endif
#ifndef HAVE_HYPOT
# define HAVE_HYPOT 1
#endif
#ifndef HAVE_LOG1P
# define HAVE_LOG1P 1
#endif

#ifndef HAVE_WORKING_LOG1P
# define HAVE_WORKING_LOG1P 1
#endif

#if !defined(HAVE_WORKING_LOG1P)
/* remap to avoid problems with getting the right entry point */
double  Rlog1p(double);
#define log1p Rlog1p
#endif


/* ----- The following constants and entry points are part of the R API ---- */

/* 30 Decimal-place constants */
/* Computed with bc -l (scale=32; proper round) */

/* SVID & X/Open Constants */
/* Names from Solaris math.h */

#ifndef M_E
#define M_E		2.718281828459045235360287471353	/* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E		1.442695040888963407359924681002	/* log2(e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E	0.434294481903251827651128918917	/* log10(e) */
#endif

#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif

#ifndef M_LN10
#define M_LN10		2.302585092994045684017991454684	/* ln(10) */
#endif

#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#ifndef M_PI_2
#define M_PI_2		1.570796326794896619231321691640	/* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4		0.785398163397448309615660845820	/* pi/4 */
#endif

#ifndef M_1_PI
#define M_1_PI		0.318309886183790671537767526745	/* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI		0.636619772367581343075535053490	/* 2/pi */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI	1.128379167095512573896158903122	/* 2/sqrt(pi) */
#endif

#ifndef M_SQRT2
#define M_SQRT2		1.414213562373095048801688724210	/* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2	0.707106781186547524400844362105	/* 1/sqrt(2) */
#endif

/* R-Specific Constants */

#ifndef M_SQRT_3
#define M_SQRT_3	1.732050807568877293527446341506	/* sqrt(3) */
#endif

#ifndef M_SQRT_32
#define M_SQRT_32	5.656854249492380195206754896838	/* sqrt(32) */
#endif

#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724	/* log10(2) */
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI	0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT_2dPI
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#endif


#ifndef M_LN_2PI
#define M_LN_2PI	1.837877066409345483560659472811	/* log(2*pi) */
#endif

#ifndef M_LN_SQRT_PI
#define M_LN_SQRT_PI	0.572364942924700087071713675677	/* log(sqrt(pi))
								   == log(pi)/2 */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi))
								 == log(2*pi)/2 */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2	0.225791352644727432363097614947	/* log(sqrt(pi/2))
								   == log(pi/2)/2 */
#endif


#ifdef MATHLIB_STANDALONE
# ifndef R_EXT_BOOLEAN_H_
/* "copy-paste" R_ext/Boolean.h if not already included: */
 #define R_EXT_BOOLEAN_H_
 #undef FALSE
 #undef TRUE
 typedef enum { FALSE = 0, TRUE } Rboolean;
# endif
#else
# include <R_ext/Boolean.h>
#endif

#define dnorm dnorm4
#define pnorm pnorm5
#define qnorm qnorm5

#ifdef  __cplusplus
extern "C" {
#endif
	/* R's versions with !R_FINITE checks */

double R_pow(double x, double y);
double R_pow_di(double, int);

	/* Random Number Generators */

double	norm_rand(void);
double	unif_rand(void);
double  R_unif_index(double);
double	exp_rand(void);
#ifdef MATHLIB_STANDALONE
void	set_seed(unsigned int, unsigned int);
void	get_seed(unsigned int *, unsigned int *);
#endif

	/* Normal Distribution */

double	dnorm(double, double, double, int);
double	pnorm(double, double, double, int, int);
double	qnorm(double, double, double, int, int);
double	rnorm(double, double);
void	pnorm_both(double, double *, double *, int, int);/* both tails */

	/* Uniform Distribution */

double	dunif(double, double, double, int);
double	punif(double, double, double, int, int);
double	qunif(double, double, double, int, int);
double	runif(double, double);

	/* Gamma Distribution */

double	dgamma(double, double, double, int);
double	pgamma(double, double, double, int, int);
double	qgamma(double, double, double, int, int);
double	rgamma(double, double);

double  log1pmx(double); /* Accurate log(1+x) - x, {care for small x} */
double  log1pexp(double); // <-- ../nmath/plogis.c
double  log1mexp(double);
double  lgamma1p(double);/* accurate log(gamma(x+1)), small x (0 < x < 0.5) */

double  logspace_add(double, double);
double  logspace_sub(double, double);
double  logspace_sum(const double *, int);

	/* Beta Distribution */

double	dbeta(double, double, double, int);
double	pbeta(double, double, double, int, int);
double	qbeta(double, double, double, int, int);
double	rbeta(double, double);

	/* Chi-squared Distribution */

double	dchisq(double, double, int);
double	pchisq(double, double, int, int);
double	qchisq(double, double, int, int);
double	rchisq(double);

	/* Non-central Chi-squared Distribution */

double	dnchisq(double, double, double, int);
double	pnchisq(double, double, double, int, int);
double	qnchisq(double, double, double, int, int);
double	rnchisq(double, double);

	/* F Distibution */

double	df(double, double, double, int);
double	pf(double, double, double, int, int);
double	qf(double, double, double, int, int);
double	rf(double, double);

	/* Student t Distibution */

double	dt(double, double, int);
double	pt(double, double, int, int);
double	qt(double, double, int, int);
double	rt(double);

	/* Binomial Distribution */

double  dbinom_raw(double x, double n, double p, double q, int give_log);
double	dbinom(double, double, double, int);
double	pbinom(double, double, double, int, int);
double	qbinom(double, double, double, int, int);
double	rbinom(double, double);

	/* Multinomial Distribution */

void	rmultinom(int, double*, int, int*);

	/* Exponential Distribution */

double	dexp(double, double, int);
double	pexp(double, double, int, int);
double	qexp(double, double, int, int);
double	rexp(double);

	/* Hypergeometric Distibution */

double	dhyper(double, double, double, double, int);
double	phyper(double, double, double, double, int, int);
double	qhyper(double, double, double, double, int, int);
double	rhyper(double, double, double);

	/* Poisson Distribution */

double	dpois_raw (double, double, int);
double	dpois(double, double, int);
double	ppois(double, double, int, int);
double	qpois(double, double, int, int);
double	rpois(double);

	/* Non-central Beta Distribution */

double	dnbeta(double, double, double, double, int);
double	pnbeta(double, double, double, double, int, int);
double	qnbeta(double, double, double, double, int, int);
double	rnbeta(double, double, double);

	/* Non-central Student t Distribution */

double	pnt(double, double, double, int, int);

	/* Studentized Range Distribution */

double	ptukey(double, double, double, double, int, int);
double	qtukey(double, double, double, double, int, int);

	/* Gamma and Related Functions */
double	gammafn(double);
double	lgammafn(double);
double	lgammafn_sign(double, int*);

double	beta(double, double);
double	lbeta(double, double);

double	choose(double, double);
double	lchoose(double, double);

	/* Bessel Functions */

double	bessel_i(double, double, double);
double	bessel_j(double, double);
double	bessel_k(double, double, double);
double	bessel_y(double, double);
double	bessel_i_ex(double, double, double, double *);
double	bessel_j_ex(double, double, double *);
double	bessel_k_ex(double, double, double, double *);
double	bessel_y_ex(double, double, double *);


	/* General Support Functions */

int	imax2(int, int);
int	imin2(int, int);
double	fmax2(double, double);
double	fmin2(double, double);
double	fsign(double, double);

/* More accurate cos(pi*x), sin(pi*x), tan(pi*x)

   These declarations might clash with system headers if someone had
   already included math.h with __STDC_WANT_IEC_60559_FUNCS_EXT__
   defined (and we try, above).
   We check for that via the value of __STDC_IEC_60559_FUNCS__
*/
#if !(defined(__STDC_IEC_60559_FUNCS__) && __STDC_IEC_60559_FUNCS__ >= 201506L)
double cospi(double);
double sinpi(double);
double tanpi(double);
#endif
double Rtanpi(double); /* our own in any case */

/* Compute the log of a sum or difference from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 * or  log (exp (logx) - exp (logy))
 *
 * without causing overflows or throwing away too much accuracy:
 */
double  logspace_add(double logx, double logy);
double  logspace_sub(double logx, double logy);


/* ----------------- Private part of the header file ------------------- */

#if defined(MATHLIB_STANDALONE) && !defined(MATHLIB_PRIVATE_H)
/* second is defined by nmath.h */

/* If isnan is a macro, as C99 specifies, the C++
   math header will undefine it. This happens on macOS */
# ifdef __cplusplus
  int R_isnancpp(double); /* in mlutils.c */
#  define ISNAN(x)     R_isnancpp(x)
# else
#  define ISNAN(x)     (isnan(x)!=0)
# endif

# define R_FINITE(x)    R_finite(x)
int R_finite(double);

# ifdef _WIN32  /* not Win32 as no config information */
#  ifdef RMATH_DLL
#   define R_EXTERN extern __declspec(dllimport)
#  else
#   define R_EXTERN extern
#  endif
R_EXTERN double NA_REAL;
R_EXTERN double R_PosInf;
R_EXTERN double R_NegInf;
R_EXTERN int N01_kind;
#  undef R_EXTERN
#else
extern int N01_kind;
# endif

#endif /* MATHLIB_STANDALONE */

#ifdef  __cplusplus
}
#endif

#endif /* RMATH_H */
