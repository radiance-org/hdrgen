/*
 *  pq2.h	- convert to/from PQ
 *
 *  Created by Greg Ward on 7/28/15.
 *  Copyright 2014 Dolby Laboratories. All rights reserved.
 *
 */

#pragma once

#define PQpeak	10000.
#define PQn	(2610./(4096.*4.))
#define PQm	(2523./4096.*128.)
#define PQc1	(3424./4096.)
#define PQc2	(2413./4096.*32.)
#define PQc3	(2392./4096.*32.)

#undef max

inline double
max(double a, double b)
{
	if (a > b) return a;
	return b;
}

inline float
toPQ(float Y)
{
	if (Y <= 0) return .0;
	if (Y >= PQpeak) return 1.;

	double	Yn = pow(Y*(1./PQpeak), PQn);
	
	return float(pow((PQc1 + PQc2*Yn)/(1. + PQc3*Yn), PQm));
}

inline float
fromPQ(float q)
{
	if (q <= 0) return .0;
	if (q >= 1.) return PQpeak;
	
	double	qm = pow(q, 1./PQm);
	
	return float(pow(max(qm - PQc1, 0.)/(PQc2 - PQc3*qm), 1./PQn) * PQpeak);
}
