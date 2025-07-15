/*
 * Compute exposure that optimizes HDR image range within specified limits
 *
 *	Jan. 6, 2022	G. Ward
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dmessage.h"
#include "pimage.h"

#define HSIZE	100
unsigned long	histo[HSIZE];		/* our histogram */

int
main(int argc, char *argv[])
{
	float		extrema[2];
	double		minv, maxv;
	int		hlower, hupper, tarlen;
	unsigned long	lost[2];
	ImgReader *	ir;
	ImgStruct	ims;
					/* exit on serious errors */
	dmessage_class_flags[DMCsystem] |= DMFexit;
	dmessage_class_flags[DMCparameter] |= DMFexit;
	dmessage_class_flags[DMCresource] |= DMFexit;
	dmessage_class_flags[DMCdata] |= DMFexit;
	dmessage_class_flags[DMCinput] |= DMFexit;
					/* check input */
	if (argc != 4 || (minv = atof(argv[1])) <= 0 ||
			(maxv = atof(argv[2])) <= minv) {
		DMESGF(DMCinput, "%s vmin vmax input.hdr", argv[0]);
		return 1;
	}
	PloadStandardReaders();		/* open HDR image */
	ir = PopenImageF(argv[3], 0, NULL);
	if (!ir) return 1;
	if (ir->cs.dtype != IDTfloat) {
		DMESG(DMCdata, "Input file must be HDR (float)");
		return 1;
	}
	ims.xres = ir->xres;		/* read/convert to grayscale */
	ims.yres = ir->yres;
	ims.csp = &ICS_Y;
	ims.img = NULL;
	if (!PrenderImageR(&ims, ir, 0))
		return 1;
	IRclose(ir);
	extrema[0] = 1; extrema[1] = 0;	/* get extrema */
	if (!PcomputeHisto(extrema, NULL, 0, &ims, PHCluminance))
		return 1;
	if (minv <= 0) {
		DMESG(DMCdata, "Some input gray levels are <= 0");
		return 1;
	}
	if ((minv <= extrema[0]) & (extrema[1] <= maxv)) {
		puts("1.");		/* we're in range */
		return 0;
	}
	if (maxv/minv >= extrema[1]/extrema[0]) {
					/* we have enough range */
		printf("%.1e\n", sqrt(maxv*minv/(extrema[1]*extrema[0])));
		return 0;
	}
					/* else need to optimize */
	if (!PcomputeHisto(extrema, histo, HSIZE, &ims, PHCluminance)) {
		DMESG(DMCdata, "Black image?");
		return 1;
	}
					/* target dyn range */
	tarlen = HSIZE * log(maxv/minv) / log(extrema[1]/extrema[0]);
	hlower = 0; hupper = HSIZE;	/* trim each end */
	lost[0] = lost[1] = 0;		/* balancing losses */
	while (hupper - hlower > tarlen)
		if (lost[0]+histo[hlower] <= lost[1]+histo[hupper-1])
			lost[0] += histo[hlower++];
		else
			lost[1] += histo[--hupper];

	printf("%.2e\n", maxv/( extrema[0] * pow(extrema[1]/extrema[0],
							hupper*(1./HSIZE)) ));
	return 0;
}
