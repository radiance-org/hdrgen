/*
 *  PQconvert.cpp
 *  
 *  Convert to/from 16-bit PQ encoding (forces P3 gamut on PQ side)
 *
 *  Created by Greg Ward on 11/24/14.
 *  Copyright 2014 Dolby Laboratories. All rights reserved.
 *
 */

#if defined(_WIN32) || defined(_WIN64)
#define	USE_FORK 0
#else
#define USE_FORK 1
#endif

#include "system.h"
#include "random.h"
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#if USE_FORK
#include <sys/wait.h>
#endif
#include "dmessage.h"
#include "panimage.h"
#include "pq2.h"

#ifndef EXACT2FLOAT
#define EXACT2FLOAT	0		// exact PQ-to-float conversion?
#endif

static ImgColorSpace	ICS_16bit_lin = {
		IDTushort, IPFrgb, true, 1, 0,
		{{.680, .320}, {.265, .690}, {.150, .060}, {.3140, .3510}}
	};

static ImgColorSpace	ICS_real = {
		IDTfloat, IPFrgb, true, 1, 0,
		{{.680, .320}, {.265, .690}, {.150, .060}, {.3140, .3510}}
	};


static float	PQtable12[4097];	// 12-bit PQ lookup table

static inline int
rand4(void)				// random value 0-15
{
	return random() >> 27;
}

// Convert from 16-bit PQ color space to floating-point output
static bool
convertFromPQ(const char *outName, PanImage *inp, float calibf)
{
	DMESGF(DMCinfo, "Converting PQ input to HDR '%s'...", outName);

	if (!inp->AssertCS(ICS_16bit_lin)) {	// enforce this!
		DMESG(DMCdata, "Input image does not have 3-component 16-bit");
		return false;
	}
#if EXACT2FLOAT
	if (!inp->ConvertCS(ICS_real)) {
		DMESG(DMCdata, "Cannot convert color space to float");
		return false;
	}
	for (int y = 0; y < inp->Height(); y++) {
		float *	pp = inp->RowFloat(y);
		for (int n = inp->NComp()*inp->Width(); n--; pp++)
			*pp = fromPQ(*pp);
	}
#else
	int	n = 4097*(PQtable12[4096] < 1.f);
	while (n-- > 0)				// initialize table first run
		PQtable12[n] = fromPQ(n*(1.f/4096.f));
	PanImage	res(inp->Width(), inp->Height(), ICS_real, PICready);
	for (int y = 0; y < inp->Height(); y++) {
		const unsigned short *	sp = inp->GetRowShort(y);
		float *			pp = res.RowFloat(y);
		for (n = inp->NComp()*inp->Width(); n--; )
			*pp++ = PQtable12[(*sp++ + rand4()) >> 4];
	}
	if (!inp->Take(&res))
		return false;
#endif
	if (calibf > 1e-5) {
		inp->Info(IIFstonits)->stonits = calibf;
	} else if (calibf < -1e-5) {
		*inp /= -calibf;
		inp->Info(IIFstonits)->stonits = 1.f;
	}
	return (inp->Write(outName) > 0);
}

// Convert from floating-point to PQ color space with P3 gamut
static bool
convertToPQ(const char *outName, PanImage *inp, int nbits = 16)
{
	float	s2nits = 1.f;
	int	y, n;

	sprintf(dmessage_buf, "Converting HDR input to %d-bit/chan PQ '%s'...",
			nbits, outName);
	DMESG(DMCinfo, dmessage_buf);

	if (inp->GetInfo(IIFstonits))
		s2nits = inp->GetInfo()->stonits;
	if (!inp->ConvertCS(ICS_real, s2nits))
		return false;
	for (y = 0; y < inp->Height(); y++) {
		float *	pp = inp->RowFloat(y);
		DASSERT(pp != NULL);
		for (n = inp->NComp()*inp->Width(); n--; pp++)
			*pp = toPQ(*pp);
	}
	if (!inp->ConvertCS(ICS_16bit_lin)) {
		DMESG(DMCdata, "Cannot convert color space to 16-bit");
		return false;
	}
	if (nbits < 16) {			// reduce bit resolution
		const int	mask = ~0<<(16-nbits);
		for (y = 0; y < inp->Height(); y++) {
			unsigned short *	pp = inp->RowShort(y);
			for (n = inp->NComp()*inp->Width(); n--; )
				*pp++ &= mask;
		}
	}
	return (inp->Write(outName) > 0);
}

// Check for %d in output format string
static bool
hasFrameFormat(const char *fmt)
{
	const char *	cp;
	while ((cp = strchr(fmt, '%')) != NULL) {
		if (cp[1] == '%') {
			fmt = cp+2;
			continue;
		}
		while (isdigit(*++cp))
			;
		return (*cp && strchr("diouXx", *cp) != NULL);
	}
	return false;
}

int
main(int argc, char *argv[])
{
	int		nproc = 1;
	const char *	av0 = argv[0];
	int		nbits = 16;
	float		calibf = 0;
	int		startf = 1;
	int		xsiz = 0, ysiz = 0;
						// turn on for debug output
	// dmessage_class_flags[DMCtrace] |= DMFstderr;
	dmessage_class_flags[DMCwarning] |= DMFstderr;
						// check arguments
	while (argc > 3 && argv[1][0] == '-') {
		switch (argv[1][1]) {
		case 's':			// start frame
			startf = atoi(argv[2]);
			argv += 2;
			argc -= 2;
			continue;
#if USE_FORK
		case 'n':			// number of processes
			nproc = atoi(argv[2]);
			if (nproc <= 0) {
				DMESG(DMCinput, "Illegal # processes");
				break;
			}
			argv += 2;
			argc -= 2;
			continue;
#endif
		case 'c':			// calibration
			calibf = atof(argv[2]);
			argv += 2;
			argc -= 2;
			continue;
		case 'o':			// color space
			if (!strcasecmp(argv[2], "P3"))
				memcpy(ICS_real.chroma, ICS_P3.chroma, sizeof(ICS_real.chroma));
			else if (!strcasecmp(argv[2], "RGB2020"))
				memcpy(ICS_real.chroma, ICS_RGB2020.chroma, sizeof(ICS_real.chroma));
			else if (!strcasecmp(argv[2], "RGB709"))
				memcpy(ICS_real.chroma, ICS_RGB709.chroma, sizeof(ICS_real.chroma));
			else {
				DMESGF(DMCinput, "Unsupported color space '%s'", argv[2]);
				break;
			}
			memcpy(ICS_16bit_lin.chroma, ICS_real.chroma, sizeof(ICS_real.chroma));
			argv += 2;
			argc -= 2;
			continue;
		case 'a':			// adjustment factor
			calibf = -atof(argv[2]);
			argv += 2;
			argc -= 2;
			continue;
		case 'd':			// dimensions
			xsiz = atoi(argv[2]);
			ysiz = atoi(argv[3]);
			argv += 3;
			argc -= 3;
			continue;
		case 'v':			// verbose mode
			dmessage_class_flags[DMCinfo] |= DMFstderr;
			++argv; --argc;
			continue;
		case 'b':			// PQ bit resolution
			nbits = atoi(argv[2]);
			argv += 2;
			argc -= 2;
			continue;
		default:
			break;
		}
		break;
	}
	if (argc != 3) {
		fputs("Usage: ", stderr);
		fputs(av0, stderr);
		fputs(" [-v][-n nproc][-d xsiz ysiz][-b #bits][-o {P3|RGB2020|RGB709}][-c stonits | -a mult][-s start#] input_spec output_spec\n", stderr);
		return 1;
	}
	const bool	looping = hasFrameFormat(argv[1]);
	if (looping ^ hasFrameFormat(argv[2])) {
		DMESG(DMCinput, "Only one file specification has format!");
		return 1;
	}
	if (!looping) nproc = 1;
						// load i/o we may need
	PloadStandardReaders();
	PaddIReaderI(&IRInterfaceMTX);
	PloadStandardWriters();
	PanImage::AddImageWriter(&IWInterfaceMTX);
	int	myStart = startf;
	int	i = nproc;
#if USE_FORK
	while (--i > 0 && fork() > 0)		// fork additional processes
		++myStart;
#endif
	for (i = myStart; (i == startf) | looping; i += nproc) {
		ImgColorSpace	mySpace;
		char		fname[512];
		sprintf(fname, argv[1], i);
		ImgReader *	ir = PopenImageF(fname, (i > startf), NULL);
		if (ir == NULL)
			break;			// normal loop exit
		if (!xsiz) xsiz = ir->xres;
		if (!ysiz) ysiz = ir->yres;
		PanImage	inpImg(xsiz, ysiz, ir->cs);
		if (!inpImg.Load(ir))
			return 1;
		IRclose(ir);
		sprintf(fname, argv[2], i);
		switch (inpImg.GetCS()->dtype) {
		case IDTfloat:
			if (calibf > 1e-7)
				inpImg.Info(IIFstonits)->stonits = calibf;
			else if (calibf < -1e-7)
				inpImg.Info(IIFstonits)->stonits *= -calibf;
			if (!convertToPQ(fname, &inpImg, nbits))
				return 1;
			break;
		case IDTubyte:
			PcopyCS(&mySpace, inpImg.GetCS());
			mySpace.dtype = IDTushort;
			if (!inpImg.ConvertCS(mySpace))
				return 1;
			// fall through
		case IDTushort:
			if (!convertFromPQ(fname, &inpImg, calibf))
				return 1;
			break;
		default:
			DMESG(DMCdata, "Unsupported input data type");
			return 1;
		}
	}
	int	rstatus = (i == startf);	// parent waits for children
#if USE_FORK
	if (nproc > 1 && myStart == startf+nproc-1) {
		int	status;
		while (--nproc && wait(&status) > 0)
			if (status)
				rstatus = 1;
	}
#endif
	return rstatus;
}
