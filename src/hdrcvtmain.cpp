/*
 *  hdrcvtmain.cpp
 *  hdrcvt
 *
 *  Convert between HDR & LDR formats
 *
 *  Created by Greg Ward on 7 Feb 2009.
 *  Copyright (c) 2009 Anyhere Software. All rights reserved.
 *
 */

#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "rtio.h"
#include "tmprivat.h"
#include "system.h"
#include "dmessage.h"
#include "panimage.h"
#include "radheader.h"

const char	USE_MSG[] = "Usage: %s [-verbose][-dim X Y | -scale F]\
[-quality Q][-expose E][-calibrate M | -limit %%][-flare]\
[-wb F][-rotate {0|90|180|270}[hv]][-blur b][-dilate r][-comment C][-p param=val[:p2=v2..]]\
[-match match.img | -tmo T][-input ICS][-output OCS][-swap] input.img output.img\n";

// Compute and remove lens flare from a high dynamic range image
extern bool	PHDremoveFlare(ImgStruct *ims,
				int (*pbar)(const char *, int) = NULL);

#define MIN_GAMMA	1.6f		// minimum expected gamma
#define DEF_GAMMA	(ICS_sRGB.gamma)

// Struct to hold our tone-mapping callbacks
struct TMOperator {
	const char	name[16];	// short name for user
	const char *	descr;		// TMO description
					// pointer to main operator
	bool		(*fptr)(PanImage *res, PanImage *inp, int ndx);
					// global map computation
	bool		(*mptr)(TMstruct *tms, int ndx);
};

extern TMOperator	TMOper[];	// forward declaration

#define ApplyTMO(dptr, sptr, n)	(*TMOper[n].fptr)(dptr, sptr, n)
#define ComputeMap(tms, n)	(*TMOper[n].mptr)(tms, n)

bool		do_human = false;	// model HVS?

// Tone-mapping object class
struct TMapObject {
	TMstruct *	tms;		// tone-mapping structure
	ImgStruct	lumImg;		// 16-bit log luminance image
	ImgStruct	rgbImg;		// 24-bit/pixel gamma'ed chroma image
	ImgInfo		inpInfo;	// input information 
			TMapObject(const PanImage &pim,
					const ImgColorSpace *dcsp = 0) {
				tms = 0; lumImg.img = 0; rgbImg.img = 0;
				Init(pim, dcsp);
			}
			~TMapObject() {
				Init();
			}
	void		Init() {
				tmDone(tms); tms = 0;
				PfreeImage(&lumImg);
				PfreeImage(&rgbImg);
			}
	bool		Init(const PanImage &pim, const ImgColorSpace *dcsp = 0);
	bool		Ready() const {
				return (tms != NULL) & (lumImg.img != NULL);
			}
	bool		MapGlobal(PanImage *pip) const;
};

// Initialize tone-mapping and create log image
bool
TMapObject::Init(const PanImage &pim, const ImgColorSpace *dcsp)
{
	Init();
	if (pim.Ready() < PICread)
		return false;
	if (pim.GetCS()->format != IPFrgb && pim.GetCS()->format != IPFxyz &&
			pim.GetCS()->format != IPFy) {
		DMESG(DMCparameter, "Tone-mapping input must be RGB, Y or XYZ");
		return false;
	}
	DASSERT(pim.GetCS()->logorig == 0);
	pim.GetInfo(&inpInfo);		// copy information header
	double	stonits = 1.;
	if (inpInfo.flags & IIFstonits)
		stonits = inpInfo.stonits;
	else
		do_human = false;
					// initialize TMO
	int	myFlags = do_human ? TM_F_HUMAN : TM_F_CAMERA;
	if (pim.GetCS()->format == IPFy ||
			(dcsp != NULL && dcsp->format == IPFy))
		myFlags |= TM_F_BW;
	if (dcsp != NULL)
		tms = tmInit(myFlags, (RGBPRIMP)dcsp->chroma, dcsp->gamma);
	else if (pim.GetCS()->format == IPFrgb)
		tms = tmInit(myFlags, (RGBPRIMP)pim.GetCS()->chroma, .0);
	else
		tms = tmInit(myFlags, (RGBPRIMP)ICS_sRGB.chroma, ICS_sRGB.gamma);
	if (tms == NULL) {
		DMESG(DMCdata, "Cannot initialize tone-mapping");
		return false;
	}
	if (tmSetSpace(tms, (pim.GetCS()->format == IPFrgb) ?
			(RGBPRIMP)pim.GetCS()->chroma : TM_XYZPRIM,
				stonits) != TM_E_OK) {
		Init();
		return false;
	}
	ImgColorSpace	logCS;
	PcopyCS(&logCS, &ICS_Y16);	// this is a bit of a lie
	logCS.gamma = 1.f;
	lumImg.xres = pim.Width(); lumImg.yres = pim.Height();
	lumImg.csp = &logCS;
	if (!PnewImage(&lumImg, .0)) {
		Init();
		return false;
	}
	if (!(tms->flags & TM_F_BW)) {
		rgbImg.xres = pim.Width(); rgbImg.yres = pim.Height();
		rgbImg.csp = &ICS_sRGB;
		if (!PnewImage(&rgbImg, .0)) {
			Init();
			return false;
		}
	}
	int	y;			// convert input from...
	switch ((pim.GetCS()->format == IPFy) << 1 |
			(pim.GetCS()->dtype == IDTushort)) {
	case 0:				// float RGB or XYZ
		for (y = 0; y < pim.Height(); y++)
			if (tmCvColors(tms, (TMbright *)ProwPtr(&lumImg,y),
					(tms->flags & TM_F_BW) ? TM_NOCHROM :
							ProwPtr(&rgbImg,y),
			(COLORV (*)[3])const_cast<float *>(pim.GetRowFloat(y)),
					pim.Width()) != TM_E_OK) {
				Init();
				return false;
			}
		break;
	case 1:				// short RGB or XYZ
		for (y = 0; y < pim.Height(); y++)
			if (tmCvRGB48(tms, (TMbright *)ProwPtr(&lumImg,y),
					(tms->flags & TM_F_BW) ? TM_NOCHROM :
							ProwPtr(&rgbImg,y),
		(uint16 (*)[3])const_cast<unsigned short *>(pim.GetRowShort(y)),
					pim.Width(),
					pim.GetCS()->gamma) != TM_E_OK) {
				Init();
				return false;
			}
		break;
	case 2:				// float Y
		for (y = 0; y < pim.Height(); y++)
			if (tmCvGrays(tms, (TMbright *)ProwPtr(&lumImg,y),
					const_cast<float *>(pim.GetRowFloat(y)),
					pim.Width()) != TM_E_OK) {
				Init();
				return false;
			}
		break;
	case 3:				// short Y
		for (y = 0; y < pim.Height(); y++)
			if (tmCvGray16(tms, (TMbright *)ProwPtr(&lumImg,y),
				const_cast<unsigned short *>(pim.GetRowShort(y)),
					pim.Width(),
					pim.GetCS()->gamma) != TM_E_OK) {
				Init();
				return false;
			}
		break;
	}
					// compute histogram
	tmAddHisto(tms, (TMbright *)lumImg.img, lumImg.xres*lumImg.yres, 1);
	return true;
}

// Put global tone-mapping results into the given image
bool
TMapObject::MapGlobal(PanImage *pip) const
{
	if (!Ready())
		return false;
	if (!pip || !pip->Ready())
		return false;
	DASSERT((pip->Width() == lumImg.xres) & (pip->Height() == lumImg.yres));
	DASSERT(pip->GetCS()->dtype == IDTubyte);
	if (!pip->Init(PICready))	// allocate destination
		return false;
	*pip->Info() = inpInfo;		// transfer source information
	pip->Info()->flags &= ~IIFstonits;
	for (int y = 0; y < lumImg.yres; y++)
		if (tms->flags & TM_F_BW) {
			uby8 *	dptr = pip->RowByte(y);
			if (tmMapPixels(tms, dptr,
					(TMbright *)ProwPtr(&lumImg,y),
					TM_NOCHROM, lumImg.xres) != TM_E_OK)
				return false;
			if (pip->GetCS()->format == IPFrgb)
				PgetRGB24fromY8(dptr, dptr, pip->Width());
		} else if (tmMapPixels(tms, pip->RowByte(y),
				(TMbright *)ProwPtr(&lumImg,y),
				ProwPtr(&rgbImg,y), lumImg.xres) != TM_E_OK)
			return false;
	return true;
}

// Add TMO descriptor to info history
static void
tagTMO(ImgInfo *infp, const char *mapDescr)
{
	const char	*preamble = "HDR mapped to 8-bit/channel using ";
	char *		cp = infp->comments;
	if (infp->flags & IIFcomments)
		while (*cp) cp++;
	while (*preamble)
		*cp++ = *preamble++;
	while (*mapDescr)
		*cp++ = *mapDescr++;
	*cp++ = '\n'; *cp = '\0';
	infp->flags |= IIFcomments;
}

// Histogram adjustment global tone-mapping operator
static bool
mapHisto(TMstruct *tms, int ndx)
{
	return (tmComputeMapping(tms, .0, .0, .0) == TM_E_OK);
}

// Look up and linearly interpolate TM value
static inline double
tmlerp(const double wlum, const float (*map)[2], const int len)
{
	static int	i = 0;
	double		x;
					// update index from last call
	if (i >= len) i = len-1;
	while (i >= 0 && wlum < map[i][0])
		--i;
	if (i < 0)
		return double(map[i=0][1]);

	while (i < len-1 && wlum >= map[i+1][0])
		++i;
	if (i >= len-1)
		return double(map[len-1][1]);

	x = (wlum - map[i][0])/(map[i+1][0] - map[i][0]);

	return (1.-x)*map[i][1] + x*map[i+1][1];
}

// Load pre-defined global curve from stdin
static bool
mapLoad(TMstruct *tms, int ndx)
{
#define MAXMAPLEN	500
	float		map[MAXMAPLEN][2];
	int		maplen = 0;
					// load map from stdin
	while (scanf("%f %f\n", &map[maplen][0], &map[maplen][1]) == 2) {
		if ((map[maplen][0] < 0) | (map[maplen][1] < 0)) {
			DMESGF(DMCdata, "Negative TM value at line %d\n",
					maplen+1);
			return false;
		}
		if (maplen && (map[maplen][0] <= map[maplen-1][0]) |
				(map[maplen][1] < map[maplen-1][1])) {
			DMESGF(DMCdata, "Input TM not monotonic at line %d",
					maplen+1);
			return false;
		}
		if (++maplen >= MAXMAPLEN && fgetc(stdin) != EOF) {
			DMESGF(DMCdata, "Input TM has more than %d entries\n",
					MAXMAPLEN);
			return false;
		}
	}
	if (!feof(stdin))
		DMESG(DMCwarning, "Extra characters at end of tone curve");
					// assign tone-mapping
	if (!tmNewMap(tms)) {
		DMESG(DMCmemory, "Cannot allocate tone map");
		return false;
	}
	for (int i = tms->mbrmax-tms->mbrmin+1; i--; ) {
		double	val =  tmlerp(tmLuminance(i+tms->mbrmin), map, maplen);
		tms->lumap[i] = TM_BRES*pow((val - map[0][1]) /
				(map[maplen-1][1] - map[0][1]), 1./tms->mongam);
	}
	return true;
#undef MAXMAPLEN
}

// Compute global gamma tone curve
static bool
mapGamma(TMstruct *tms, int ndx)
{
	const int	maxV = (1L<<(8*sizeof(TMAP_TYP))) - 1;
	const int	horig = HISTI(tms->hbrmin);
	const int	hlen = HISTI(tms->hbrmax) + 1 - horig;
	HIST_TYP	htot = 0;
	long		cnt;
	int		i;
	TMbright	lmin, lmax;
					// get histogram total
	for (i = hlen; i-- > 0; )
		htot += tms->histo[i];
	if (!htot) {
		DMESG(DMCdata, "Cannot tone-map black image");
		return false;
	}
					// ignore bottom 2%
	cnt = htot/50;
	for (i = 0; (i < hlen-1) & (cnt > 0); i++)
		cnt -= tms->histo[i];
	lmin = HISTV(horig+i);
					// ignore top 0.5%
	cnt = htot/200;
	for (i = hlen; (--i > 0) & (cnt > 0); )
		cnt -= tms->histo[i];
	lmax = HISTV(horig+i);
					// assign tone-mapping
	if (!tmNewMap(tms)) {
		DMESG(DMCmemory, "Cannot allocate tone map");
		return false;
	}
	const double	mult = log(DEFLDDYN) / (tms->mongam*(lmax-lmin));
	const double	offs = mult * (lmax - lmin);
	for (i = tms->mbrmax-tms->mbrmin+1; i--; ) {
		double	val = TM_BRES * exp((tms->mbrmin+i-lmin)*mult - offs);
		tms->lumap[i] = (val >= double(maxV)) ? maxV : (int)val;
	}
	return true;
}

// Reinhard global photographic tone-mapping operator
static bool
mapPhoto(TMstruct *tms, int ndx)
{
	const int	horig = HISTI(tms->hbrmin);
	const int	hlen = HISTI(tms->hbrmax) + 1 - horig;
	HIST_TYP	htot = 0;
	long		cnt;
	double		lavg = 0;
	int		i;
	TMbright	lmin, lmax;
					// get histogram total & mean
	for (i = hlen; i-- > 0; ) {
		htot += tms->histo[i];
		lavg += (double)HISTV(horig+i)*tms->histo[i];
	}
	if (!htot) {
		DMESG(DMCdata, "Cannot tone-map black image");
		return false;
	}
	lavg /= (double)htot;
					// ignore bottom 2%
	cnt = htot/50;
	for (i = 0; (i < hlen-1) & (cnt > 0); i++)
		cnt -= tms->histo[i];
	lmin = HISTV(horig+i);
					// ignore top 0.5%
	cnt = htot/200;
	for (i = hlen; (--i > 0) & (cnt > 0); )
		cnt -= tms->histo[i];
	lmax = HISTV(horig+i);
					// assign tone-mapping
	if (!tmNewMap(tms)) {
		DMESG(DMCmemory, "Cannot allocate tone map");
		return false;
	}
					/* compute Reinhard parameters */
	const double	alph = .18*pow(4., (2.*lavg-lmin-lmax)/(lmax-lmin));
	const double	gsca = alph/tmLuminance(lavg);
	for (i = tms->mbrmax-tms->mbrmin+1; i--; ) {
		double	val = gsca*tmLuminance(tms->mbrmin+i);
		val = TM_BRES*pow(val/(1. + val), 1./tms->mongam);
		// tms->lumap[i] = (val >= double(0xffff)) ? 0xffff : (int)val;
		tms->lumap[i] = (int)val;
	}
	return true;
}

// Generic global tone-mapping operation
static bool
TMOglobal(PanImage *res, PanImage *inp, int ndx)
{
	DMESGF(DMCtrace, "Applying %s", TMOper[ndx].descr);

	TMapObject	myMap(*inp, res->GetCS());

	if (!myMap.Ready())
		return false;

	if (!inp->ReadOnly())		// sign that it's OK to free source
		inp->Init(PICfree);

	if (!ComputeMap(myMap.tms, ndx))
		return false;

	if (!myMap.MapGlobal(res))
		return false;
	
	tagTMO(res->Info(), TMOper[ndx].descr);
	return true;
}

// Bilateral filter local tone-mapping operator
static bool
TMObilat(PanImage *res, PanImage *inp, int ndx)
{
	if (!res || !res->Ready())
		return false;

	DASSERT((res->Width() == inp->Width()) & (res->Height() == inp->Height()));
	DASSERT(res->GetCS()->dtype == IDTubyte);

	DMESGF(DMCtrace, "Applying %s", TMOper[ndx].descr);

	TMapObject	myMap(*inp, res->GetCS());

	if (!myMap.Ready())
		return false;

	if (!inp->ReadOnly())		// sign that it's OK to free source
		inp->Init(PICfree);
					// compute our base luminance mapping
	if (!ComputeMap(myMap.tms, ndx))
		return false;
	ImgColorSpace	logSpace;	// bring LogL into positive domain
	int		x, y;
	for (y = 0; y < myMap.lumImg.yres; y++) {
		TMbright *	lp = (TMbright *)ProwPtr(&myMap.lumImg, y);
		for (x = myMap.lumImg.xres; x--; )
			if (*lp <= MINBRT)
				*lp++ = 0;
			else
				*lp++ -= MINBRT;
	}
	PcopyCS(&logSpace, &ICS_Y16);
	logSpace.logorig = MINLUM;
	logSpace.gamma = double(1<<16)/TM_BRTSCALE/M_LN10;
	myMap.lumImg.csp = &logSpace;
					// compute bilateral filter on LogL
	PanImage	BLFimg(myMap.lumImg.xres, myMap.lumImg.yres, logSpace);
	if (!PbilatFilter(BLFimg.Img(), &myMap.lumImg,
			sqrt((double)myMap.lumImg.xres*myMap.lumImg.yres)/64.,
			log(1.8)*TM_BRTSCALE/(1<<16), NULL))
		return false;

	if (!res->Init(PICready))	// allocate destination
		return false;

	*res->Info() = myMap.inpInfo;	// transfer information
	res->Info()->flags &= ~IIFstonits;
					// apply bilateral filter TMO
	const int	plen = res->NComp();
	const double	gamcor = 1./TM_BRTSCALE/myMap.tms->mongam;
	if (!(myMap.tms->flags & TM_F_BW))
		DASSERT(plen == 3);
	for (y = 0; y < res->Height(); y++) {
		const unsigned short *	ilp = (unsigned short *)
						ProwPtr(&myMap.lumImg, y);
		const unsigned short *	blp = BLFimg.GetRowShort(y);
		const uby8 *		cs = (myMap.tms->flags & TM_F_BW) ?
						TM_NOCHROM :
						ProwPtr(&myMap.rgbImg, y);
		uby8 *			dp = res->RowByte(y);
		for (x = res->Width(); x--; ilp++, blp++) {
			int	li, j;
			if ((j = *blp + MINBRT) < myMap.tms->mbrmin) {
				li = 0;
			} else {
				if (j > myMap.tms->mbrmax)
					j = myMap.tms->mbrmax;
				li = myMap.tms->lumap[j - myMap.tms->mbrmin];
					// restore detail layer
				li = int( li*exp((*ilp - *blp)*gamcor) + .5 );
			}
			if (cs == TM_NOCHROM) {
				if (li > 255) li = 255;
				for (j = plen; j--; )
					*dp++ = li;
			} else {
				j = *cs++ * li / myMap.tms->cdiv[RED];
				*dp++ = (j>255) ? 255 : j;
				j = *cs++ * li / myMap.tms->cdiv[GRN];
				*dp++ = (j>255) ? 255 : j;
				j = *cs++ * li / myMap.tms->cdiv[BLU];
				*dp++ = (j>255) ? 255 : j;
			}
		}
	}
	tagTMO(res->Info(), TMOper[ndx].descr);
	return true;
}

#define	LVLRATIO	2		// ratio between pyramid levels
#define MINLSIZE	LVLRATIO	// minimum width/height in apex
#define	LTLEN		(TM_BRES*2)	// length of gamma lookup table

// Find scale factor for maximum value in image
static float
scaleMax(const PanImage &im)
{
	const float	pctl = 99.5f;	// ignore top 0.5%
	unsigned short	smax;
	float		fmax;

	if (!PcomputePercentiles((im.GetCS()->dtype==IDTfloat) ?
			(void *)&fmax : (void *)&smax, &pctl, 1,
			im.GetImg(), PHCluminance))
		return 1.f;
	if (im.GetCS()->dtype==IDTfloat)
		return 1.f/fmax;
	return 1.f/(im.GetCS()->logorig*exp(M_LN10/(1<<16)*im.GetCS()->gamma*smax));
}

// Log image multiplication (no checks)
static void
logMultIm(PanImage *dstp, const PanImage &srci, const double gm = 1.)
{
	const double	offset = (1<<16) * log10(dstp->GetCS()->logorig) /
						dstp->GetCS()->gamma;
	int		y, i, v;

	for (y = 0; y < dstp->Height(); y++) {
		const unsigned short *	rp = srci.GetRowShort(y);
		unsigned short *	cp = dstp->RowShort(y);
		for (i = dstp->Width(); i--; ) {
			v = int(gm*(double(*cp) + offset) + double(*rp++) + .5);
			*cp++ = v*(v > 0);
		}
	}
}

// Compute gamma for a particular level in pyramid
// Gammas summed over all levels should equal 1
static double
msGamma(int i, int N)
{
	return i * 2 / (double)(N*(N+1));
}

// Recursive call to compute multiscale local TMO from ratio image
static bool
compMultiScale(PanImage *cur, int i, const int N)
{
	if (i <= 1) {
		const double	ofs = (1<<16) * log10(cur->GetCS()->logorig) /
						cur->GetCS()->gamma;
		const double	gm = msGamma(i, N);
		for (int y = 0; y < cur->Height(); y++) {
			unsigned short *	cp = cur->RowShort(y);
			for (int x = cur->Width(); x--; cp++)
				*cp = int(gm*(double(*cp) + ofs) - ofs + .5);
		}
		return true;
	}
	PanImage	nextI(cur->Width()/LVLRATIO, cur->Height()/LVLRATIO,
					*cur->GetCS());
					// downsample image
	if (!nextI.Load(*cur))
		return false;
					// compute next level(s)
	if (!compMultiScale(&nextI, i-1, N))
		return false;
					// upsample result
	PanImage	upsizeI(cur->Width(), cur->Height(), *cur->GetCS());
	if (!PsizeImage(upsizeI.Img(), nextI.GetImg(), PSlinear))
		return false;
	nextI.Init();
					// multiply into this level
	logMultIm(cur, upsizeI, msGamma(i,N));

	return true;
}

// Multiscale local tone-mapping operator
static bool
TMOmultis(PanImage *res, PanImage *inp, int ndx)
{
	if (!res || !res->Ready())
		return false;

	DASSERT((res->Width() == inp->Width()) & (res->Height() == inp->Height()));
	DASSERT(res->GetCS()->dtype == IDTubyte);

	DMESGF(DMCtrace, "Applying %s", TMOper[ndx].descr);

	TMapObject	myMap(*inp, res->GetCS());

	if (!myMap.Ready())
		return false;

	if (!inp->ReadOnly())		// sign that it's OK to free source
		inp->Init(PICfree);
					// compute our base luminance mapping
	if (!ComputeMap(myMap.tms, ndx))
		return false;
	ImgColorSpace	logSpace;	// bring LogL into positive domain
	int		x, y;
	for (y = 0; y < myMap.lumImg.yres; y++) {
		TMbright *	lp = (TMbright *)ProwPtr(&myMap.lumImg, y);
		for (x = myMap.lumImg.xres; x--; )
			if (*lp <= MINBRT)
				*lp++ = 0;
			else
				*lp++ -= MINBRT;
	}
	PcopyCS(&logSpace, &ICS_Y16);
	logSpace.logorig = MINLUM;
	logSpace.gamma = float(1<<16)/TM_BRTSCALE/M_LN10;
	myMap.lumImg.csp = &logSpace;
					// create log(ratio) base image
	TMAP_TYP	lumTab[LTLEN];
	int		i = LTLEN;
	while (i--) {
		int	li = TM_BRTSCALE*myMap.tms->mongam*log((i+.5)*(1./TM_BRES))
				+ (.5-MINBRT);
		lumTab[i] = li * (li > 0);
	}
	int		nlevels;
	PanImage	curImg, baseRatio;
	if (!curImg.Take(&myMap.lumImg))
		return false;
	x = curImg.Width(); y = curImg.Height();
	for (nlevels = 0; (x > MINLSIZE) & (y > MINLSIZE); nlevels++)
		{ x /= LVLRATIO; y /= LVLRATIO; }
	if (nlevels <= 0) {
		DMESGF(DMCparameter, "Image resolution too low for %s",
				TMOper[ndx].descr);
		return false;
	}
	if (!baseRatio.Init(curImg.Width(), curImg.Height(), logSpace))
		return false;
	for (y = 0; y < curImg.Height(); y++) {
		const unsigned short *	sp = curImg.GetRowShort(y);
		unsigned short *	dp = baseRatio.RowShort(y);
		for (x = curImg.Width(); x--; sp++) {
			i = *sp + MINBRT;
			if (i < myMap.tms->mbrmin)
				i = myMap.tms->mbrmin;
			else if (i > myMap.tms->mbrmax)
				i = myMap.tms->mbrmax;
			i = myMap.tms->lumap[i - myMap.tms->mbrmin];
			if (i >= LTLEN)
				i = LTLEN-1;
			i = lumTab[i] - int(*sp);
			if (i <= MINBRT)
				*dp++ = 0;
			else
				*dp++ = i - MINBRT;
		}
	}
					// comptue multiscale local TMO
	if (!compMultiScale(&baseRatio, nlevels, nlevels))
		return false;
					// multiply new ratio image
	logMultIm(&curImg, baseRatio);
	baseRatio.Init();		// don't need this anymore
					// transfer information
	*res->Info() = myMap.inpInfo;
	res->Info()->flags &= ~IIFstonits;
	tagTMO(res->Info(), TMOper[ndx].descr);
					// compute scale factor
	float		sf = scaleMax(curImg);
					// check for grayscale
	if (res->GetCS()->format == IPFy)
		return PmapImage(res->Img(), curImg.GetImg(), sf);
					// convert LogL result to gamma'ed space
	ImgColorSpace	gamSpace;
	if (myMap.tms->flags & TM_F_BW) {
		PcopyCS(&gamSpace, &ICS_Y8);
	} else {
		PcopyCS(&gamSpace, &ICS_Y16);
		sf *= pow(double(TM_BRES)/double(1<<16), myMap.tms->mongam);
	}
	gamSpace.gamma = myMap.tms->mongam;
	if (!curImg.ConvertCS(gamSpace, sf))
		return false;
	DASSERT(res->NComp() == 3);
	if (!res->Init(PICready))	// allocate destination
		return false;
	if (myMap.tms->flags & TM_F_BW) {
		for (y = 0; y < res->Height(); y++)
			PgetRGB24fromY8(res->RowByte(y), curImg.GetRowByte(y),
						res->Width());
		return true;
	}
					// put color back into result
	for (y = 0; y < res->Height(); y++) {
		const unsigned short *	ilp = curImg.GetRowShort(y);
		const uby8 *		cs = ProwPtr(&myMap.rgbImg, y);
		uby8 *			dp = res->RowByte(y);
		for (x = res->Width(); x--; ilp++) {
			i = *cs++ * *ilp / myMap.tms->cdiv[RED];
			*dp++ = (i>255) ? 255 : i;
			i = *cs++ * *ilp / myMap.tms->cdiv[GRN];
			*dp++ = (i>255) ? 255 : i;
			i = *cs++ * *ilp / myMap.tms->cdiv[BLU];
			*dp++ = (i>255) ? 255 : i;
		}
	}
	return true;
}

#undef LTLEN
#undef MINLSIZE
#undef LVLRATIO

// Table of available tone-mapping operators
#define TMONONE	(-1)		// no TMO
#define IHISTO	0		// index for histogram TMO (must be zero)
#define	IBILAT	4		// index for bilateral filter TMO
#define IMULTS	5		// index for multiscale local TMO
#define NTMOS	6		// number of tone-mapping operators
#define RHISTO	(-2)		// reader does histogram TMO

TMOperator	TMOper[NTMOS] = {
	{"histogram", "global histogram adjustment", &TMOglobal, &mapHisto},
	{"map", "global curve from stdin", &TMOglobal, &mapLoad},
	{"photographic", "global photographic TMO", &TMOglobal, &mapPhoto},
	{"gamma", "global gamma curve", &TMOglobal, &mapGamma},
	{"bilateral", "bilateral filter TMO", &TMObilat, &mapGamma},
	{"multiscale", "multiscale local TMO", &TMOmultis, &mapGamma},
};

// Print out TMO usage message
static void
tonemapUsage(FILE *fout)
{
	fputs("\tRecognized TMOs:\n", fout);
	for (int i = 0; i < NTMOS; i++)
		fprintf(fout, "\t\t%-16s - %-32s\n",
				TMOper[i].name, TMOper[i].descr);
	fputs("\tBilateral and multiscale may be combined with other operators,\n", fout);
	fputs("\te.g, 'bilateral+histogram' or 'photo+multis'\n", fout);
	fputs("\tModifier 'human' may be applied as well, e.g. 'human+histo'\n", fout);
}

// Set tone-mapping operator based on descriptive string
static bool
tonemapMethod(int *tmop, const char *name)
{
	const char *	cp;
					// handle concatenated specifications
	while ((cp = strpbrk(name, "+/-,: ")) != NULL) {
		char	buf[16];
		strncpy(buf, name, cp-name);
		buf[cp-name] = '\0';
		if (!tonemapMethod(tmop, buf))
			return false;
		name = cp+1;
	}
	const int	nlen = strlen(name);
	int		i;
	if (nlen < 1)
		goto nonesuch;
	if (!strncasecmp(name, "human", nlen)) {
		do_human = true;	// modifier
		return true;
	}
	for (i = 0; i < NTMOS; i++)
		if (!strncasecmp(name, TMOper[i].name, nlen)) {
			if ((i != IBILAT) & (i != IMULTS))
				TMOper[IBILAT].mptr = TMOper[IMULTS].mptr =
						TMOper[i].mptr;
			if ((*tmop != IBILAT) & (*tmop != IMULTS))
				*tmop = i;
			return true;
		}
nonesuch:
	DMESGF(DMCinput, "Unknown tone-mapping operator: %s", name);
	tonemapUsage(stderr);
	return false;
}

// Read chromaticity coordinates from inp string
static bool
getChroma(float chr[2], const char *inp, const char *name)
{
	inp += strlen(name)+1;

	if (sscanf(inp, "%f,%f", &chr[0], &chr[1]) == 2 &&
			(-2.f <= chr[0]) & (chr[0] <= 3.f) &
			(-2.f <= chr[1]) & (chr[1] <= 3.f))
		return true;

	sprintf(dmessage_buf, "Illegal %s chromaticity: %s", name, inp);
	DMESG(DMCinput, dmessage_buf);
	return false;
}

// Print usage message for color space options
static void
colorUsage(FILE *fout)
{
	fputs("\tRecognized color modes:\n", fout);
	fputs("\t\tsRGB           - Standard 24-bit RGB\n", fout);
	fputs("\t\tAdobeRGB       - Adobe 1998 RGB\n", fout);
	fputs("\t\tP3             - P3 gamut RGB\n", fout);
	fputs("\t\tXYZ            - CIE XYZ\n", fout);
	fputs("\t\tgray           - grayscale (Y)\n", fout);
	fputs("\t\tfloat          - floating-point (HDR)\n", fout);
	fputs("\t\tshort          - unsigned short (16-bit)\n", fout);
	fputs("\t\tbyte           - unsigned byte (8-bit)\n", fout);
	fputs("\t\tgamma=VAL      - set output gamma value\n", fout);
	fputs("\t\tlinear         - synonym for 'gamma=1'\n", fout);
	fputs("\t\tR=x,y          - set red chromaticity\n", fout);
	fputs("\t\tG=x,y          - set green chromaticity\n", fout);
	fputs("\t\tB=x,y          - set blue chromaticity\n", fout);
	fputs("\t\tW=x,y          - set white chromaticity\n", fout);
	fputs("\tCombine with space or '+', e.g., 'short+AdobeRGB'\n", fout);
	fputs("\tDefault is linear float with CCIR-709 (RGB) primaries\n", fout);
}

// Alter color space according to specified pixel type
static bool
colorSpace(ImgColorSpace *csp, const char *typ)
{
	const char *	cp;
					// handle concatenated specifications
	while ((cp = strpbrk(typ, "+/: ")) != NULL) {
		char	buf[32];
		strncpy(buf, typ, cp-typ);
		buf[cp-typ] = '\0';
		if (!colorSpace(csp, buf))
			return false;
		typ = cp+1;
	}
	if (!strcasecmp(typ, "sRGB")) {
		PcopyCS(csp, &ICS_sRGB);
		return true;
	}
	if (!strcasecmp(typ, "AdobeRGB")) {
		csp->format = IPFrgb;
		memcpy(csp->chroma, ICS_RGB98Adobe.chroma, sizeof(csp->chroma));
		return true;
	}
	if (!strcasecmp(typ, "P3")) {
		csp->format = IPFrgb;
		memcpy(csp->chroma, ICS_P3.chroma, sizeof(csp->chroma));
		return true;
	}
	if (!strcasecmp(typ, "XYZ")) {
		csp->format = IPFxyz;
		memcpy(csp->chroma, ICS_XYZ.chroma, sizeof(csp->chroma));
		return true;
	}
	if (!strcasecmp(typ, "gray") || !strcasecmp(typ, "grey") ||
			!strcasecmp(typ, "Y")) {
		csp->format = IPFy;
		return true;
	}
	if (!strcasecmp(typ, "float") || !strcasecmp(typ, "HDR")) {
		csp->dtype = IDTfloat;
		if (csp->gamma >= MIN_GAMMA)
			csp->gamma = 1.f;
		return true;
	}
	if (!strcasecmp(typ, "short") || !strcasecmp(typ, "16-bit")) {
		csp->dtype = IDTushort;
		return true;
	}
	if (!strcasecmp(typ, "byte") || !strcasecmp(typ, "8-bit")) {
		csp->dtype = IDTubyte;
		if (csp->gamma < MIN_GAMMA)
			csp->gamma = DEF_GAMMA;
		return true;
	}
	if (!strcasecmp(typ, "linear")) {
		csp->gamma = 1.f;
		return true;
	}
	if (!strncasecmp(typ, "gamma=", 6)) {
		csp->gamma = atof(typ+6);
		if ((csp->gamma <= 0) | (csp->gamma > 10)) {
			DMESGF(DMCinput, "Illegal gamma: %s", typ+6);
			return false;
		}
		return true;
	}
	if (!strncasecmp(typ, "R=", 2)) {
		csp->format = IPFrgb;
		return getChroma(csp->chroma[0], typ, "R");
	}
	if (!strncasecmp(typ, "G=", 2)) {
		csp->format = IPFrgb;
		return getChroma(csp->chroma[1], typ, "G");
	}
	if (!strncasecmp(typ, "B=", 2)) {
		csp->format = IPFrgb;
		return getChroma(csp->chroma[2], typ, "B");
	}
	if (!strncasecmp(typ, "W=", 2))
		return getChroma(csp->chroma[3], typ, "W");

	DMESGF(DMCinput, "Unknown conversion option: %s", typ);
	colorUsage(stderr);
	return false;
}

// Compute percentile maximum in all channels
static float
percentileMax(const ImgStruct *img, const float pct)
{
	ImgColorSpace	myCS;
	PixelVal	pres;

	PcopyCS(&myCS, img->csp);
	myCS.format = IPFy;
	myCS.chroma[0][0] = myCS.chroma[0][1] = 1./3.;
	pres.csp = &myCS;
	if (!PcomputePercentiles(pres.v.b, &pct, 1, img, PHCrgb))
		return 1.f;
	return PgetY(pres);
}

// Check that view settings are consistent
static void
checkView(ImgInfo *ip)
{
	char	pView[256];
	if (GetImgInfoParam(ip, "VIEW", pView)) {
		if (!(ip->flags & IIFview)) ip->view[0] = '\0';
		if (pView[0] != ' ') strcat(ip->view, " ");
		strcat(ip->view, pView);
		ip->flags |= IIFview;
		SetImgInfoParam(ip, "VIEW", NULL);
	}
	if ((ip->flags & (IIFhvangle|IIFview)) == (IIFhvangle|IIFview)) {
		char *	vs = strstr(ip->view, "-vh ");
		if (vs != NULL)
			ip->hvangle = atof(vs+4);
		if (strstr(ip->view, "-vtl") != NULL)
			ip->hvangle = 0;
	}
}

// Swap bytes in a short, int, or long word
static inline void
swapBytes(uby8 *bytes, int n)
{
	uby8	t;
	switch (n) {
	case 2:
		t = bytes[0]; bytes[0] = bytes[1]; bytes[1] = t;
		break;
	case 4:
		t = bytes[0]; bytes[0] = bytes[3]; bytes[3] = t;
		t = bytes[1]; bytes[1] = bytes[2]; bytes[2] = t;
		break;
	case 8:
		t = bytes[0]; bytes[0] = bytes[7]; bytes[7] = t;
		t = bytes[1]; bytes[1] = bytes[6]; bytes[6] = t;
		t = bytes[2]; bytes[2] = bytes[5]; bytes[5] = t;
		t = bytes[3]; bytes[3] = bytes[4]; bytes[4] = t;
		break;
	default:
		DMESG(DMCparameter, "Illegal data size in swapBytes");
	}
}

// Load a RAW image, swapping bytes if necessary
static bool
loadRAW(PanImage *imp, const char *fname, bool need2swap = false)
{
	int	fd = open(fname, O_RDONLY);

	if (fd < 0) {
		DMESGF(DMCresource, "Cannot open '%s' for reading", fname);
		return false;
	}
	SET_FD_BINARY(fd);
	const int	wsiz = ImgDataSize[imp->GetCS()->dtype];
	const int	rowLen = wsiz*imp->NComp()*imp->Width();
	off_t		flen = lseek(fd, 0, SEEK_END);
	off_t		fpos = (off_t)rowLen*imp->Height();
	if ((flen < fpos) | (flen - fpos > 4096)) {
		sprintf(dmessage_buf,
	"RAW input '%s' unexpected length (%ld rather than %ld bytes)",
				fname, (long)flen, (long)fpos);
		DMESG(DMCdata, dmessage_buf);
		close(fd);
		return false;
	}
	if ((fpos = flen - fpos)) {		// size of header to skip
		sprintf(dmessage_buf, "Skipping first %ld bytes in RAW input '%s'",
				(long)fpos, fname);
		DMESG(DMCwarning, dmessage_buf);
	}
	if (!imp->Init(PICready)) {		// allocate destination
		close(fd);
		return false;
	}
	if (lseek(fd, fpos, SEEK_SET) < 0) {	// read the RAW input
		DMESGF(DMCresource, "Seek error on input '%s'", fname);
		return false;
	}
	need2swap &= (wsiz > 1);
	sprintf(dmessage_buf, "Loading %dx%d %s%s RAW image from '%s'",
			imp->Width(), imp->Height(),
			need2swap ? "swapped " : "",
			PdescribeCS(imp->GetCS(),NULL), fname);
	DMESG(DMCtrace, dmessage_buf);
	for (int y = 0; y < imp->Height(); y++) {
		uby8 *	dst = imp->Row(y);
		if (read(fd, dst, rowLen) != rowLen) {
			DMESGF(DMCresource, "Error reading input '%s'", fname);
			close(fd);
			return false;
		}
		if (need2swap) {
			uby8 *	bp = dst + rowLen;
			while (bp > dst)
				swapBytes(bp -= wsiz, wsiz);
		}
	}
	close(fd);
	return true;
}

// Write a RAW image, swapping bytes first if needed
static bool
writeRAW(const char *fname, const PanImage &im, bool need2swap = false)
{
	int			fd = open(fname, O_WRONLY|O_CREAT|O_TRUNC, 0666);
	
	if (fd < 0) {
		DMESGF(DMCresource, "Cannot open '%s' for writing", fname);
		return false;
	}
	SET_FD_BINARY(fd);
	const int	wsiz = ImgDataSize[im.GetCS()->dtype];
	const int	rowLen = wsiz*im.NComp()*im.Width();
	uby8 *		buf = NULL;
	int		y;
	if ((need2swap &= (wsiz > 1)))
		buf = new uby8 [rowLen];
	sprintf(dmessage_buf, "Writing %dx%d %s%s RAW image to '%s'",
			im.Width(), im.Height(),
			need2swap ? "swapped " : "",
			PdescribeCS(im.GetCS(),NULL), fname);
	DMESG(DMCinfo, dmessage_buf);
	for (y = 0; y < im.Height(); y++)
		if (need2swap) {
			memcpy(buf, im[y], rowLen);
			uby8 *	bp = buf + rowLen;
			while (bp > buf)
				swapBytes(bp -= wsiz, wsiz);

			if (write(fd, buf, rowLen) != rowLen)
				break;
		} else if (write(fd, im[y], rowLen) != rowLen)
			break;

	delete [] buf;
	close(fd);
	if (y < im.Height()) {
		DMESGF(DMCresource, "Error writing output '%s'", fname);
		return false;
	}
	return true;
}

// Convert between image types and formats
int
main(int argc, char * argv[])
{
	bool		rem_flare = false;
	float		expos = 1.f;
	float		maxvalue = -101.f;
	float		set_s2nits = -1.f;
	int		qual = -1;
	double		scale = 0;
	float		dilate = 0;
	int		xres = 0, yres = 0;
	int		tmIndex = TMONONE;
	float		white_bal = 0;
	int		orient = -1;
	const char *	match_histo = NULL;
	bool		need2swap = false;
	bool		setInpCS = false;
	bool		matchInput = false;
	char		outComments[2048];
	char *		commp = outComments;
	char		outParams[2048];
	char *		paramp = outParams;
	ImgColorSpace	inpCS, outCS;
	float		blur = 0;
						// set progname
	fixargv0(argv[0]);
						// Initialize Pancine i/o
	PloadStandardReaders();
	PaddIReaderI(&IRInterfaceDPT);
	PaddIReaderI(&IRInterfaceMTX);
	PloadStandardWriters();
	PanImage::AddImageWriter(&IWInterfaceDPT);
	PanImage::AddImageWriter(&IWInterfaceMTX);
	PanImage::defMathFlags = PIMstrictcolor|PIMwarn2error;
	PcopyCS(&inpCS, &ICS_RGB709);
	PcopyCS(&outCS, &ICS_RGB709);
						// turn on warnings
	dmessage_class_flags[DMCinfo] |= DMFstderr;
	dmessage_class_flags[DMCwarning] |= DMFstderr;
						// get options
	while (argc > 3 && argv[1][0] == '-') {
		const int	alen = strlen(argv[1]);
		if (alen < 2)
			break;
		if (argc > 3 && !strncasecmp(argv[1], "-verbose", alen)) {
			dmessage_class_flags[DMCtrace] |= DMFstderr;
		} else if (argc > 4 && !strncasecmp(argv[1], "-quality", alen)) {
			qual = atoi(argv[2]);
			if (qual < 0) qual = 0;
			else if (qual > 100) qual = 100;
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-blur", alen)) {
			blur = atof(argv[2]);
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-scale", alen)) {
			scale = atof(argv[2]);
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-rotate", alen)) {
			switch (atoi(argv[2])) {
			case 0: orient = 0; break;
			case 90: orient = PIOcw90; break;
			case -90: case 270: orient = PIOccw90; break;
			case 180: case -180: orient = PIO180; break;
			default:
				DMESG(DMCinput, "Illegal rotation angle");
				return 1;
			}
			const char *	cp = argv[2];
			cp += (*cp == '-');
			while (*cp && isdigit(*cp)) cp++;
			if (toupper(*cp) == 'H')	// horizontal flip
				switch (orient) {
				case 0: orient = PIOhflip; break;
				case PIOcw90: orient = PIOcw90hflip; break;
				case PIOccw90: orient = PIOcw90vflip; break;
				case PIO180: orient = PIOvflip; break;
				}
			else if (toupper(*cp) == 'V')	// vertical flip
				switch (orient) {
				case 0: orient = PIOvflip; break;
				case PIOcw90: orient = PIOcw90vflip; break;
				case PIOccw90: orient = PIOcw90hflip; break;
				case PIO180: orient = PIOhflip; break;
				}
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-expose", alen)) {
			expos = atof(argv[2]);
			if (argv[2][0] == '+' || expos <= 0)
				expos = pow(2., expos);
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-calibrate", alen)) {
			maxvalue = atof(argv[2]);
			if (maxvalue <= .0f) {
				DMESG(DMCinput, "Illegal calibration value");
				return 1;
			}
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-limit", alen)) {
			maxvalue = -atof(argv[2]);
			if ((maxvalue <= -100.f) | (maxvalue > .0f)) {
				DMESG(DMCinput, "Illegal percentile limit");
				return 1;
			}
			++argv; --argc;
		} else if (argc > 5 && !strncasecmp(argv[1], "-dim", alen)) {
			xres = atoi(argv[2]);
			yres = atoi(argv[3]);
			argv += 2; argc -= 2;
		} else if (argc > 4 && !strncasecmp(argv[1], "-dilate", alen)) {
			dilate = atof(argv[2]);
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-erode", alen)) {
			dilate = -atof(argv[2]);
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-xres", alen)) {
			xres = atoi(argv[2]);
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-yres", alen)) {
			yres = atoi(argv[2]);
			++argv; --argc;
		} else if (argc > 3 && !strncasecmp(argv[1], "-swap", alen)) {
			need2swap = true;
		} else if (argc > 4 && !strncasecmp(argv[1], "-input", alen)) {
			if (!colorSpace(&inpCS, argv[2]))
				return 1;
			setInpCS = true;
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-output", alen)) {
			matchInput = !strncasecmp(argv[2], "in", 2);
			if (!matchInput && !colorSpace(&outCS, argv[2]))
				return 1;
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-tmo", alen)) {
			if (!tonemapMethod(&tmIndex, argv[2]))
				return 1;
			++argv; --argc;
		} else if (argc > 3 && !strncasecmp(argv[1], "-flare", alen)) {
			rem_flare = true;
		} else if (argc > 4 && !strncasecmp(argv[1], "-wb", alen)) {
			white_bal = atof(argv[2]);
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-match", alen)) {
			match_histo = argv[2];
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-comment", alen)) {
			strcpy(commp, argv[2]);
			while (*commp) ++commp;
			*commp++ = '\n';
			++argv; --argc;
		} else if (argc > 4 && !strncasecmp(argv[1], "-parameter", alen)
					&& isalpha(argv[2][0])
					&& strchr(argv[2], '=')) {
			strcpy(paramp, argv[2]);
			while (*paramp) {
				if (*paramp == ';') *paramp = '\n';
				++paramp;
			}
			if (paramp[-1] != '\n') *paramp++ = '\n';
			++argv; --argc;
		} else
			break;
		++argv; --argc;
	}
	*commp = '\0';
	*paramp = '\0';
	if (argc != 3) {
		fprintf(stdout, USE_MSG, progname);
		fputc('\n', stdout);
		tonemapUsage(stdout);
		fputc('\n', stdout);
		colorUsage(stdout);
		return 1;
	}
	if (tmIndex != TMONONE) {
		if (match_histo != NULL) {
			DMESG(DMCinput, "Cannot match histogram while tone-mapping");
			return 1;
		}
		expos = 1.f;
		maxvalue = -101.f;
	}
	const char *	inpfn = argv[1];
	const char *	outfn = argv[2];
	const char *	inpsfx = PgetSuffix(inpfn);
	const char *	outsfx = PgetSuffix(outfn);
	if (inpsfx == NULL) {
		DMESGF(DMCinput, "Input file name '%s' missing suffix", inpfn);
		return 1;
	}
	if (outsfx == NULL) {
		DMESGF(DMCinput, "Output file name '%s' missing suffix", outfn);
		return 1;
	}
	ImgReader *	irp = NULL;		// check input image
	int		in_xres = xres;
	int		in_yres = yres;
	if (!strcasecmp(inpsfx, "RAW")) {	// RAW input?
		if ((in_xres <= 0) | (in_yres <= 0)) {
			DMESG(DMCinput, "Need dimensions for RAW input");
			return 1;
		}
		set_s2nits = maxvalue*(maxvalue > .0f);
	} else {				// else open image reader
		irp = PopenImageF(inpfn, false, NULL);
		if (irp == NULL)
			return 1;
		in_xres = irp->xres;
		in_yres = irp->yres;
		if (irp->pixAspect > 1.01)
			in_yres = (int)(in_yres/irp->pixAspect + .5);
		else if (irp->pixAspect < 0.99)
			in_xres = (int)(in_xres*irp->pixAspect + .5);
		if (setInpCS) {			// munge input color space?
			if ((ImgPixelLen[inpCS.format] != ImgPixelLen[irp->cs.format]) |
					(inpCS.dtype != irp->cs.dtype)) {
				DMESG(DMCinput, "Cannot munge input color space -- incompatible pixel type");
				return 1;
			}
			PcopyCS(&irp->cs, &inpCS);
		} else
			PcopyCS(&inpCS, &irp->cs);
		if (maxvalue > .0f) {		// calibrate output
			ImgInfo	myInfo;
			if (IRgetInfo(irp, &myInfo) == IREnone &&
					myInfo.flags & IIFstonits)
				expos *= myInfo.stonits / maxvalue;
			else
				set_s2nits = maxvalue;
		}
	}
	if (scale > 0) {			// adjust output dimensions
		xres = (int)(scale*in_xres + .5);
		yres = (int)(scale*in_yres + .5);
	} else if ((xres <= 0) & (yres <= 0)) {
		xres = in_xres;
		yres = in_yres;
	} else if (yres <= 0 || (xres > 0 &&
			(xres+1)*in_yres < yres*in_xres)) {
		yres = in_yres*xres/in_xres;
	} else if (xres <= 0 || xres*in_yres > (yres+1)*in_xres) {
		xres = in_xres*yres/in_yres;
	}
	if (tmIndex >= 0) {			// convert primaries & gamma?
		if (inpCS.dtype == IDTubyte) {
			DMESG(DMCinput, "Tone-mapping requires HDR input");
			return 1;
		}
		if (matchInput) {
			DMESG(DMCinput, "Cannot match color space with tone-mapping");
			return 1;
		}
		outCS.dtype = IDTubyte;
		if (outCS.gamma < MIN_GAMMA)
			outCS.gamma = DEF_GAMMA;
	} else if (matchInput) {		// match input color space
		PcopyCS(&outCS, &inpCS);
	} else if ((irp != NULL) & (inpCS.dtype == outCS.dtype)) {
		inpCS.gamma = outCS.gamma;	// do some stuff on load
		if (inpCS.format == outCS.format)
			memcpy(inpCS.chroma, outCS.chroma, sizeof(inpCS.chroma));
	}
						// check tone-mapping
	if (rem_flare | (tmIndex > IHISTO)) {
		if (inpCS.dtype == IDTushort) {
			rem_flare &= (irp != NULL);
			if (rem_flare) {	// upconvert in Load()
				inpCS.dtype = IDTfloat;
				inpCS.gamma = 1.f;
			}
		} else if (inpCS.dtype != IDTfloat) {
			DMESGF(DMCinput, "%s requires HDR input",
					rem_flare ? "Flare removal" : "TMO");
			return 1;
		}
	} else if (tmIndex == IHISTO && (xres == in_xres) & (yres == in_yres)
			&& white_bal <= .01f
			&& irp != NULL && irp->ri->ToneMapping != NULL) {
		DMESG(DMCinfo, "Using image reader's tone-mapping");
		PcopyCS(&inpCS, &outCS);
		IRtoneMapping(irp, .0, .0, do_human);
		tmIndex = RHISTO;		// use as reminder below
	}
	PanImage	matchImg;		// load image to match?
	if (match_histo != NULL && !matchImg.Load(match_histo))
		return 1;
	if (irp != NULL) {			// resample during load?
		if ((inpCS.dtype != IDTubyte) | (outCS.dtype != IDTfloat)) {
			in_xres = xres; in_yres = yres;
		} else if ((xres < in_xres) & (yres < in_yres)) {
			in_xres = xres; in_yres = yres;
			if (!PmatchColorSpace(matchImg.GetCS(), &inpCS, PICMptype))
				PcopyCS(&inpCS, &outCS); // convert during load, too
		}
	} else if (tmIndex >= 0 && (xres != in_xres) | (yres != in_yres)) {
		DMESG(DMCparameter, "Cannot simultaneously resize and tone-map RAW input");
		return 1;
	}
						// load our image
	PanImage	inpImg(in_xres, in_yres, inpCS, PICfree);
	if (irp != NULL) {
		if (!inpImg.Load(irp))
			return 1;
		IRclose(irp);			// done with reader
		if ((in_xres == xres) & (in_yres == yres)) {
			xres = in_xres = inpImg.Width();
			yres = in_yres = inpImg.Height();
		}
	} else if (!loadRAW(&inpImg, inpfn, need2swap))
		return 1;
						// reorient image?
	if (orient > 0 || (!orient && inpImg.GetInfo(IIForientation) &&
				(orient = inpImg.GetInfo()->orientation) > 1)) {
		if (!inpImg.Reorient(orient)) {
			DMESG(DMCparameter, "Image orientation failure");
			return 1;
		}
		#ifndef YMAJOR
		#define YMAJOR 4		// swapped x & y axes?
		#endif
		if ((RHortab[0] ^ RHortab[orient-1]) & YMAJOR) {
			int	t = xres;
			xres = yres; yres = t;
		}
	}
						// base exposure on limit?
	if ((-100.f < maxvalue) & (maxvalue <= .0f))
		expos /= percentileMax(inpImg.GetImg(), 100.f+maxvalue);
						// convert 8-bit to float?
	if ((inpCS.dtype == IDTubyte) & (outCS.dtype == IDTfloat)) {
		PanImage	tmpImg(inpImg.Width(), inpImg.Height(), outCS);
		DMESG(DMCinfo, "Dequantizing 8-bit input");
#if 1
		if (!PdequantizeImage2(tmpImg.Img(), inpImg.GetImg(), expos, 0, 0))
			return 1;
#else
		if (!Pdecontour(&tmpImg, inpImg, expos))
			return 1;
#endif
		inpImg.GetInfo(tmpImg.Info());
		if (tmpImg.GetInfo(IIFstonits))
			tmpImg.Info()->stonits /= expos;
		expos = 1.f;
		if (!inpImg.Take(&tmpImg))
			return 1;
		PcopyCS(&inpCS, inpImg.GetCS());
	}
						// perform white balance
	if (white_bal > 0.01f && !PautoWhiteBal(inpImg.Img(), white_bal, true))
		return 1;
						// perform flare removal
	if (rem_flare)
		PHDremoveFlare(inpImg.Img());
						// final TMO/resize/conversion
	PanImage	outImg(xres, yres, outCS, PICfree);
	if (tmIndex >= 0) {
		if (!ApplyTMO(&outImg, &inpImg, tmIndex))
			return 1;
	} else {
		int	when2expos = (0.999f <= expos) & (expos <= 1.001f) ? 0 :
					(ImgDataSize[inpCS.dtype] >
					    ImgDataSize[outCS.dtype]) ? -1 : 1;
		int	when2match = !match_histo ? 0 :
				(inpCS.dtype == IDTfloat ||
				 PmatchColorSpace(matchImg.GetCS(), &inpCS,
							PICMptype)) ? -1 : 1;
		when2expos *= !when2match;
		if (when2expos < 0)
			inpImg *= expos;
		if (when2match < 0 && !PmatchHisto(inpImg.Img(), matchImg.GetImg()))
			return 1;
		if ((xres == in_xres) & (yres == in_yres) &&
				PmatchColorSpace(&inpCS, &outCS, PICMall)) {
			if (!outImg.Take(&inpImg))
				return 1;
		} else if (!outImg.Load(inpImg))
			return 1;

		inpImg.Init();
		if (when2expos > 0)
			outImg *= expos;
		if (when2match > 0) {
			if (outCS.dtype != IDTfloat &&
					!PmatchColorSpace(matchImg.GetCS(), &outCS, PICMptype) &&
					!matchImg.ConvertCS(outCS))
				return 1;
			if (!PmatchHisto(outImg.Img(), matchImg.GetImg()))
				return 1;
		}
		matchImg.Init();
		if (tmIndex == RHISTO)		// reminder flag from above
			tagTMO(outImg.Info(), TMOper[IHISTO].descr);
		expos = 1.f;
#if 0
// hack 6-bit output	(use 8 for threshold in PdequantizeImage2 as well!)
for (int y=outImg.Height()*(outCS.dtype==IDTubyte); y--; ) {
uby8 *	pp = outImg.RowByte(y);
for (int n=outImg.NComp()*outImg.Width(); n--; pp++)
*pp = (*pp & 0xfc) | 0x2;
}
#endif
	}
	if (match_histo != NULL)		// matched histogram?
		sprintf(strchr(outImg.Info(IIFcomments)->comments,'\0'),
				"Matched histogram to '%s'\n", match_histo);

	if (dilate != 0)			// dilate or erode image?
		outImg <<= dilate;

	if (blur != 0)				// blur or sharpen image?
		outImg ^= blur;

	if (set_s2nits > 0)			// override sample-to-nits?
		outImg.Info(IIFstonits)->stonits = set_s2nits;

	if (outComments[0])			// add output comments?
		strlcat(outImg.Info(IIFcomments)->comments, outComments,
						sizeof(defImgInfo.comments));

	if (outParams[0])			// add output parameters?
		strlcat(outImg.Info(IIFparams)->params, outParams,
						sizeof(defImgInfo.params));
	if (outCS.dtype == IDTfloat &&		// hack for MS in JPEG-HDR
			!strcasecmp(outsfx, IWInterfaceJPEG.suffix))
		strlcat(outImg.Info(IIFparams)->params, "JHtmo=multiscale\n",
						sizeof(defImgInfo.params));
	checkView(outImg.Info());		// make sure view is consistent

	if ( strcasecmp(outsfx, "RAW") ?	// write out final image
			!outImg.Write(outfn, qual) :
			!writeRAW(outfn, outImg, need2swap) )
		return 1;
	return 0;
}
