/*
 *  warpmain.c
 *
 *  Warp an image (or sequence) according to provided source grid.
 *
 *  Created by Greg Ward on Fri June 6, 2023.
 *  Copyright (c) 2023 Anyhere Software. All rights reserved.
 *
 */

#include <stdlib.h>
#include <ctype.h>
#include "system.h"
#include "rtio.h"
#include "pimage.h"
#include "imgwriter.h"
#include "dmessage.h"

void
usage_error(const char *pnam)
{
	fputs("Usage: ", stderr);
	fputs(pnam, stderr);
	fputs(" [-v][-q qual][-x hres][-y vres][-dim hres vres][-s start][-e end] \\\n", stderr);
	fputs("\t{ wgrid.txt | - } wg_width wg_height inspec.img outspec.img\n", stderr);
	fputs("Where grid has uv source pairs, English order, normalized to [0,1] range\n", stderr);
	exit(1);
}

/* Check for %d in output format string */
int
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
	return 0;
}

/* Load grid positions */
float *
loadGrid(FILE *fp, const int wd, const int ht)
{
	float	*wgarr, *wgp, dummy;
	int	noutside = 0;
	int	n;

	if (!fp | (wd < 2) | (ht < 2))
		return NULL;
	wgarr = (float *)malloc(sizeof(float)*2*wd*ht);
	if (!wgarr) {
		DMESG(DMCmemory, "Out of memory in loadGrid()");
		return NULL;
	}
	for (n = ht*wd*2, wgp = wgarr; n--; wgp++) {
	    	if (fscanf(fp, "%f", wgp) != 1) {
	    		DMESGF(DMCdata, "Error reading warp grid at byte %ld", ftell(fp));
	    		free(wgarr);
	    		return NULL;
		}
		noutside += (*wgp < 0) | (*wgp > 1);
	}
	DTEST(fscanf(fp, "%f", &dummy) == 1, DMCwarning,
			"Unexpected value(s) at end of grid data");
	DTESTF(noutside > ht*wd*2/25, DMCwarning, "%.1f%% grid values outside [0,1] range",
			100./(ht*wd*2)*noutside);
	return wgarr;
}

/* Apply a normalized warp grid to image */
int
warpImage(ImgStruct *ib, const ImgStruct *ia, float *grid, const int wd, const int ht)
{
	float		*ourGrid, *wgp;
	const float	*igp;
	int		i, j;

	if (!ib || !ib->csp || !ia || !ia->csp | !ia->img || !grid)
		return 0;
	ourGrid = (float *)malloc(sizeof(float)*2*wd*ht);
	if (!ourGrid) {
		DMESG(DMCmemory, "Cannot allocate temp array in warpImage()");
		return 0;
	}
	wgp = ourGrid;			/* convert normalized source positions to pixels */
	igp = grid;
	for (j = 0; j < ht; j++)
	    for (i = 0; i < wd; i++) {
	    	*wgp++ = *igp++ * (float)ia->xres;
	    	*wgp++ = *igp++ * (float)ia->yres;
	    }
	i = PwarpImage(ib, ia, (float (*)[2])ourGrid, wd, ht, PScubic, Pblack);
	free(ourGrid);
	return i;
}

/* Load the specified image file */
int
loadImage(ImgStruct *ib, const char *fname, int quiet, ImgInfo *infop)
{
	ImgReader	*ir = PopenImageF(fname, quiet, NULL);
	int		ok;

	if (!ir) return 0;

	if (!ib->csp) ib->csp = &ir->cs;
	if ((ib->xres <= 0) | (ib->yres <= 0)) {
		ib->xres = ir->xres;
		ib->yres = ir->yres;
	}
	ok = PrenderImageR(ib, ir, quiet);
	if (infop) {
		*infop = defImgInfo;
		if (ok) IRgetInfo(ir, infop);
	}
	IRclose(ir);
	return ok;
}

/* Write out the specified image file */
int
saveImage(const ImgStruct *ia, const char *fname, int qual, const ImgInfo *infop)
{
	extern const ImgWriterInterface	IWInterfaceJPEG;
	extern const ImgWriterInterface	IWInterfaceTIFF;
	extern const ImgWriterInterface	IWInterfaceRad;
	extern const ImgWriterInterface IWInterfaceEXR;
	const ImgWriterInterface	*wip = NULL;
	const char			*cp;
	ImgWriteBuf			iwb;

	if (!ia || !ia->csp | !ia->img | !fname)
		return 0;
	if (!(cp = PgetSuffix(fname))) {	/* figure out what we're writing */
		DMESGF(DMCinput, "Output image name '%s' missing suffix", fname);
		return 0;
	}
	if (!strcasecmp(IWInterfaceJPEG.suffix, cp))
		wip = &IWInterfaceJPEG;
	else if (!strcasecmp(IWInterfaceTIFF.suffix, cp))
		wip = &IWInterfaceTIFF;
	else if (!strcasecmp(IWInterfaceRad.suffix, cp))
		wip = &IWInterfaceRad;
	else if (!strcasecmp(IWInterfaceEXR.suffix, cp))
		wip = &IWInterfaceEXR;
	if (!wip) {
		DMESGF(DMCinput, "Unsupported output image type '%s'", cp);
		return 0;
	}
	if ((qual < 0) | (qual > 100)) {	/* get output quality setting */
		if (infop && infop->flags & IIFquality)
			qual = infop->quality;
		else
			qual = DEF_IQUALITY;
	}
	cp = (*wip->SupportedCS)(ia->csp, qual);
	if (!cp) {
		sprintf(dmessage_buf, "Output color space %s not supported by %s writer",
				PdescribeCS(ia->csp,NULL), wip->suffix);
		DMESG(DMCparameter, dmessage_buf);
		return 0;
	}
	sprintf(dmessage_buf, "Writing %s image '%s'", cp, fname);
	DMESG(DMCtrace, dmessage_buf);
	PsetWriteBuf(&iwb, ia);
	if (infop) iwb.info = *infop;
	if (qual >= 0) {
		iwb.info.quality = qual;
		iwb.info.flags |= IIFquality;
	}
	if ((*wip->WriteImage)(fname, &iwb) <= 0)
		DMESGF(DMCresource, "Error writing output image '%s'", fname);
	return 1;
}

/* Warp input images and write out results */
int
main(int argc, char *argv[])
{
	int	oQual = -1;
	int	fstart=1, fend=0;
	int	xres=0, yres=0;
	int	gwidth=0, gheight=0;
	float	*grid = NULL;
	FILE	*fp = NULL;
	int	doSeq, fn;
	int	a;
					/* most errors call exit() */
	dmessage_class_flags[DMCsystem] |= DMFexit;
	dmessage_class_flags[DMCparameter] |= DMFexit;
	dmessage_class_flags[DMCresource] |= DMFexit;
	dmessage_class_flags[DMCdata] |= DMFexit;
	PloadStandardReaders();
					/* process options */
	for (a = 1; a < argc && argv[a][0] == '-'; a++) {
		switch (argv[a][1]) {
		case 'v':			/* verbose mode */
			dmessage_class_flags[DMCinfo] |= DMFstderr;
			dmessage_class_flags[DMCtrace] |= DMFstderr;
			continue;
		case 'q':			/* output quality */
			oQual = atoi(argv[++a]);
			continue;
		case 'x':			/* output width */
			xres = atoi(argv[++a]);
			continue;
		case 'y':			/* output height */
			yres = atoi(argv[++a]);
			continue;
		case 'd':			/* width and height */
			xres = atoi(argv[++a]);
			yres = atoi(argv[++a]);
			continue;
		case 's':			/* start frame # */
			fstart = atoi(argv[++a]);
			continue;
		case 'e':			/* end frame # */
			fend = atoi(argv[++a]);
			continue;
		case '\0':			/* grid from stdin */
			fp = stdin;
			break;
		default:
			DMESGF(DMCinput, "Unknown option: %s", argv[a]);
			usage_error(argv[0]);
		}
		break;
	}
	if (a != argc-5)
		usage_error(argv[0]);
	gwidth = atoi(argv[a+1]);	/* load warp grid */
	gheight = atoi(argv[a+2]);
	if ((gwidth <= 1) | (gheight <= 1))
		usage_error(argv[0]);
	if (!fp && !(fp = fopen(argv[a], "r"))) {
		DMESGF(DMCresource, "Cannot open '%s'", argv[a]);
		return 1;
	}
	grid = loadGrid(fp, gwidth, gheight);
	if (!grid)
		return 1;
	if (fp != stdin)
		fclose(fp);
	a += 3;				/* check i/o image specs */
	doSeq = hasFrameFormat(argv[a]);
	if (doSeq ^ hasFrameFormat(argv[a+1])) {
		DMESG(DMCinput, "Input spec has %d iff output has %d");
		usage_error(argv[0]);
	}
	if (!doSeq) fend = fstart+1;	/* process image(s) */
	for (fn = fstart; (fend <= 0) | (fn < fend); fn++) {
		ImgStruct	imgIn, imgOut;
		ImgInfo		iInfo;
		char		fname[256];
		memset(&imgIn, 0, sizeof(imgIn));
		sprintf(fname, argv[a], fn);
		if (!loadImage(&imgIn, fname, fn>fstart, &iInfo)) {
			DTEST((fend>0)&(fn>fstart),
				DMCwarning, "Ran out of input frames");
			break;
		}
		imgOut = imgIn;		/* warp image */
		imgOut.img = NULL;
		if (xres > 0) imgOut.xres = xres;
		if (yres > 0) imgOut.yres = yres;
		if (!warpImage(&imgOut, &imgIn, grid, gwidth, gheight))
			return 1;
		PfreeImage(&imgIn);
		iInfo.flags &= ~(IIFhdensity|IIFhvangle|IIFcrop|IIForientation|IIFview);
		sprintf(fname, argv[a+1], fn);
		if (!saveImage(&imgOut, fname, oQual, &iInfo))
			return 1;
		PfreeImage(&imgOut);
	}
	free(grid);
	return (fn == fstart);
}
