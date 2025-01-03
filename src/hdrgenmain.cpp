/*
 * Pancine-based program for generating high dynamic-range image
 *
 * Created 9/14/01 by Greg Ward
 */

#include "pancine.h"
#include "phdrimg.h"
#include "panimage.h"
#include <fcntl.h>

static int
progressBar(const char *act, int pctdone)
{
	static const char *	prev_msg = "No Message";
	const int	actwidth = 24, barwidth = 54;
	int	i;
	
	if (act == NULL)
		act = prev_msg;
	else
		prev_msg = act;
	const int	actlen = strlen(act);
	for (i = actwidth-actlen; i--; )
		cerr << ' ';
	cerr << act << ' ';
	for (i = barwidth*pctdone/100; i--; )
		cerr << '*';
	if (pctdone < 100)
		cerr << '\r';
	else
		cerr << '\n';
	return 1;
}

static const ImgColorSpace *
csSpec(const char *cname)
{
	if (!strcasecmp(cname, "sRGB"))
		return &ICS_RGB709;
	if (!strcasecmp(cname, "XYZ"))
		return &ICS_XYZ;
	if (!strcasecmp(cname, "AdobeRGB"))
		return &ICS_RGB98Adobe;
	if (!strcasecmp(cname, "P3"))
		return &ICS_P3;
	DMESGF(DMCinput,
"Unknown color space '%s'\n\t(accepted: sRGB, XYZ, AdobeRGB, P3)", cname);
	return NULL;
}

int
main(int argc, const char * argv[])
{
	const char *		outfn = NULL;
	const char *		varfn = NULL;
	const char *		resfn = NULL;
	float			stonits = -1.f;
	bool			removeFlare = false;
	PCacheImage *		cip = NULL;
	int			qual = -1;
	bool			force = false;
	int32			crect[4] = {0,0,0,0};
	char			outParams[2048];
	char *			paramp = outParams;
	PHDImageMaker		mkhdr;
	ImgStruct		outimg, varimg;
	const PCacheImage *	cip1;
	DBRecord		rcomm;
	int			i;
						// Initialize Pancine image i/o
	PloadStandardReaders();
	PloadStandardWriters();
	PanImage::AddImageWriter(&IWInterfaceMTX);
						// turn on progress reporting
	mkhdr.reportProgress = &progressBar;
						// increase cache size
	gCacheSize = 1000L*1024L*1024L;
	outimg.csp = &ICS_RGB709;
						// get arguments
	for (i = 1; i < argc; i++)
		if (argv[i][0] == '-')		// option
			switch (argv[i][1]) {
			case 'v':		// verbose reporting
				dmessage_class_flags[DMCtrace] |= DMFstderr;
				dmessage_class_flags[DMCinfo] |= DMFstderr;
				dmessage_class_flags[DMCwarning] |= DMFstderr;
				break;
			case 'o':		// output image
				outfn = argv[++i];
				break;
			case 'c':		// output color space
				outimg.csp = csSpec(argv[++i]);
				if (outimg.csp == NULL)
					return 1;
				break;
			case 'k':		// variance image
				varfn = argv[++i];
				break;
			case 'r':		// response function file
				resfn = argv[++i];
				break;
			case 's':		// stonits for next image
				stonits = atof(argv[++i]);
				break;
			case 'a':		// toggle alignment
				mkhdr.solveFlags ^= PHDFalignment;
				break;
			case 'e':		// toggle exposure adjustment
				mkhdr.solveFlags ^= PHDFexposure;
				break;
			case 'g':		// toggle ghost removal
				mkhdr.solveFlags ^= PHDFghosting;
				break;
			case 'x':		// toggle skip exposures
				mkhdr.solveFlags ^= PHDFskipexp;
				break;
			case 'f':		// toggle flare removal
				removeFlare = !removeFlare;
				break;
			case 'm':		// cache size (Mbytes)
				gCacheSize = 1024L*1024L*atoi(argv[++i]);
				break;
			case 'q':		// output quality setting
				qual = atoi(argv[++i]);
				if (qual < 0) qual = 0;
				else if (qual > 100) qual = 100;
				break;
			case 'p':		// output parameter
				if (!strchr(argv[++i], '=')) {
					DMESGF(DMCinput,
					"Parameter '%s' missing '=' sign\n",
							argv[i]);
					break;
				}
				strcpy(paramp, argv[i]);
				while (*paramp) ++paramp;
				*paramp++ = '\n';
				break;
			case 'K':		// crop rectangle
				crect[0] = atoi(argv[++i]);
				crect[1] = atoi(argv[++i]);
				crect[2] = atoi(argv[++i]);
				crect[3] = atoi(argv[++i]);
				break;
			case 'F':		// force overwrite?
				force = !force;
				break;
			default:
				DMESGF(DMCinput, "Unknown option: %s", argv[i]);
				return 1;
			}
		else {				// input image
			cip = mkhdr.NewCacheImage();
			if (cip == NULL) {
				DMESG(DMCparameter, "Too many input images");
				return 1;
			}
			if (!cip->SetImage(argv[i]))
				return 1;
			if (stonits > .0f) {
				PDBsetField(&cip->ircd, PDBFstonits, stonits);
				stonits = -1.f;
			}
			if ((crect[2] > crect[0]) & (crect[3] > crect[1]))
				PDBarrayField(&cip->ircd, PDBFcrop, crect, 4);
		}
	*paramp = '\0';
	if (outfn == NULL) {
		DMESG(DMCinput, "No -o output specified");
		return 1;
	}
	if (cip == NULL) {
		DMESG(DMCinput, "No input images specified");
		return 1;
	}
	if (!force) {				// don't clobber output
		int	fd = open(outfn, O_WRONLY|O_CREAT|O_EXCL, 0666);
		if (fd < 0) {
			DMESGF(DMCresource, "Cannot open output '%s'",
						outfn);
			return 1;
		}
		close(fd);
	}
						// read response file
	if (resfn != NULL && mkhdr.ReadResponse(resfn))
		resfn = NULL;
						// make hdr image
	outimg.xres = outimg.yres = 0;
	outimg.img = NULL;
	varimg = outimg;
	if (!(varfn == NULL ? mkhdr.GetHDImage(&outimg) :
			mkhdr.GetHDImage(&outimg, &varimg))) {
		DMESG(DMCdata, "Failure generating HDR image");
		if (!force) unlink(outfn);
		return 1;
	}
						// write response if appropriate
	if (resfn != NULL && !mkhdr.WriteResponse(resfn)) {
		DMESGF(DMCresource, "Cannot write response to '%s'", resfn);
		if (!force) unlink(outfn);
		return 1;
	}
						// get output information
	for (i = 0; (cip1 = mkhdr.GetCacheImage(i)) != NULL; i++)
		rcomm &= cip1->ircd;
	if (qual >= 0)
		PDBsetField(&rcomm, PDBFquality, (int32)qual);
	else
		PDBclearField(&rcomm, PDBFquality);
	PDBsetField(&rcomm, PDBFcapdate, cip->ircd);
	PDBsetField(&rcomm, PDBFstonits, mkhdr.GetStoNits());

	PanImage	wrtIm;			// transfer to PanImage
	if (!wrtIm.Take(&outimg))
		return 1;
	PsetInfo(wrtIm.Info(), rcomm);
	if (*outParams)				// add output parameters
		strcat(wrtIm.Info(IIFparams)->params, outParams);
						// history comments
	const char *	fn;
	char *		comm = wrtIm.Info(IIFcomments)->comments;
	while (*comm) ++comm;
	sprintf(comm, "%s created HDR image from", argv[0]);
	while (*comm) ++comm;
	for (i = mkhdr.GetNI(); i--; )
		if (PDBgetField(mkhdr.GetCacheImage(i)->ircd,PDBFfile).Get(&fn)) {
			sprintf(comm, " '%s'", fn);
			while (*comm) ++comm;
		}
	*comm++ = '\n'; *comm = '\0';
						// done with mkhdr info.
	mkhdr.Done();
						// remove flare?
	if (removeFlare && PHDremoveFlare(wrtIm.Img(), &progressBar)) {
		strcpy(comm, "Removed lens flare\n");
		while (*comm) ++comm;
	}
						// hack for MS in JPEG-HDR
	if (!strcasecmp(PgetSuffix(outfn),IWInterfaceJPEG.suffix))
		strcat(wrtIm.Info(IIFparams)->params, "JHtmo=multiscale\n");
	if (!wrtIm.Write(outfn))		// write hdr image
		return 1;
	if (varimg.img == NULL)			// all done if no -k
		return 0;
	if (!wrtIm.Take(&varimg))
		return 1;
	if (qual >= 0)
		wrtIm.Info(IIFquality)->quality = qual;
	if (*outParams)				// add output parameters
		strcat(wrtIm.Info(IIFparams)->params, outParams);
	if (!wrtIm.Write(varfn))		// write variance image
		return 1;
	return 0;
}
