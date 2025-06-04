/*
 * Load and operate on same-sized 2-D bitmaps
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "system.h"
#include "pimage.h"
#include "abitmap.h"
#include "dmessage.h"

extern "C" int	isflt(char *s);

enum CompareOp {coNone, coLT, coGT, coLE, coGE, coEQ, coNE};

CompareOp	revcomp[] = {coNone, coGE, coLE, coGT, coLT, coEQ, coNE};

float		epsilon = 0;	// comparison epsilon in normalized units

// Print usage message and exit with bad status
void
usage()
{
	fputs("Usage: ", stderr);
	fputs(progname, stderr);
	fputs(" [-f][-v][-e epsilon] result.bmp = expr1 [ LOGIC_OP expr2 .. ]\n", stderr);
	fputs("\twhere an expression may be a bitmap file prefixed by an optional '~',\n", stderr);
	fputs("\tor a comparison between either a pair of grayscale images or\n", stderr);
	fputs("\ta grayscale image and a normalized gray value.\n", stderr);
	fputs("\tThe epsilon (-e) setting for comparisons defaults to 0.\n", stderr);
	fputs("Comparison operators (synonyms): < > <= >= == != (LT GT LE GE EQ NE)\n", stderr);
	fputs("Logical operators (synonyms): & | ^ + - (AND OR XOR)\n", stderr);
	exit(1);
}

// Turn string comparison operator to enum
CompareOp
iscomparison(const char *str)
{
	if (!strcmp(str, "<") | !strcasecmp(str, "LT"))
		return coLT;
	if (!strcmp(str, ">") | !strcasecmp(str, "GT"))
		return coGT;
	if (!strcmp(str, "<=") | !strcasecmp(str, "LE"))
		return coLE;
	if (!strcmp(str, ">=") | !strcasecmp(str, "GE"))
		return coGE;
	if (!strcmp(str, "==") | !strcasecmp(str, "EQ"))
		return coEQ;
	if (!strcmp(str, "!=") | !strcasecmp(str, "NE"))
		return coNE;
	return coNone;
}

// Perform comparison for floats
static inline bool
fcompare(float a, CompareOp co, float b)
{
	switch (co) {
	case coLT:
		return (a < b-epsilon);
	case coGT:
		return (a > b+epsilon);
	case coLE:
		return (a <= b+epsilon);
	case coGE:
		return (a >= b-epsilon);
	case coEQ:
		return (a >= b-epsilon) & (a <= b+epsilon);
	case coNE:
		return (a < b-epsilon) | (a > b+epsilon);
	}
	return false;
}

// Perform comparison for ints
static inline bool
icompare(int a, CompareOp co, int b, int ieps)
{
	switch (co) {
	case coLT:
		return (a < b-ieps);
	case coGT:
		return (a > b+ieps);
	case coLE:
		return (a <= b+ieps);
	case coGE:
		return (a >= b-ieps);
	case coEQ:
		return (a >= b-ieps) & (a <= b+ieps);
	case coNE:
		return (a < b-ieps) | (a > b+ieps);
	}
	return false;
}

// Compare values in two buffers and write results to bitmap
bool
compare_buf(ABitMap2 *mp, const ImgReadBuf *bleft,
			const CompareOp co, const ImgReadBuf *bright)
{
	const uby8 *	pbr = bright->buf;
	int		x, y;
					// may need to adjust start
	if (bright->subsample && !PequalRect(&bleft->r, &bright->r))
		pbr = bright->buf + ImgDataSize[bright->cs.dtype] *
				(bleft->r.ytop*mp->Width() + bleft->r.xleft);

	switch (bleft->cs.dtype) {
	case IDTfloat: {
		const float *	pfl = (const float *)bleft->buf;
		const float *	pfr = (const float *)pbr;
		for (y = bleft->r.ytop; y < bleft->r.ybottom; y++)
			for (x = bleft->r.xleft; x < bleft->r.xright; x++,
						pfr += bright->subsample)
				if (fcompare(*pfl++, co, *pfr))
					mp->Set(x, y);
		} break;
	case IDTushort: {
		const int		seps = int(epsilon*65536.);
		const unsigned short *	psl = (const unsigned short *)bleft->buf;
		const unsigned short *	psr = (const unsigned short *)pbr;
		for (y = bleft->r.ytop; y < bleft->r.ybottom; y++)
			for (x = bleft->r.xleft; x < bleft->r.xright; x++,
						psr += bright->subsample)
				if (icompare(*psl++, co, *psr, seps))
					mp->Set(x, y);
		} break;
	case IDTubyte: {
		const int		beps = int(epsilon*256.);
		const uby8 *	pbl = bleft->buf;
		for (y = bleft->r.ytop; y < bleft->r.ybottom; y++)
			for (x = bleft->r.xleft; x < bleft->r.xright; x++,
						pbr += bright->subsample)
				if (icompare(*pbl++, co, *pbr, beps))
					mp->Set(x, y);
		} break;
	default:
		DMESG(DMCparameter, "Unsupported data type for comparison");
		return false;
	}
	return true;
}

// See if two input images are identically encoded 16-bit depth maps
bool
matchDepthEncoding16(ImgReader *ir1, ImgReader *ir2)
{
	static const char	refDepth[] = "REFDEPTH";
	ImgInfo			iInfo;
	char			rd1[256];
	const char *		rd2;

	if ((ir1->ri != &IRInterfaceDPT) | (ir2->ri != &IRInterfaceDPT))
		return false;
	if (IRgetInfo(ir1, &iInfo) != IREnone ||
			!GetImgInfoParam(&iInfo, refDepth, rd1))
		return false;
	if (IRgetInfo(ir2, &iInfo) != IREnone ||
			!(rd2 = FindImgInfoParam(&iInfo, refDepth)))
		return false;
	rd2 += sizeof(refDepth);	// skips '=' as well
	if (!strcasecmp(rd1, rd2))
		return true;		// easy match
	double	drat = atof(rd2);
	if (drat <= 1e-7)
		return false;
	drat = atof(rd1)/drat;
	return (0.995 <= drat) & (drat <= 1.005);
}

// Perform comparison between two specifications and return as bitmap
bool
make_comparison(ABitMap2 *rmp, char *left_spec, CompareOp co, char *right_spec)
{
	ImgReadBuf	rbLeft, rbRight;
					// figure out what we're comparing
	if (isflt(left_spec)) {
		if (isflt(right_spec)) {
			DMESG(DMCinput, "Comparison between two constants");
			return false;
		}
		co = revcomp[co];	// reverse our comparison
		char *	ts = left_spec;
		left_spec = right_spec;
		right_spec = ts;
		rbRight.subsample = 0;
	} else
		rbRight.subsample = !isflt(right_spec);
					// open input image(s)
	ImgReader *	irLeft = PopenImageF(left_spec, false, NULL);
	if (!irLeft)
		return false;
	if (ImgPixelLen[irLeft->cs.format] != 1) {
		DMESGF(DMCdata, "Input '%s' not single-channel", left_spec);
		return false;
	}
	rbLeft.cs = irLeft->cs;
	rbLeft.subsample = 1;
	rbLeft.r = irLeft->fr;		// always read left by scanlines
	ImgReader *	irRight = NULL;
	PixelVal	Kpval;
	if (rbRight.subsample) {	// load right image if one
		irRight = PopenImageF(right_spec, false, NULL);
		if (!irRight)
			return false;
		if (ImgPixelLen[irRight->cs.format] != 1) {
			DMESGF(DMCdata, "Input '%s' not single-channel", right_spec);
			return false;
		}
		if ((irLeft->xres != irRight->xres) | (irLeft->yres != irRight->yres)) {
			sprintf(dmessage_buf, "'%s' and '%s' dimensions do not match",
					left_spec, right_spec);
			DMESG(DMCdata, dmessage_buf);
			return false;
		}
		if (irRight->cs.dtype != irLeft->cs.dtype) {
			sprintf(dmessage_buf, "'%s' and '%s' have different data types",
					left_spec, right_spec);
			DMESG(DMCdata, dmessage_buf);
			return false;
		}
		rbRight.cs = irRight->cs;
					// special case for 16-bit encoded depth
		if (matchDepthEncoding16(irLeft, irRight))
			rbLeft.cs.dtype = rbRight.cs.dtype = IDTushort;
		if (PequalRect(&irLeft->fr, &irRight->fr)) {
			DMESG(DMCtrace, "Comparing images by scanline");
			rbRight.r = irRight->fr;
		} else {		// load entire right image
			DMESG(DMCtrace, "Comparing full image");
			rbRight.r.xleft = rbRight.r.ytop = 0;
			rbRight.r.xright = irRight->xres;
			rbRight.r.ybottom = irRight->yres;
		}
		rbRight.buf = (uby8 *)malloc(ImgReadBufLen(&rbRight));
	} else {			// else assign constant value
		rbRight.cs = rbLeft.cs;
		Kpval.csp = &ICS_Y;
		Kpval.v.f[0] = atof(right_spec);
		if (rbRight.cs.dtype != IDTfloat &&
				(Kpval.v.f[0] < 0) | (Kpval.v.f[0] >= 1)) {
			DMESGF(DMCinput,
		"Comparison value %f not in [0,1) range for non-float image",
					Kpval.v.f[0]);
			return false;
		}			// convert constant value
		Kpval = PconvertPixel(Kpval, &rbRight.cs);
		rbRight.buf = Kpval.v.b;
		DMESG(DMCtrace, "Comparing image to constant");
		memset(&rbRight.r, 0, sizeof(ImgRect));
	}
	rbLeft.buf = (uby8 *)malloc(ImgReadBufLen(&rbLeft));
	if (!rbLeft.buf | !rbRight.buf) {
		DMESG(DMCmemory, "Cannot allocate image buffer");
		return false;
	}
	rmp->NewBitMap(irLeft->xres, irLeft->yres);
	do {
		if (IRreadRec(irLeft, &rbLeft) != IREnone) {
			PreportReaderErr(irLeft);
			return false;
		}
		bool	moreRight = irRight && IRmoreRec(irRight);
		if (moreRight && IRreadRec(irRight, &rbRight) != IREnone) {
			PreportReaderErr(irRight);
			return false;
		}
		if (!compare_buf(rmp, &rbLeft, co, &rbRight))
			return false;

		rbLeft.r = irLeft->nr;
		if (moreRight) rbRight.r = irRight->nr;

	} while (IRmoreRec(irLeft));

	free(rbLeft.buf);
	IRclose(irLeft);
	if (irRight) {
		free(rbRight.buf);
		IRclose(irRight);
	}
	return true;
}

// Operate on input bitmaps or compare two grayscale images
int
main(int argc, char *argv[])
{
	char *		outbmp = NULL;
	bool		force = false;
	bool		checkOp = false;
	int		op = 0;
	ABitMap2	curMap;
	int		a;

	fixargv0(argv[0]);				// sets progname
	dmessage_class_flags[DMCwarning] |= DMFstderr;
							// supported inputs
	PaddIReaderI(&IRInterfaceBMP);
	PaddIReaderI(&IRInterfaceTIFF);
	// PaddIReaderI(&IRInterfaceEXR);
	PaddIReaderI(&IRInterfaceDPT);
	PaddIReaderI(&IRInterfaceMTX);

	for (a = 1; a < argc; a++) {
		if (!outbmp & (argv[a][0] == '-'))	// option(s)
			switch (argv[a][1]) {
			case 'v':
				dmessage_class_flags[DMCinfo] |= DMFstderr;
				dmessage_class_flags[DMCtrace] |= DMFstderr;
				continue;
			case 'f':
				force = true;
				continue;
			case 'e':
				if (a >= argc-1) 
					usage();
				epsilon = atof(argv[++a]);
				epsilon *= float(epsilon > 0);
				continue;
			default:
				usage();
			}
		if (!outbmp) {				// need output bmp
			outbmp = argv[a++];
			if (a >= argc || strcmp(argv[a], "="))
				usage();
			if (!force && access(outbmp, 0) == 0) {
				DMESGF(DMCinput,
		"Output file '%s' already exists -- use '-f' to overwrite",
						outbmp);
				return 1;
			}
			continue;
		}
		if (checkOp) {				// expecting operator?
			switch (argv[a][0]) {
			case '^':
			case '|':
			case '+':
			case '&':
			case '-':
				if (argv[a][1])
					usage();
				op = argv[a][0];
				break;
			default:
				if (!strcasecmp(argv[a], "and"))
					op = '&';
				else if (!strcasecmp(argv[a], "or"))
					op = '|';
				else if (!strcasecmp(argv[a], "xor"))
					op = '^';
				else
					usage();
				break;
			}
			if (a+1 >= argc)		// missing final expr?
				usage();
			checkOp = false;
			continue;
		}
		ABitMap2	newMap;			// get (next) input
		CompareOp	co;			// check for comparison
		if (a+2 < argc && (co = iscomparison(argv[a+1]))) {
			if (!make_comparison(&newMap, argv[a], co, argv[a+2]))
				return 1;
			a += 2;
		} else {				// load 2-D bitmap
			bool	uNot = (argv[a][0] == '~');
			if (!ReadBitMap2(&newMap, argv[a]+uNot))
				return 1;
			if (uNot) {
				DMESG(DMCtrace, "INVERT operation");
				newMap.Invert();
			}
		}
		if (op && (newMap.Width() != curMap.Width()) |
				(newMap.Height() != curMap.Height())) {
			DMESG(DMCdata, "Input bitmap dimensions do not match");
			return 1;
		}
		switch (op) {
		case 0:				// first input
			curMap.Take(&newMap);
			break;
		case '^':
			DMESG(DMCtrace, "XOR operation");
			curMap ^= newMap;
			break;
		case '|':
		case '+':
			DMESG(DMCtrace, "OR operation");
			curMap |= newMap;
			break;
		case '&':
			DMESG(DMCtrace, "AND operation");
			curMap &= newMap;
			break;
		case '-':
			DMESG(DMCtrace, "SUBTRACT operation");
			curMap -= newMap;
			break;
		default:
			DMESG(DMCparameter, "Botched operation!");
			return 1;
		}
		checkOp = true;				// on to next (if any)
	}
	if (!outbmp | !curMap.Width())
		usage();
	DMESGF(DMCinfo, "%.4f%% bits set in final map",
			curMap.SumTotal()*100./curMap.Width()/curMap.Height());
	return !WriteBitMap2(curMap, outbmp);
}
