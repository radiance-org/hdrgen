# hdrgen repository
Although this repository builds only a few tools, the panlib library that underpins
these tools is extensive.  Its main focus is on high dynamic-range (HDR) image manipulation,
but most operations support 16-bit/channel and 8-bit/channel as well.  Much of the library
is tailored to support the Photosphere image cataloger and HDR builder app, which is not
available for 64-bit targets and needs updating.  Volunteers on getting it working again are encouraged to step up.

These are the tools we currently build:

# hdrgen 
- takes a set of low dynamic-range (LDR) images in a common 8-bit/channel format and merges them into an HDR format

# hdrcvt 
- translates between various HDR and LDR formats and optionally performs numerous operations during conversion

# PQconvert 
- batch-converts a set of HDR images into 16-bit/channel TIFFs appropriate for input to HDR video encoders such as ffmpeg

# expose2range 
- computes and prints optimal multiplier to get input HDR image pixels within the given value bounds

# warpimage 
- warps an HDR or LDR image according to the specified input grid vertex locations

# bitmapop 
- performs a variety of comparisons and so forth between images or bitmaps to compute an output bitmap as a bilevel BMP

Most tools will give a usage message if executed without arguments.
If you have specific questions, please write to Greg Ward at <GJWard@lbl.gov> or post to the appropriate Discourse forum,
<a href="https://discourse.radiance-online.org/c/hdri/5">here</a>.
