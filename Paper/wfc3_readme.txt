+------------------------------------------------------------------------------+

Mosaic images of M83 from the Hubble Space Telescope 
WFC3 Early Release Science (ERS) program

  9 February 2011

  Max Mutchler 
  Research and Instrument Scientist
  Space Telescope Science Institute
  3700 San Martin Drive, Baltimore, MD 21218

Table of Contents:

     1. Introduction and observations
     2. Data calibration and registration
     3. Rejection of cosmic rays and other artifacts
     4. Image combination
     5. Data products and filenaming convention
     6. Publications
     7. References
     
+------------------------------------------------------------------------------+

1. Introduction and observations

The new Wide Field Camera 3 (WFC3) was installed by astronauts during Hubble
Servicing Mission 4 (SM4) in May 2009. In August 2009, mosaic images of M83 were
obtained using both the visible (UVIS) and infrared (IR) channels of WFC3, as
part of the Early Release Science (ERS) program conducted by the WFC3 Science
Oversight Committee (SOC):

  Bruce Balick, Howard E. Bond, Daniela Calzetti, C. Marcella Carollo, 
  Michael J. Disney, Michael A. Dopita, Jay A. Frogel, Donald N.B. Hall, 
  Jon A. Holtzman, Randy A. Kimble, Patrick J. McCarthy, Robert W. O'Connell 
  (chair), Francesco Paresce, Abhijit Saha, Joseph I. Silk, John T. Trauger, 
  Alistair R. Walker, Bradley C. Whitmore, Rogier A. Windhorst, Erick T. Young

If you utilize these prepared images or catalogs for scientific or educational
purposes, please acknowledge the WFC3 Science Oversight Committee.

The new WFC3 observations include two adjacent pointings which form a small
mosaic in the UVIS channel, but which are discontinuos in the IR channel. Small
dithers within each pointing were used to fill the gaps between the two UVIS
CCD chips, and allow for the rejection of cosmic rays and detector artifacts. A
mosaic image was also obtained for a parallel field, using the Advanced Camera
for Surveys (ACS, which was repaired during SM4). In summary, new images were
obtained using the following detectors and filters:

  detector   filter  
  ---------  ------

  WFC3/UVIS  F225W  Wide UV
  WFC3/UVIS  F336W  U-band
  WFC3/UVIS  F373N  [O II]
  WFC3/UVIS  F438W  B-band
  WFC3/UVIS  F487N  H-beta
  WFC3/UVIS  F502N  [O III]
  WFC3/UVIS  F555W  V-band, South field
  WFC3/UVIS  F547M  V-band, North field
  WFC3/UVIS  F657N  Wide H-alpha and [N II]
  WFC3/UVIS  F673N  [S II]
  WFC3/UVIS  F814W  I-band

  WFC3/IR    F110W  Wide YJ
  WFC3/IR    F128N  Paschen-beta 
  WFC3/IR    F160W  H-band
  WFC3/IR    F164N  [Fe II]

  ACS/WFC    F555W  V-band parallel field
  ACS/WFC    F625W  R-band parallel field (SDSS r)

For more observation details, see visits 61-63 and B1-B3 in 
HST program 11360 (PI Bob O'Connell):

http://www.stsci.edu/observing/phase2-public/11360.pro

Since our carefully prepared mosaics represent a significant investment of
expert processing beyond the standard archival products, we are releasing this 
drizzle-combined FITS data as a High-Level Science Product via the Multimission
Archive at STScI (MAST):

http://archive.stsci.edu/prepds/wfc3ers

The resulting color composite image was released on 5 November 2009, 
and more general information can be found via the press release:

http://hubblesite.org/newscenter/archive/releases/2009/29/

+------------------------------------------------------------------------------+

2. Data calibration and registration

The archival raw data was retrieved from MAST, and although it is automatically
calibrated on-the-fly, the latest available calibration reference files were 
only available from the WFC3 website, and were applied by running the
calibration pipeline (calwf3) offline. As of this writing, good inflight
WFC3/UVIS bias and dark calibrations are automatically applied, but no flat
fielding. The calibration reference files also propagate data quality flags
which identify many types of detector artifacts, most of which were excluded
during the MultiDrizzle processing described below. These include bad CCD
columns, hot pixels, and saturated pixels. As input to the drizzle-combination
described below, we generated and used our own flt images. 

The first and last four rows of each CCD where manually masked for rejection by
setting the data quality value for all those pixels to 8192 in the [DQ]
extentions of the input flt images.

Masking for the IR channel "blobs" (probably dust motes) is now in place,
but since our dithers were not large enough to span them, these pixels
were not excluded in our final output.

Mosaic datasets involving shifts larger than about 100 arcseconds necessarily
involve different guide stars in different visits, and so guide star positional
errors render the WCS unreliable for adequate image registration. Therefore, we
used stars (or small clusters) in the images to measure relative positions and
apply delta-shifts to register the frames. We generally achieved internal image
registrations to within 0.1 pixels. No attempt to align the images with any
other astrometric catalog was made, and so there is likely a small offset from
catalog positions, on the order of ~0.5 arcsec. The following is a sample shift
file for UVIS F814W, showing columns for x-shift, y-shift, rotation and scale
corrections:

  # frame: output
  # refimage: ib6w62tvq_wcs.fits
  # form: delta 
  # units: pixels 

  ib6w62tvq_flt.fits     0.00    0.00    0.0    1.0
  ib6w62txq_flt.fits     0.00   -0.16    0.0    1.0
  ib6w62tzq_flt.fits    -0.16    0.00    0.0    1.0

  ib6wb2b6q_flt.fits    -6.74  -15.69    0.0    1.0
  ib6wb2b8q_flt.fits    -6.95  -15.86    0.0    1.0
  ib6wb2bzq_flt.fits    -6.72  -16.00    0.0    1.0

+------------------------------------------------------------------------------+

3. Rejection of cosmic rays and other artifacts

Although the rejection of cosmic rays and most detector artifacts is embedded in
the MultiDrizzle processing described below, the following is an overview of how
it is accomplished for the WFC3 UVIS data.

A median image is constructed from undistorted single-drizzled
(*single_sci.fits) images. This median image -- or the appropriate sections of
it -- are blotted back to the distorted space of the input images, where it can
be used to identify cosmic rays. The dither package tasks of "deriv" and
"driz_cr" are used to compare this blotted image and its derivative image with
the original input image, and generate a cosmic ray mask. Finally, all the
input images, together with their newly created cosmic ray masks, are drizzled
onto a single output mosaic, which has units of electrons/second in each pixel.

MultiDrizzle also produces exposure weight maps, which indicate the background
and instrumental noise for each pixel in the science data. Because there is some
overlap between adjacent tiles in the mosaic, and all the bad pixels which are
rejected, the total exposure time in the final mosaic varies as a function of
position. The weight maps were visually inspected to ensure that they showed the
pixels affected by cosmic rays and detector artifacts to have lower weight. It
is equally important to verify that real objects (e.g. the cores of stars, or
bright nebular features) are not being rejected. The only irregular rejection
was related to a partially-rejected window ghost (caused by the bright nucleus
of M83) in the overlap area of the F438W and F814W mosaics (the extent of it is
evident in the weight maps).

For the WFC3 IR data, the flt images are already free of cosmic rays, due to the
up-the-ramp readouts and pipeline processing. So the rejection step (drz_cr) was
turned off.

Some artifacts remain in the final images. In the UVIS mosaics, there are some
residual cosmic rays in the gap overlap strips. In the IR mosaics, some bad
detector pixels are evident near the edges of the image -- this includes the
large elliptical "death star" artifact. Also, the IR "blobs" are not masked out.

+------------------------------------------------------------------------------+

4. Image combination

The IRAF/STSDAS MultiDrizzle task (Koekemoer et al., 2002) was used to combine
the mosaics. MultiDrizzle is a PyRAF script which performs, on a list of input
images: bad pixel identification, sky subtraction*, rejection of cosmic rays and
other artifacts (as described above), and a final drizzle combination (Fruchter
et al., 2002) with the cosmic ray masks, into a clean image. MultiDrizzle also
applies the latest filter-specific geometric distortion corrections to each
image, as specified in the IDCTAB reference tables. The output images are
oriented with North up and East left. For more on MultiDrizzle see:

http://www.stsci.edu/hst/wfc3/tools/MultiDrizzle/              (WFC3 drizzling)
http://www.stsci.edu/hst/HST_overview/documents/multidrizzle/  (Handbook)

The following key MultiDrizzle parameters were used for the WFC3 UVIS data:

multidrizzle.mdriztab = no
multidrizzle.coeffs = 'header'
multidrizzle.context = no
multidrizzle.clean = no
multidrizzle.ra = '204.26884'
multidrizzle.dec = '-29.83995' 
multidrizzle.build = no
multidrizzle.shiftfile = 'shifts_f814w.txt' 
multidrizzle.static = yes  
multidrizzle.static_sig = 5.0
multidrizzle.skysub = no
multidrizzle.driz_separate = yes
multidrizzle.driz_sep_outnx = 5000
multidrizzle.driz_sep_outny = 8500
multidrizzle.driz_sep_kernel = 'square'
multidrizzle.driz_sep_wt_scl = 'exptime'
multidrizzle.driz_sep_scale = 0.0396
multidrizzle.driz_sep_pixfrac = 1.0
multidrizzle.driz_sep_rot = 0.0
multidrizzle.driz_sep_fillval = 9999.9
multidrizzle.driz_sep_bits = 4352
multidrizzle.median = yes
multidrizzle.combine_type = 'minmed'
multidrizzle.combine_nlow = 0
multidrizzle.combine_nhigh = 1
multidrizzle.combine_hthresh = 8888.8  # exclude fill values
multidrizzle.blot = yes
multidrizzle.driz_cr = yes
multidrizzle.driz_cr_snr = '4.0 3.0' 
multidrizzle.driz_cr_scale = '1.2 0.7'  
multidrizzle.driz_combine = yes
multidrizzle.final_wht_type = 'EXP'
multidrizzle.final_outnx = 5000
multidrizzle.final_outny = 8500
multidrizzle.final_kernel = 'square'
multidrizzle.final_wt_scl = 'exptime'
multidrizzle.final_scale = 0.0396
multidrizzle.final_pixfrac = 1.0
multidrizzle.final_rot = 0.0
multidrizzle.final_fillval = 0.0
multidrizzle.final_bits = 4352
multidrizzle.crbit = 0

The following key MultiDrizzle and imcombine parameters were used for 
the WFC3 IR data:

multidrizzle.mdriztab = no
multidrizzle.workinplace = no
multidrizzle.coeffs = 'header'
multidrizzle.context = no  
multidrizzle.clean = no    
multidrizzle.ra = '204.26855'
multidrizzle.dec = '-29.84029'   
multidrizzle.build = no   
multidrizzle.shiftfile = 'shifts_f160w.txt' 
multidrizzle.static = no
multidrizzle.skysub = no  
multidrizzle.driz_separate = yes
multidrizzle.driz_sep_outnx = 2000  
multidrizzle.driz_sep_outny = 3400  
multidrizzle.driz_sep_kernel = 'square'   # better than turbo
multidrizzle.driz_sep_wt_scl = 'exptime'
multidrizzle.driz_sep_scale = 0.1283
multidrizzle.driz_sep_pixfrac = 1.0
multidrizzle.driz_sep_rot = 0.0   
multidrizzle.driz_sep_fillval = 9999.9   #  arbitrary high value
multidrizzle.driz_sep_bits = 4864  #  include saturation and blobs
multidrizzle.median = yes
multidrizzle.combine_type = 'median'   
multidrizzle.combine_nlow = 0
multidrizzle.combine_nhigh = 0    # zero helps suppress motes, persistence
multidrizzle.combine_hthresh = 8888.8  # exclude fill values
multidrizzle.blot = yes
multidrizzle.driz_cr = no
multidrizzle.driz_combine = yes
multidrizzle.final_wht_type = 'EXP'
multidrizzle.final_outnx = 2000
multidrizzle.final_outny = 3400
multidrizzle.final_kernel = 'gaussian'  
multidrizzle.final_wt_scl = 'exptime'
multidrizzle.final_scale = 0.0792   # enhanced scale, 2X UVIS scale
multidrizzle.final_pixfrac = 0.8
multidrizzle.final_rot = 0.0  # rotate North up, East left
multidrizzle.final_fillval = 0.0
multidrizzle.final_bits = 4864 
multidrizzle.crbit = 0  

* Note that although MultiDrizzle substracts sky in order to reject cosmic rays,
that sky value was not subtracted from our final drizzled mosaics.

A detailed description of drizzling a WFC3 mosaic dataset can be found in the
2010 HST Calibration Workshop (Mutchler, 2010), which was also re-published as
section 4.3 of the WFC3 Data Handbook (Rajan, 2010).

+------------------------------------------------------------------------------+

5. Data products and filenaming convention

The drizzled mosaics described here are available as High-Level Science Products
(HLSP) via the Multimission Archive at STScI (MAST). The filenames and
descriptions of our drizzled science mosaics are listed below. This filenaming
convention makes these HLSP files compatible and queryable with MAST and the
Hubble Legacy Archive (HLA). Note that we include the early (version 1 or "v1"
in filenames) images used to produce some of the early publications listed
below, and these files contain only the south field. But we later combined both
fields using the best calibrations available at the time of this HLSP data
release. We recommend the use of these version 2 ("v2" in filenames) files for
most subsequent analyses.

The WFC3 UVIS mosaics have the native detector scale of 0.0396 arcsec/pixel. The
enhanced-resolution IR mosaics have a scasle of 0.0792 arcsec/pixel,  which is
62% of the IR detector scale (and exactly twice the UVIS scale). For each
drizzled science mosaic (drz_sci), there is a corresponding exposure weight map
(drz_wht), and also small preview JPEG images available on the MAST website. 

WFC3 UVIS mosaics (native scale 0.0396 arcsec/pixel):

  hlsp_wfc3ers_hst_wfc3_m83_f225w_v2_drz_sci.fits      UV wide
  hlsp_wfc3ers_hst_wfc3_m83_f336w_v2_drz_sci.fits      U-band (Stromgren u)
  hlsp_wfc3ers_hst_wfc3_m83_f373n_v2_drz_sci.fits      [O II] 
  hlsp_wfc3ers_hst_wfc3_m83_f438w_v2_drz_sci.fits      B-band (WFPC2 B)
  hlsp_wfc3ers_hst_wfc3_m83_f487n_v2_drz_sci.fits      H-beta
  hlsp_wfc3ers_hst_wfc3_m83_f502n_v2_drz_sci.fits      [O III] 
  hlsp_wfc3ers_hst_wfc3_m83-n_f547m_v2n_drz_sci.fits   V-band north(Stromgren y)
  hlsp_wfc3ers_hst_wfc3_m83-s_f555w_v2s_drz_sci.fits   V-band south (WFPC2 V)
  hlsp_wfc3ers_hst_wfc3_m83_f673n_v2_drz_sci.fits      Wide H-alpha and [N II]
  hlsp_wfc3ers_hst_wfc3_m83_f673n_v2_drz_sci.fits      [S II] 
  hlsp_wfc3ers_hst_wfc3_m83_f814w_v2_drz_sci.fits      I-band (WFPC2 I)

WFC3 UVIS photometric catalogs (derived from version 1 images):

  hlsp_wfc3ers_hst_wfc3_m83_whitelight_v1.fits        UBVI detection image
  hlsp_wfc3ers_hst_wfc3_m83_cat_all_v1.txt            Catalog of all sources
  hlsp_wfc3ers_hst_wfc3_m83_cat_cluster_auto_v1.txt   Automatic cluster catalog
  hlsp_wfc3ers_hst_wfc3_m83_cat_cluster_manual_v1.txt Manual cluster catalog

WFC3 IR mosaics (discontinuous, enhanced scale 0.0792 arcsec/pixel):

  hlsp_wfc3ers_hst_wfc3_m83-x_f110w_v2x_drz_sci.fits   Wide YJ
  hlsp_wfc3ers_hst_wfc3_m83-x_f128n_v2x_drz_sci.fits   Paschen-beta
  hlsp_wfc3ers_hst_wfc3_m83-x_f160w_v2x_drz_sci.fits   H-band
  hlsp_wfc3ers_hst_wfc3_m83-x_f164n_v2x_drz_sci.fits   [Fe II]

The ACS parallel field images (F555W and F625W) are not included in this HLSP
collection, but pipeline-drizzled images are available from MAST and the Hubble
Legacy Archive (HLA):

http://hla.stsci.edu

+------------------------------------------------------------------------------+

6. SOC papers resulting from this data (as of Feb 2011):

"The Luminosity, Mass, and Age Distribution of Compact Star Clusters in M83 
based on HST/WFC3 Observations", R. Chandar et al., 2010, ApJ, 719, 966 

"Supernova Remnants and the Interstellar Medium of M83: Imaging and Photometry 
with the WFC3 on HST", M. Dopita et al., 2010, ApJ 710, 964

"Large-scale Shock-ionized and Photo-ionized Gas in M83: The Impact of 
Star Formation," S. Hong et al., submitted (Dec, 2010)

"Resolved Stars in M83 Based on HST/WFC3 ERS observations", 
H. Kim et al. (in preparation)

"Using Halpha Morphology and Surface Brightness Fluctuations to Age-Date 
Star Clusters in M83", B. Whitmore, ApJ, 729, 1

See a complete list of SOC publications on the HLSP webpage:
http://archive.stsci.edu/prepds/wfc3ers

+------------------------------------------------------------------------------+

7. References

Fruchter, A. S., Hook, R. N., 2002, PASP 114, 144
http://arxiv.org/abs/astro-ph/9808087

Koekemoer, A. M., Fruchter, A. S., Hook, R. N., Hack, W. 2002,
HST Calibration Workshop, Ed. S. Arribas, A. M. Koekemoer,
B. Whitmore (STScI: Baltimore), p.337
http://www.stsci.edu/hst/HST_overview/documents/calworkshop/workshop2002

Mutchler, M., 2010, HST Calibration Workshop (in press)
http://www.stsci.edu/institute/conference/cal10/proceedings

Rajan, A., et al., 2010, WFC3 Data Handbook, version 2.0 (Baltimore: STScI)

+------------------------------------------------------------------------------+
