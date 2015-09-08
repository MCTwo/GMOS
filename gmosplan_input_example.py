'''
This file provides an example of how to interact and run obsplan.py

All position angles are defined +CCW from North towards East
'''
from __future__ import division
import numpy
import tools
import obsplan
import gmosplan

###########################################################################
## USER INPUT
###########################################################################

## General inputs

#Prefix for all output files
prefix = 'zwcl1856_mask1_revA'
# Hour angle of the target for the mask (float; unit:hours)
HA = 1.0

## Star and Galaxy catalog inputs

catalog = '/Users/dawson/OneDrive/Research/Clusters/ZwCL1856/catalogs/Gemini/unmosaic/zwcl1856_join_r_ind.txt'
objid_ttype = 'NUMBER'
ra_ttype = 'ALPHA_J2000'
dec_ttype = 'DELTA_J2000'
mag_ttype = 'MAG_AUTO'
xccd_ttype = 'X_IMAGE'
yccd_ttype = 'Y_IMAGE'
passband = 'r'
equinox = '2000'
slitwidth = 1.0 # arcsec
pixscale = 0.146 # pixel scale of the image

# Go ahead and read in the catalog since this will be needed to create the
# galaxy mask later in the user input section
cat = tools.readcatalog(catalog)
key = tools.readheader(catalog)

# read in the separate g-band catalog
catalog_g = '/Users/dawson/OneDrive/Research/Clusters/ZwCL1856/catalogs/Gemini/unmosaic/zwcl1856_join_g_ind.txt'
cat_gband = tools.readcatalog(catalog_g)
key_gband = tools.readheader(catalog_g)

# Clean up the catalogs, getting rid of junk/noisy objects
mask_temp1 = numpy.logical_and(cat[:,key['MAGERR_AUTO']] > 0.00001,cat_gband[:,key_gband['MAGERR_AUTO']] > 0.00001)
mask_temp2 = numpy.logical_and(cat[:,key['MAGERR_AUTO']] < 1.0, cat_gband[:,key_gband['MAGERR_AUTO']] < 1.0)
mask_magERR = numpy.logical_and(mask_temp1, mask_temp2)
mask_mag99 = numpy.logical_and(cat[:,key['MAG_AUTO']] != 99, cat_gband[:,key_gband['MAG_AUTO']] != 99)
mask_clean = mask_mag99*mask_magERR
cat = cat[mask_clean,:]
cat_gband = cat_gband[mask_clean,:]

## Slitmask ds9 region input

# Path/Name of the ds9 region file defining the bound and orientation of the
# slit mask. Note region should be defined using Coordinate/WCS/Degrees and
# Size/WCS/Arcmin options, with the Size 5 by 16.1 arcmin, Angle will then
# correspond to the slitmask's parallactic angle (i.e. +CCW from north towards
# east) with the guider camera in the North-east quadrent at Angle=0.
regfile = 'mask_north.reg'

## Slit size inputs

#The amount of sky on either side of the galaxy to include in slit (arcsec)
sky = (0.75,0.75)
# The ttype index of the galaxy size. If one value is entered then the galaxy
# will be assumed circular, of three values are entered the the galaxy will be
# assumed elliptical with major axis radius (a), and minor axis radius (b), and
# position angle (pa_gal) of the major axis measured +CCW from North towards
# east
A_gal_ttype = 'A_IMAGE'
B_gal_ttype = 'B_IMAGE' # if None then galaxy assumed circular
pa_ga_ttype = 'THETA_IMAGE' # if None then galaxy assumed circular

## Guide and alignment star inputs

# Enter lists of the object ids for each type of star, these will be matched to
# the object ids in the catalog

# # Guide star id's
# gs_ids = (32013)
# # Alignment star id's
# as_ids = (29400,31060,31086,34795,39587,39822)

## Exclusion list input

# ttype catalog of galaxies to exclude from mask (excludes matching ttype = objid). exobjid_ttype is ['string'] ttype name of the objid column in the exfile, the objid's should correspond to some of the objid's in the objid array

exfile = None
exobjid_ttype = 'NUMBER'
#exfile = '/Users/dawson/SkyDrive/Observing/Keck2013a/MACSJ1752/macs1752_Mask1_rev0_maskcat.txt' #a string (e.g. 'exclusion.txt') or None
#exobjid_ttype = 'objID'

## Priority code (i.e. selection weight) input

# Currently obsplan is only setup to calculate an objects priority code based on
# its photo-z relative to the cluster redshift
#z_cluster = 0.19
#photo_z_ttype = 'z_phot'
#photo_z_err_ttype = 'z_phot_Err'

## Create a sample definition

# Currently obsplan is only setup to break samples according to one object
# variable (e.g. magnitude).  The sampel_param_ttype is the ttype name of the
# vairable in the catalog to be used to make the sample division (e.g.
# magnitude). samplebounds defines the min and max of sample_param for each
# sample: e.g. (sample1 lowerbound, sample1 upperbound, sample2 lower bound,
# sample2 upper bound, etc., etc.).
samplecut = ('MAG_AUTO',20.0)
## Define hard catalog cuts/limits
R_bounds = (12,22.5)

## Preselected list input

# ttype catalog of galaxies to preselect in dsim.
psfile = None #a string (e.g. 'preselect.txt') or None
psobjid_ttype = None

## Create galaxy selection mask

# Some how the user at this point needs to create a boolean type mask for the
# catalog that will filter out all objects that are not galaxies
mask_galaxy = numpy.logical_and(cat_gband[:,key['CLASS_STAR']]<0.985,cat_gband[:,key['FWHM_IMAGE']]>4)

## Apply zero point corrections
cat_gband[:,key_gband['MAG_AUTO']] += 0
#cat_gband[:,key_gband['MAG_APER']] += 0
cat[:,key['MAG_AUTO']] += 0
#cat[:,key['MAG_APER']] += 0

## Create a magnitude mask

# It is likely that the a brigth and faint end mask should be used
mask_mag = numpy.logical_and(cat[:,key[mag_ttype]] >= R_bounds[0],
                             cat[:,key[mag_ttype]] <= R_bounds[1])

###########################################################################
## Automated Portion
###########################################################################

# create basic 1D arrays from catalog
objid = cat[:,key[objid_ttype]]
ra = cat[:,key[ra_ttype]]
dec = cat[:,key[dec_ttype]]
mag = cat[:,key[mag_ttype]]
x_ccd = cat[:,key[xccd_ttype]]
y_ccd = cat[:,key[yccd_ttype]]

# Create the slitmask mask
mask_slitmask = obsplan.createSlitmaskMask(regfile,ra,dec)

# Create the exclusion list mask
if exfile == None:
    # then make a numpy.array of just True's
    mask_ex = numpy.ones(numpy.size(objid)) == 1
elif exfile != None and exobjid_ttype != None:
    mask_ex = obsplan.createExclusionMask(objid,exfile,exobjid_ttype)

####################################################
## Determine the "sample" assignment for each object
####################################################
# Define red sequence bounds
# (x1,y1,x2,y2) coordinates that define the line
rs_high = (12, 2.25, 21, 1.75)
rs_low = (12, 1.25, 21, 0.75)

# Select cluster members based on a red sequence cut
# add the assumed redsequence bounds
x = numpy.arange(R_bounds[0],R_bounds[1])
m_high = (rs_high[3]-rs_high[1])/(rs_high[2]-rs_high[0])
m_low = (rs_low[3]-rs_low[1])/(rs_low[2]-rs_low[0])
y_high = m_high*(x-rs_high[0])+rs_high[1]
y_low = m_low*(x-rs_low[0])+rs_low[1]

# Create a red sequence selection mask
v_g = cat_gband[:,key_gband['MAG_AUTO']]
r_g = cat[:,key['MAG_AUTO']]
mask_rs = numpy.logical_and(v_g-r_g>=m_low*(cat[:,key[mag_ttype]]-rs_low[0])+rs_low[1], v_g-r_g<=m_high*(cat[:,key[mag_ttype]]-rs_high[0])+rs_high[1])

# create a blue cloud selection mask
mask_bc = v_g-r_g<m_low*(cat[:,key[mag_ttype]]-rs_low[0])+rs_low[1]

# Determine the sample assignments
#0) alignment stars
#1) red sequence - red sequence, with R <= 21.5
#2) blue cloud - red sequence, with 21.5 < R <=22.5
#3) everything else - blue cloud, with R <= 21.5
sample = 1*mask_rs + 2*mask_bc + 3*~mask_rs*~mask_bc
sample *= mask_galaxy

# Determine the selection flag for each galaxy. If non-zero then the object is
# preselected
selectflag = obsplan.assignSelectionFlag(objid,psfile,psobjid_ttype)

# determine object declination and the mask PA from the regfile
box = obsplan.readMaskRegion(regfile)
delta = box[1]
pa_mask = box[4]

# Determine the optimal slit PA
pa_slit = obsplan.optimalPA(pa_mask,HA,delta)

# Determine the slit size for each object
A_gal = numpy.sqrt(cat[:,key[A_gal_ttype]])*0.2 # suprimecam pixel size to arcsec size
if B_gal_ttype == None:
    B_gal = None
else:
    B_gal = cat[:,key[B_gal_ttype]]*0.2
if pa_ga_ttype == None:
    pa_gal = None
else:
    pa_gal = cat[:,key[pa_ga_ttype]]
len1, len2 = obsplan.slitsize(pa_slit,sky,A_gal,B_gal,pa_gal)

# Filter the galaxy catalog before creating dsim input
mask_temp = mask_mag*mask_ex*mask_slitmask
# but need to include any preselected galaxies that might be excluded by the
# above masks
mask = numpy.logical_or(mask_temp,selectflag)

# Create the output fits file
gmosplan.create_obstab(prefix, pixscale, objid[mask].astype(int), ra[mask],
                       dec[mask], x_ccd[mask], y_ccd[mask], mag[mask],
                       sample[mask], slitwidth, slitsize_y=len1[mask]+len2[mask],
                       slittilt=pa_slit, slitpos_y='Default')

# Create the target galaxy slit region file
length = len1+len2
obsplan.makeSlitmaskRegion(prefix,ra[mask],dec[mask],pa_slit,length[mask],sample[mask])
