from __future__ import division
import numpy
import pylab
import tools
import sys

def write_dsim_header(F,regfile,prefix):
    '''
    F = an opened file (e.g. F=open(filename,'w')
    box =
    '''
    box = readMaskRegion(regfile)
    F.write('#This catalog was created by obsplan.py and is intended to be used \n')
    F.write('#as input to the deimos slitmask software following the format \n')
    F.write('#outlined at http://www.ucolick.org/~phillips/deimos_ref/masks.html\n')
    F.write('#Note that the automatic generation of this file does not include\n')
    F.write('#guide or alignment stars.\n')
    F.write('#ttype1 = objID\n')
    F.write('#ttype2 = ra\n')
    F.write('#ttype3 = dec\n')
    F.write('#ttype4 = equinox\n')
    F.write('#ttype5 = magnitude\n')
    F.write('#ttype6 = passband\n')
    F.write('#ttype7 = priority_code\n')
    F.write('#ttype8 = sample\n')
    F.write('#ttype9 = selectflag\n')
    F.write('#ttype10 = pa_slit\n')
    F.write('#ttype11 = len1\n')
    F.write('#ttype12 = len2\n')
    #Write in the Slitmask information line
    F.write('{0}\t{1:0.6f}\t{2:0.6f}\t2000\tPA={3:0.2f}\n'
            .format(prefix,box[0]/15.,box[1],box[4]))

def write_guide_stars(F,gs_ids,objid,ra,dec,magnitude,equinox='2000',passband='R'):
    '''
    F = dsim file
    gs_ids = (list of integers)

    '''
    N = numpy.size(gs_ids)
    # Possibly reformat radii array to enable single float input
    if N == 1:
        gs_ids = numpy.reshape(gs_ids,(1,))
    for i in gs_ids:
        mask_s = objid == i
        if numpy.sum(mask_s) ==0:
            print 'obsplan.write_guide_stars: no objects in catalog match the guide star id:{0}, please check your input, exiting'.format(i)
            sys.exit()
        ra_i = tools.deg2ra(ra[mask_s],':')
        dec_i = tools.deg2dec(dec[mask_s],':')
        mag_i = magnitude[mask_s][0]
        F.write('{0}  {1}  {2}  {3}  {4:0.2f}  {5}  -1  0  1\n'
                .format(i,ra_i,dec_i,equinox,mag_i,passband))

def write_align_stars(F,as_ids,objid,ra,dec,magnitude,equinox='2000',passband='R'):
    N = numpy.size(as_ids)
    # Possibly reformat radii array to enable single float input
    if N == 1:
        as_ids = numpy.reshape(as_ids,(1,))
    for i in as_ids:
        mask_s = objid == i
        if numpy.sum(mask_s) ==0:
            print 'obsplan.write_align_stars: no objects in catalog match the alignment star id:{0}, please check your input, exiting'.format(i)
            sys.exit()
        ra_i = tools.deg2ra(ra[mask_s],':')
        dec_i = tools.deg2dec(dec[mask_s],':')
        mag_i = magnitude[mask_s][0]
        F.write('{0}  {1}  {2}  {3}  {4:0.2f}  {5}  -2  0  1\n'
                .format(i,ra_i,dec_i,equinox,mag_i,passband))

def write_galaxies_to_dsim(F,objid,ra,dec,magnitude,priority_code,sample,selectflag,pa_slit,len1,len2,equinox='2000',passband='R'):
    '''

    '''
    from math import floor
    for i in numpy.arange(numpy.size(objid)):
        #convert deg RA to sexadec RA
        ra_i = ra[i]/15.0
        rah = floor(ra_i)
        res = (ra_i-rah)*60
        ram = floor(res)
        ras = (res-ram)*60.
        #convert deg dec to sexadec dec
        dec_i = dec[i]
        if dec_i<0:
            sign = -1.
            dec_i = abs(dec_i)
        else:
            sign = 1.
        decd = floor(dec_i)
        res = (dec_i-decd)*60.
        decm = floor(res)
        decs = (res-decm)*60.
        if numpy.size(pa_slit) == 1:
            pa_slit_i = pa_slit
        else:
            pa_slit_i = pa_slit[i]
        if sign==-1:
            F.write('{0:0d}\t{1:02.0f}:{2:02.0f}:{3:06.3f}\t-{4:02.0f}:{5:02.0f}:{6:06.3f}\t{7}\t{8:0.2f}\t{9}\t{10:0.0f}\t{11:0.0f}\t{12:0.0f}\t{13:0.2f}\t{14:0.1f}\t{15:0.1f}\n'
                    .format(objid[i],rah,ram,ras,decd,decm,decs,equinox,magnitude[i],passband,priority_code[i],sample[i],selectflag[i],pa_slit_i,len1[i],len2[i]))
        else:
            F.write('{0:0d}\t{1:02.0f}:{2:02.0f}:{3:06.3f}\t{4:02.0f}:{5:02.0f}:{6:06.3f}\t{7}\t{8:0.2f}\t{9}\t{10:0.0f}\t{11:0.0f}\t{12:0.0f}\t{13:0.2f}\t{14:0.1f}\t{15:0.1f}\n'
                    .format(objid[i],rah,ram,ras,decd,decm,decs,equinox,magnitude[i],passband,priority_code[i],sample[i],selectflag[i],pa_slit_i,len1[i],len2[i]))

def makeSlitmaskRegion(prefix,ra,dec,pa_slit,length,sample,width=1):
    '''
    create a region file that maps the suggested slit of each galaxy
    '''
    outputname = prefix+'_slits.reg'
    out = open(outputname,'w')
    out.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+'\n')
    out.write('fk5'+'\n')
    for i in numpy.arange(numpy.size(ra)):
        ra_i = ra[i]
        dec_i = dec[i]
        length_i = length[i]
        if numpy.size(pa_slit) == 1:
            pa_slit_i = pa_slit
        else:
            pa_slit_i = pa_slit[i]
        if sample[i] == 1:
            color = 'green'
        else:
            color = 'blue'
        out.write('box({0:1.5f},{1:1.5f},{2:1.1f}",{3:1.1f}",{4:0.0f}) # color={5}'.format(ra_i,dec_i,width,length_i,pa_slit_i,color)+'\n')
    out.close()

def plotcoverage(redshift,lambda_central,filename=None):
    '''
    Creates a plot of common spectral features in the observed frame. Along with
    a spectra coverage region for the DEIMOS 1200 line grating for the given
    central wavelength.
    INPUT:
    redshift = [float] redshift of the cluster
    lambda_central = [float] central wavelength for the DEIMOS 1200 line grating
    filename = [None or string] if not None then an image will be saved with
       filename
    OUTPUT:
    -- displays a figure
    -- if filename given will save figure as a png file
    '''
    fig = pylab.figure(figsize=(20,4.5))
    pylab.xlim((lambda_central-2300, lambda_central+2300))
    xl = pylab.xlim()
    yl = (0,1)

    lcl_low = lambda_central-1300
    lcl_up = lambda_central+1300

    #Plot the wavelength coverage
    pylab.plot((lambda_central,lambda_central),yl,'-b',linewidth=3)
    pylab.plot((lcl_low,lcl_low),yl,
               '-b',alpha=0.5,linewidth=3)
    pylab.plot((lcl_up,lcl_up),yl,
               '-b',alpha=0.5,linewidth=3)
    pylab.fill_between(pylab.arange(lcl_low,lcl_up+100,100),yl[0],yl[1],
                       facecolor='blue',alpha=0.25)

    #At the ends of the coverage window show the +/- range due to where the
    #slit is placed on the mask
    pylab.plot((lcl_low+411,lcl_low+411),yl,
               '--b',alpha=0.5,linewidth=3)
    pylab.plot((lcl_low-411,lcl_low-411),yl,
               '-.b',alpha=0.5,linewidth=3)
    pylab.plot((lcl_up+411,lcl_up+411),yl,
               '--b',alpha=0.5,linewidth=3)
    pylab.plot((lcl_up-411,lcl_up-411),yl,
               '-.b',alpha=0.5,linewidth=3)

    #Plot common spectral lines
    x_Lyb = 1025.7*(1+redshift)
    x_Lya = 1215.7*(1+redshift)
    x_CIV = 1549.1*(1+redshift)
    x_AlIII = 1858.7*(1+redshift)
    x_FeII = 2600*(1+redshift)
    x_MgII = 2799.8*(1+redshift)
    x_MgI = 2852*(1+redshift)
    x_OII = 3727.61*(1+redshift)
    x_CalK = 3933.667*(1+redshift)
    x_CalH = 3968.472*(1+redshift)
    x_Hd = 4101.74*(1+redshift)
    x_Gband = 4305*(1+redshift)
    x_Hg = 4340.47*(1+redshift)
    x_Hb = 4861.33*(1+redshift)
    x_OIII_1 = 4960.3*(1+redshift)
    x_OIII_2 = 5008.24*(1+redshift)
    x_Mgb = 5176*(1+redshift)
    x_FeI = 5269*(1+redshift)
    x_NaD = 5893*(1+redshift)
    x_NII_1 = 6548.06*(1+redshift)
    x_Ha = 6562.799*(1+redshift)
    x_NII_2 = 6585.2*(1+redshift)
    x_SII = 6725.5*(1+redshift)

    # Plot dashed lines at the respective redshifts
    pylab.plot((x_Lyb,x_Lyb),yl,'--k')
    pylab.plot((x_Lya,x_Lya),yl,'--k')
    pylab.plot((x_CIV,x_CIV),yl,'--k')
    pylab.plot((x_AlIII,x_AlIII),yl,'--k')
    pylab.plot((x_FeII,x_FeII),yl,'--k')
    pylab.plot((x_MgII,x_MgII),yl,'--k')
    pylab.plot((x_MgI,x_MgI),yl,'--k')
    pylab.plot((x_OII,x_OII),yl,'--k')
    pylab.plot((x_CalK,x_CalK),yl,'--k')
    pylab.plot((x_CalH,x_CalH),yl,'--k')
    pylab.plot((x_Hd,x_Hd),yl,'--k')
    pylab.plot((x_Gband,x_Gband),yl,'--k')
    pylab.plot((x_Hg,x_Hg),yl,'--k')
    pylab.plot((x_Hb,x_Hb),yl,'--k')
    pylab.plot((x_OIII_1,x_OIII_1),yl,'--k')
    pylab.plot((x_OIII_2,x_OIII_2),yl,'--k')
    pylab.plot((x_Mgb,x_Mgb),yl,'--k')
    pylab.plot((x_FeI,x_FeI),yl,'--k')
    pylab.plot((x_NaD,x_NaD),yl,'--k')
    pylab.plot((x_NII_1,x_NII_1),yl,'--k')
    pylab.plot((x_Ha,x_Ha),yl,'--k')
    pylab.plot((x_NII_2,x_NII_2),yl,'--k')
    pylab.plot((x_SII,x_SII),yl,'--k')

    labeloff = 0.5
    pylab.text(lambda_central, labeloff*(yl[0]+yl[1]),
               '$\lambda_\mathrm{central}$'+'={0}'.format(lambda_central),
               horizontalalignment='right',verticalalignment='center',
               rotation='vertical',fontsize=16)

    if x_Lyb > xl[0] and x_Lyb < xl[1]:
        pylab.text(x_Lyb, labeloff*(yl[0]+yl[1]), 'Ly-beta', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_Lya > xl[0] and x_Lya < xl[1]:
        pylab.text(x_Lya, labeloff*(yl[0]+yl[1]), 'Ly-alpha', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_CIV > xl[0] and x_CIV < xl[1]:
        pylab.text(x_CIV, labeloff*(yl[0]+yl[1]), 'C IV', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_AlIII > xl[0] and x_AlIII < xl[1]:
        pylab.text(x_AlIII, labeloff*(yl[0]+yl[1]), 'Al III', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_FeII > xl[0] and x_FeII < xl[1]:
        pylab.text(x_FeII, labeloff*(yl[0]+yl[1]), 'Fe II', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_MgII > xl[0] and x_MgII < xl[1]:
        pylab.text(x_MgII, labeloff*(yl[0]+yl[1]), 'Mg II', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_MgI > xl[0] and x_MgI < xl[1]:
        pylab.text(x_MgI, labeloff*(yl[0]+yl[1]), 'Mg I', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_OII > xl[0] and x_OII < xl[1]:
        pylab.text(x_OII, labeloff*(yl[0]+yl[1]), '[O II]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_CalK > xl[0] and x_CalK < xl[1]:
        pylab.text(x_CalK, labeloff*(yl[0]+yl[1]), 'Cal K', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_CalH > xl[0] and x_CalH < xl[1]:
        pylab.text(x_CalH, labeloff*(yl[0]+yl[1]), 'Cal H', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_Hd > xl[0] and x_Hd < xl[1]:
        pylab.text(x_Hd, labeloff*(yl[0]+yl[1]), 'Hd', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_Gband > xl[0] and x_Gband < xl[1]:
        pylab.text(x_Gband, labeloff*(yl[0]+yl[1]), 'G-band', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_Hg > xl[0] and x_Hg < xl[1]:
        pylab.text(x_Hg, labeloff*(yl[0]+yl[1]), 'Hg', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_Hb > xl[0] and x_Hb < xl[1]:
        pylab.text(x_Hb, labeloff*(yl[0]+yl[1]), 'Hb', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_OIII_1 > xl[0] and x_OIII_1 < xl[1]:
        pylab.text(x_OIII_1, labeloff*(yl[0]+yl[1]), '[OIII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_OIII_2 > xl[0] and x_OIII_2 < xl[1]:
        pylab.text(x_OIII_2, labeloff*(yl[0]+yl[1]), '[OIII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_Mgb > xl[0] and x_Mgb < xl[1]:
        pylab.text(x_Mgb, labeloff*(yl[0]+yl[1]), 'Mg I(b)', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_FeI > xl[0] and x_FeI < xl[1]:
        pylab.text(x_FeI, labeloff*(yl[0]+yl[1]), 'Fe I', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_NaD > xl[0] and x_NaD < xl[1]:
        pylab.text(x_NaD, labeloff*(yl[0]+yl[1]), 'Na I (D)', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_NII_1 > xl[0] and x_NII_1 < xl[1]:
        pylab.text(x_NII_1, (labeloff-0.1)*(yl[0]+yl[1]), '[NII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_Ha > xl[0] and x_Ha < xl[1]:
        pylab.text(x_Ha, labeloff*(yl[0]+yl[1]), 'Ha', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_NII_2 > xl[0] and x_NII_2 < xl[1]:
        pylab.text(x_NII_2, (labeloff+0.1)*(yl[0]+yl[1]), '[NII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')
    if x_SII > xl[0] and x_SII < xl[1]:
        pylab.text(x_SII, labeloff*(yl[0]+yl[1]), '[SII]', horizontalalignment='right',verticalalignment='center', rotation='vertical')

    pylab.xlim(xl)
    frame1 = pylab.gca()
    frame1.axes.get_yaxis().set_visible(False)
    pylab.xlabel('$\lambda_{observed}$',fontsize=18)
    if filename:
        pylab.savefig(filename)
    pylab.show()
