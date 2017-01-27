#! /usr/bin/env python
# Do not support miriad

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys, glob
import uvdata.uvdata as uvd

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_apply_fhd.py [options] obsid(do not include .uvfits) or zen.jds.pol.uv')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,pol=True)
o.add_option('--xtalk',dest='xtalk',default=False,action='store_true',
            help='Toggle: apply xtalk solutions to data. Default=False')
o.add_option('--fit',dest='fit',default=False,action='store_true',
             help='Toggle: do a 7th order polyfit to sols over the band. Default=False')
o.add_option('--amp',dest='amp',default=False,action='store_true',
             help='Toggle: apply only amplitude solutions')
o.add_option('--omnipath',dest='omnipath',default='%s.npz',type='string',
            help='Format string (e.g. "path/%s.npz", where you actually type the "%s") which converts the input file name to the omnical npz path/file.')
o.add_option('--npz',dest='npz',default=None,type='string',
             help='specify npz file names, (format: path/name, without .pol.npz), otherwise find npz according to obsid and omnipath')
o.add_option('--outtype', dest='outtype', default='uvfits', type='string',
             help='Type of the output file, .uvfits, or miriad, or fhd')
o.add_option('--intype', dest='intype', default=None, type='string',
             help='Type of the input file, .uvfits or fhd')
o.add_option('--instru', dest='instru', default='mwa', type='string',
             help='instrument type. Default=mwa')
opts,args = o.parse_args(sys.argv[1:])


def fit_func(x,z):
    sum = numpy.zeros((x.size))
    for ii in range(z.size):
        sum *= x
        sum += z[ii]
    return sum

#File Dictionary
pols = opts.pol.split(',')
files = {}

#create a dictionary of file lists
for filename in args:
    if opts.intype == 'fhd':
        files[filename] = []
        obs = filename + '*'
        filelist = glob.glob(obs)
        files[filename] = filelist
    elif opts.intype == 'uvfits':
        files[filename] = filename + '.uvfits'
    elif opts.intype == 'miriad':
        files[filename] = {}
        for p in pols:
            fn = filename.split('.')
            fn[-2] = p
            files[filename][p] = '.'.join(fn)
    else:
        raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')

#start processing
for f,filename in enumerate(args):
    
    #create an out put filename
    if opts.outtype == 'uvfits':
        suffix = 'O'
        if opts.fit:
            suffix = 'fit' + suffix
        if opts.amp:
            suffix = suffix + 'amp'
        newfile = filename + '_' + suffix + '.uvfits'
    if os.path.exists(newfile):
        print '    %s exists.  Skipping...' % newfile
        continue

    #read in the file
    print '  Reading', files[filename]
    uvi = uvd.UVData()
    if opts.intype == 'fhd':
        uvi.read_fhd(files[filename])
    elif opts.intype == 'uvfits':
        uvi.read_uvfits(files[filename])
    Nblts = uvi.Nblts
    Nfreqs = uvi.Nfreqs
    Nbls = uvi.Nbls
    pollist = uvi.polarization_array
    freqs = uvi.freq_array[0]

    #find npz for each pol, then apply
    for ip,p in enumerate(pols):
        if not opts.npz == None:
            omnifile = opts.npz + '.' + p + '.npz'
        else:
            omnifile = opts.omnipath % (filename.split('/')[-1]+'.'+p)
        print '  Reading and applying:', omnifile
        _,gains,_,xtalk = capo.omni.from_npz(omnifile) #loads npz outputs from omni_run
#********************** if choose to fit a smooth function to sols ***************************
        if opts.fit:
            for key in gains[p[0]].keys():
                if opts.instru == 'mwa':
                    fqs,amp = [],[]
                    gamp = numpy.abs(numpy.mean(gains[p[0]][key][1:53],axis=0))
                    for nn in range(0,384):
                        if nn%16 in [0,1,2,13,14,15]: continue
                        fqs.append(freqs[nn])
                        amp.append(gamp[nn])
                    fqs = numpy.array(fqs)
                    amp = numpy.array(amp)
                else:
                    fqs = freqs
                    amp = numpy.abs(numpy.mean(gains[p[0]][key],axis=0))
                fit_coeff = numpy.polyfit(fqs,amp,7)
                gains[p[0]][key] = gains[p[0]][key]/numpy.abs(gains[p[0]][key])*fit_func(freqs,fit_coeff)
#*********************************************************************************************
        pid = numpy.where(pollist == aipy.miriad.str2pol[p])[0][0]
        for ii in range(0,Nblts):
            a1 = uvi.ant_1_array[ii]
            a2 = uvi.ant_2_array[ii]
            p1,p2 = p
            ti = ii/Nbls
                #for jj in range(0,Nfreqs):
            if opts.xtalk:
                try: uvi.data_array[:,0][:,:,pid][ii] -= xtalk[p][(a1,a2)]
                except(KeyError):
                    try: uvi.data_array[:,0][:,:,pid][ii] -= xtalk[p][(a2,a1)].conj()
                    except(KeyError): pass
            if opts.amp:
                try: uvi.data_array[:,0][:,:,pid][ii] /= numpy.abs(gains[p1][a1][ti])
                except(KeyError): pass
                try: uvi.data_array[:,0][:,:,pid][ii] /= numpy.abs(gains[p2][a2][ti])
                except(KeyError): pass
            else:
                try: uvi.data_array[:,0][:,:,pid][ii] /= gains[p1][a1][ti]
                except(KeyError): pass
                try: uvi.data_array[:,0][:,:,pid][ii] /= gains[p2][a2][ti].conj()
                except(KeyError): pass

    #write file
#uvi.history = ''
    if opts.outtype == 'uvfits':
        print 'writing:' + newfile
        uvi.write_uvfits(newfile,spoof_nonessential=True)
        print 'saving ' + newfile


