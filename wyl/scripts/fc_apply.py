import numpy as np,aipy
import pyuvdata.uvdata as uvd
import sys,optparse

o = optparse.OptionParser()
o.set_usage('fc_apply.py [options] obsid(do not include .uvfits)')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,pol=True)
o.add_option('--omnipath',dest='omnipath',default='%s.npz',type='string',
            help='Format string (e.g. "path/%s.fc.npz", where you actually type the "%s") which converts the input file name to the omnical npz path/file.')
o.add_option('--outtype', dest='outtype', default='uvfits', type='string',
             help='Type of the output file, .uvfits, or miriad, or fhd')
o.add_option('--intype', dest='intype', default='uvfits', type='string',
             help='Type of the input file, .uvfits or fhd')
opts,args = o.parse_args(sys.argv[1:])
pols = opts.pol.split(',')
files = {}
for filename in args:
    if opts.intype == 'fhd':
        files[filename] = []
        filelist = glob.glob(opts.fhdpath+'/vis_data/'+filename+'*')+glob.glob(opts.fhdpath+'/metadata/'+filename+'*')
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

for f,filename in enumerate(args):
    newfile = filename + '_fc' + '.uvfits'
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
    for ip,p in enumerate(pols):
        pid = np.where(pollist == aipy.miriad.str2pol[p])[0][0]
        omnifile = opts.omnipath % (filename.split('/')[-1]+'.'+p)
        print '  Reading and applying:', omnifile
        fcfile = np.load(omnifile)
        gains = {p[0]:{}}
        for k in fcfile.keys():
            if k[0].isdigit():
                a = int(k[:-1])
                gains[p[0]][a] = fcfile[k][0]
        for ii in range(0,Nblts):
            a1 = uvi.ant_1_array[ii]
            a2 = uvi.ant_2_array[ii]
            p1,p2 = p
            try: uvi.data_array[:,0][:,:,pid][ii] /= gains[p1][a1]
            except(KeyError): pass
            try: uvi.data_array[:,0][:,:,pid][ii] /= gains[p2][a2].conj()
            except(KeyError): pass
    if opts.outtype == 'uvfits':
        print 'writing:' + newfile
        uvi.write_uvfits(newfile,spoof_nonessential=True)
        print 'saving ' + newfile
    
