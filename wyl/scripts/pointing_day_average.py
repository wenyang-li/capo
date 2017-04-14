import numpy as np
import sys, glob, optparse
import capo.omni as omni
from astropy.io import fits

o = optparse.OptionParser()
o.set_usage('pointing_day_average.py [options]')
o.set_description(__doc__)
o.add_option('--metafits', dest='metafits', default='/users/wl42/data/wl42/Nov2016EoR0/', type='string',
             help='path to metafits files, look for pointings')
opts,args = o.parse_args(sys.argv[1:])

delays = {
'0,5,10,15,1,6,11,16,2,7,12,17,3,8,13,18':-5,
'0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15':-4,
'0,3,6,9,0,3,6,9,0,3,6,9,0,3,6,9':-3,
'0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6':-2,
'0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3':-1,
'0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0':0,
'3,2,1,0,3,2,1,0,3,2,1,0,3,2,1,0':1,
'6,4,2,0,6,4,2,0,6,4,2,0,6,4,2,0':2,
'9,6,3,0,9,6,3,0,9,6,3,0,9,6,3,0':3,
'12,8,4,0,13,9,5,1,14,10,6,2,15,11,7,3':4,
'15,10,5,0,16,11,6,1,17,12,7,2,18,13,8,3':5,
}

opts,args = o.parse_args(sys.argv[1:])
#p = sys.argv[1]
pols = ['xx','yy']
for p in pols:
    fn=glob.glob('./*'+p+'.npz')
    g = {}
    nfiles = {}
    for f in fn:
        obs = f.split('/')[-1].split('.')[0]
        if not obs[0].isdigit(): continue
        meta, gains, vismdl, xtalk = omni.from_npz(f)
        metafits = '../'+obs+'.metafits'
        hdu = fits.open(metafits)
        day = int(obs)/86400
        suffix = str(day)+'_'+str(delays[hdu[0].header['DELAYS']])
        if not g.has_key(suffix): g[suffix]={p[0]:{}}
        if not nfiles.has_key(suffix): nfiles[suffix]=0
        nfiles[suffix]+=1
        for a in gains[p[0]].keys():
            if not g[suffix][p[0]].has_key(a): g[suffix][p[0]][a] = gains[p[0]][a]
            else: g[suffix][p[0]][a] += gains[p[0]][a]
    print nfiles
    for suf in g.keys():
        for a in g[suf][p[0]].keys():
            g[suf][p[0]][a] /= nfiles[suf] 
        outfn = 'omniave_'+str(suf)+'.'+p+'.npz'
        omni.to_npz(outfn, meta, g[suf], vismdl, xtalk)

