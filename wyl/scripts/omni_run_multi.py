#! /usr/bin/env python

import omnical, aipy, numpy, capo
import optparse, os, sys, glob
from astropy.io import fits
import pickle
from multiprocessing import Pool
from scipy.io.idl import readsav
#from IPython import embed

o = optparse.OptionParser()
o.set_usage('omni_run_multi.py [options] *uvcRRE/obsid')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--calpar',dest='calpar',type='string',default=None,
            help='Path and name of calpar file (txt or npz).')
o.add_option('--redinfo',dest='redinfo',type='string',default='',
            help='Path and name of .bin redundant info file.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .npz files. Include final / in path.')
o.add_option('--ba',dest='ba',default=None,
            help='Antennas to exclude, separated by commas.')
o.add_option('--flength',dest='flength',default=None,
             help='a threshold for baseline lengths to use, in meters')
o.add_option('--ftype', dest='ftype', default='', type='string',
            help='Type of the input file, .uvfits, or miriad, or fhd, to read fhd, simply type in the path/obsid')
o.add_option('--tave', dest='tave', default=False, action='store_true',
             help='choose to average data over time before calibration or not')
o.add_option('--gave', dest='gave', default=False, action='store_true',
             help='choose to average solution over time after calibration or not')
o.add_option('--iftxt', dest='iftxt', default=False, action='store_true',
            help='A switch to write the npz info to a ucla format txt file or not')
o.add_option('--iffits', dest='iffits', default=False, action='store_true',
            help='A switch to write the npz info to a ucla format fits file or not')
o.add_option('--removedegen',dest='removedegen',default=False,action='store_true',
             help='A switch to turn remove degen on')
o.add_option('--instru', dest='instru', default='mwa', type='string',
             help='instrument type')
opts,args = o.parse_args(sys.argv[1:])

#Dictionary of calpar gains and files
pols = opts.pol.split(',')
files = {}
#files=[]
g0 = {} #firstcal gains
g_scale = {'x': 1.0, 'y': 1.0}
if opts.calpar != None: #create g0 if txt file is provided
    fname = opts.calpar
    if fname.endswith('.txt'):
        f = open(fname,'r')
        Ntimes = []
        Nfreqs = []
        for line in f:
            temp = line.split(',')[:7]
            if temp[0].startswith('#'): continue
            temp2 = []
            for ii, s in enumerate(temp):
                if ii == 0: continue
                elif s.strip() == 'EE': s = 'xx' #need to check the convension
                elif s.strip() == 'NN': s = 'yy'
                elif s.strip() == 'EN': s = 'xy'
                elif s.strip() == 'NE': s = 'yx'
                temp2.append(s)
            if not temp2[2].strip() in pols: continue
            temp3 = [temp2[2], int(temp2[0]), float(temp2[3]), float(temp2[1]), float(temp2[4]), float(temp2[5])]  #temp3=[pol,ant,jds,freq,real,imag]
            if not temp3[2] in Ntimes: Ntimes.append(temp3[2])
            if not temp3[3] in Nfreqs: Nfreqs.append(temp3[3])
            if not g0.has_key(temp3[0][0]):
                g0[temp3[0][0]] = {}
            if not g0[temp3[0][0]].has_key(temp3[1]):
                g0[temp3[0][0]][temp3[1]] = []
            gg = complex(temp3[4],temp3[5])
            g0[temp3[0][0]][temp3[1]].append(gg.conjugate()/abs(gg))
        for pp in g0.keys():
            for ant in g0[pp].keys():
                g0[pp][ant] = numpy.array(g0[pp][ant])
                g0[pp][ant] = g0[pp][ant].reshape(len(Ntimes),len(Nfreqs))
    elif fname.endswith('.npz'):
        for pp,p in enumerate(pols):
            g0[p[0]] = {}   #obs(or jds).pol.fc.npz
            fpname = fname.split('.')
            fpname[-3] = p
            fpname = '.'.join(fpname)
            print '   Reading: ', fpname
            cp = numpy.load(fpname)
            for i in cp.keys():
                if i[0].isdigit():
                    g0[p[0]][int(i[:-1])] = cp[i] / numpy.abs(cp[i])
    elif fname.endswith('.fits'):
        g0 = capo.omni.fc_gains_from_fits(opts.calpar)
        for key1 in g0:
            for key2 in g0[key1]:
                g0[key1][key2] /= numpy.abs(g0[key1][key2])
    elif fname.endswith('.sav'):
        cal = readsav(opts.calpar,python_dict=True)
        fqfl = numpy.zeros((128,384),dtype=bool)
        for ff in range(384):
            if ff%16==0 or ff%16==15: fqfl[:,ff]=True
        g = cal['cal']['GAIN'][0]
        g0['x'] = {}
        g0['y'] = {}
        gx = numpy.ma.masked_array(g[0],fqfl,fill_value=1.0)
        gy = numpy.ma.masked_array(g[1],fqfl,fill_value=1.0)
        gnan = numpy.where(numpy.isnan(numpy.mean(gx,axis=1)))[0]
        g_scale['x'] = numpy.nanmean(numpy.abs(gx[56:-1]))
        g_scale['y'] = numpy.nanmean(numpy.abs(gy[56:-1]))
        for nn in gnan:
            gx[nn] = g_scale['x']
            gy[nn] = g_scale['y']
        for ii in range(0,cal['cal']['N_TILE'][0]):
            g0['x'][ii] = gx[ii].filled()
            g0['y'][ii] = gy[ii].filled()
    else:
        raise IOError('invalid calpar file')

#if not provided, will initiate g0 with units in the reading file part

for filename in args:
    files[filename] = {}
    if opts.ftype == 'miriad':
        for p in pols:
            fn = filename.split('.')
            fn[-2] = p
            files[filename][p] = '.'.join(fn)
    elif opts.ftype == 'uvfits':
        files[filename][filename] = filename + '.uvfits'
    elif opts.ftype == 'fhd':
        obs = filename + '*'
        filelist = glob.glob(obs)
        if len(pols) == 1:
            p = pols[0]
            if p == 'xx':
                try: filelist.remove(filename + '_vis_YY.sav')
                except: pass
            elif p == 'yy':
                try: filelist.remove(filename + '_vis_XX.sav')
                except: pass
            else:
                raise IOError('do not support cross pol')
            files[filename][p] = filelist
        else:
            files[filename][filename] = filelist
    else:
        raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')

#################################################################################################
def calibration(infodict):#dict=[filename, g0, timeinfo, d, f, ginfo, freqs, polar, cal, calpar]
    filename = infodict['filename']
    g0 = infodict['g0']
    pos = infodict['position']
    p = infodict['pol']
    d = infodict['data']
    f = infodict['flag']
    ginfo = infodict['ginfo']
    freqs = infodict['freqs']
    timeinfo = infodict['timeinfo']
    calpar = infodict['calpar']
    ex_ants = infodict['ex_ants']
    print 'Getting reds from calfile'
    print 'generating info:'
    filter_length = None
    if not opts.flength == None: filter_length = float(opts.flength)
    info = capo.omni.pos_to_info(pos, pols=list(set(''.join([p]))), filter_length=filter_length, ex_ants=ex_ants, crosspols=[p])

    ### Omnical-ing! Loop Through Compressed Files ###

    print '   Calibrating ' + p + ': ' + filename
    
    #if txt file or first cal is not provided, g0 is initiated here, with all of them to be 1.0
    if calpar == None:
        if not g0.has_key(p[0]): g0[p[0]] = {}
        for iant in range(0, ginfo[0]):
            if opts.tave: g0[p[0]][iant] = numpy.ones((1,ginfo[2]))
            else: g0[p[0]][iant] = numpy.ones((ginfo[1],ginfo[2]))
    elif calpar.endswith('.npz') or calpar.endswith('.sav'):
        for key in g0[p[0]].keys():
            g0_temp = g0[p[0]][key]
            g0[p[0]][key] = numpy.resize(g0_temp,(ginfo[1],ginfo[2]))


    t_jd = timeinfo['times']
    t_lst = timeinfo['lsts']

    data,wgts,xtalk = {}, {}, {}
    m2,g2,v2 = {}, {}, {}
    data = d #indexed by bl and then pol (backwards from everything else)

    wgts[p] = {} #weights dictionary by pol
    for bl in f:
        i,j = bl
        wgts[p][(j,i)] = wgts[p][(i,j)] = numpy.logical_not(f[bl][p]).astype(numpy.int)
    print '   Logcal-ing' 
    m1,g1,v1 = capo.omni.redcal(data,info,gains=g0, removedegen=opts.removedegen) #SAK CHANGE REMOVEDEGEN
    print '   Lincal-ing'
    m2,g2,v2 = capo.omni.redcal(data, info, gains=g1, vis=v1, uselogcal=False, removedegen=opts.removedegen)
    if opts.tave:
        for a in g2[p[0]].keys():
            g2[p[0]][a] = numpy.resize(g2[p[0]][a],(ginfo[1],ginfo[2]))
        for bl in v2[p].keys():
            v2[p][bl] = numpy.resize(v2[p][bl],(ginfo[1],ginfo[2]))
    if opts.gave:
        for a in g2[p[0]].keys():
            gmean = numpy.mean(g2[p[0]][a],axis=0)
            g2[p[0]][a] = numpy.resize(gmean,(ginfo[1],ginfo[2]))
    xtalk = capo.omni.compute_xtalk(m2['res'], wgts) #xtalk is time-average of residual
    ############# correct the center of each coarse band if instrument is mwa ######################
    if opts.instru == 'mwa':
        for a in g2[p[0]].keys():
            for ff in range(0,384):
                if ff%16==8:
                    g2[p[0]][a][:,ff] = (g2[p[0]][a][:,ff-1]+g2[p[0]][a][:,ff+1])/2
    ############# To rescale solutions if not remove degen ####################
    if not opts.removedegen:
        g_rescale = 0
        ncount = 0
        flag = numpy.zeros((ginfo[1],ginfo[2]),dtype=bool)
        if opts.calpar.endswith('.sav'):
            blacklist = [0,15]
            for ki in range(0,384):
                if ki%16 in blacklist: flag[:,ki] = True
            flag[0] = True
            flag[53] = True
            flag[54] = True
            flag[55] = True
        for ka in g2[p[0]].keys():
            gd = numpy.ma.masked_array(numpy.abs(g2[p[0]][ka]),flag)
            g_rescale += numpy.mean(gd)
            ncount += 1
        g_rescale /= ncount
        g_rescale /= g_scale[p[0]]
        v_rescale = g_rescale*g_rescale
        refn = min(g2[p[0]].keys())
        for ka in g2[p[0]].keys(): ### take tile 1001 as reference tile ###
            g2[p[0]][ka] /= g_rescale
            g2[p[0]][ka] /= (g2[p[0]][refn]/numpy.abs(g2[p[0]][refn]))
        for kb in v2[p].keys():
            v2[p][kb] *= v_rescale
    if opts.calpar.endswith('.sav'):
        for key in g0[p[0]].keys():
            if not g2[p[0]].has_key(key):
                g2[p[0]][key] = g0[p[0]][key]
    ###########################################################################################
    m2['history'] = 'OMNI_RUN: '+''.join(sys.argv) + '\n'
    m2['jds'] = t_jd
    m2['lsts'] = t_lst
    m2['freqs'] = freqs

    if opts.ftype == 'miriad':
        npzname = opts.omnipath+'.'.join(filename.split('/')[-1].split('.')[0:4])+'.npz'
    else:
        npzname = opts.omnipath+filename.split('/')[-1]+'.'+ p +'.npz'

    print '   Saving %s'%npzname
    capo.omni.to_npz(npzname, m2, g2, v2, xtalk)
    return npzname
#######################################################################################################

exec('from %s import antpos as _antpos'% opts.cal)
for f,filename in enumerate(args):

    npzlist = []
    infodict = {}
    filegroup = files[filename]
    info_dict = []
    print "  Reading data: " + filename
    if opts.ftype == 'miriad':
        for p in pols:
            dict0 = capo.wyl.uv_read_omni([filegroup[p]], filetype = 'miriad', antstr='cross', p_list=[p], tave=opts.tave)
            infodict[p] = dict0[p]
            infodict[p]['filename'] = filegroup[p]
            infodict['name_dict'] = dict0['name_dict']
    else:
        infodict = capo.wyl.uv_read_omni([filegroup[key] for key in filegroup.keys()], filetype=opts.ftype, antstr='cross', p_list=pols, tave=opts.tave)
        for p in pols:
            infodict[p]['filename'] = filename
    print "  Finish reading."
    for p in pols:
        if opts.calpar == None:
            infodict[p]['g0'] = {}
        else:
            infodict[p]['g0'] = {}
            infodict[p]['g0'][p[0]] = g0[p[0]]
        infodict[p]['calpar'] = opts.calpar
        infodict[p]['position'] = _antpos
        if opts.ba:
            for a in opts.ba.split(','):
                if not int(a) in infodict[p]['ex_ants']:
                    infodict[p]['ex_ants'].append(int(a))
        ex_ants = sorted(infodict[p]['ex_ants'])
        print '   Excluding antennas:', ex_ants

        info_dict.append(infodict[p])
    print "  Start Parallelism:"
    par = Pool(2)
    npzlist = par.map(calibration, info_dict)
    par.close()
    name_dict = infodict['name_dict']

    if opts.iftxt: #if True, write npz gains to txt files
        scrpath = os.path.abspath(sys.argv[0])
        pathlist = os.path.split(scrpath)[0].split('/')
        repopath = '/'.join(pathlist[0:-1])+'/'
        print '   Writing to txt:'
        capo.wyl.writetxt(npzlist, repopath, ex_ants=ex_ants, name_dict=name_dict)
        print '   Finish'

    if opts.iffits: #if True, write npz gains to fits files
        scrpath = os.path.abspath(sys.argv[0])
        pathlist = os.path.split(scrpath)[0].split('/')
        repopath = '/'.join(pathlist[0:-1])+'/'
        print '   Writing to fits:'
        capo.wyl.writefits(npzlist, repopath, ex_ants=ex_ants, name_dict=name_dict)
        print '   Finish'



