import capo.omni as omni
import capo.wyl as wyl
from scipy.io.idl import readsav
import numpy as np, omnical, aipy
from multiprocessing import Pool
import optparse, os, sys, glob

o = optparse.OptionParser()
o.set_usage('FHD_lin.py [options] obsid')
aipy.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--fhdout',dest='fhdout',type='string',default='/users/wl42/scratch/FHD_out/',
             help='Path to FHD output.')
o.add_option('--vstr',dest='vstr',type='string',default='PhaseII_EoR0_2',
             help='FHD firstpass version string')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
             help='Path to save .npz files. Include final / in path.')
opts,args = o.parse_args(sys.argv[1:])

obsid = args[0]
fhd_path = opts.fhdout + 'fhd_' + opts.vstr + '/'
pols = opts.pol.split(',')

#************************** Generate Models **************************************************
print "   Generating Model"
sol = readsav(fhd_path+'calibration/'+obsid+'_cal.sav',python_dict=True)
g0 = {}
g0['x'] = {}
g0['y'] = {}
for jj in range(0,384):
    if jj%16 in [0,15]:
        sol['cal']['GAIN'][0][0][:,jj] = 1
        sol['cal']['GAIN'][0][1][:,jj] = 1
for ii in range(sol['cal']['GAIN'][0][0].shape[0]):
    g0['x'][ii] = sol['cal']['GAIN'][0][0][ii]
    g0['y'][ii] = sol['cal']['GAIN'][0][1][ii]
vis_list = glob.glob(fhd_path+'vis_data/'+obsid+'*')
metalist = glob.glob(fhd_path+'metadata/'+obsid+'*')
filelist = vis_list + metalist
uv_model = wyl.data_fhd()
uv_model.read_data_only(filelist,use_model=True)
exec('from %s import antpos'% opts.cal)
vis_model = {}
vis_data = uv_model.data_array.reshape(uv_model.Ntimes,uv_model.Nbls,uv_model.Nfreqs,uv_model.Npols)
a1 = uv_model.ant_1_array[:uv_model.Nbls]
a2 = uv_model.ant_2_array[:uv_model.Nbls]

#********************************* lincal function *******************************************
def lincal(datadict):
    p = datadict['pol']
    g1 = {}
    g1[p[0]] = g0[p[0]]
    d = datadict['data']
    f = datadict['flag']
    ginfo = datadict['ginfo']
    freqs = datadict['freqs']
    timeinfo = datadict['timeinfo']
    ex_ants = datadict['ex_ants']
    info = omni.pos_to_info(antpos, pols=list(set(''.join([p]))), ex_ants=ex_ants, crosspols=[p])
    for key in g1[p[0]].keys():
        g1_temp = g1[p[0]][key]
        g1[p[0]][key] = np.resize(g1_temp,(ginfo[1],ginfo[2]))
    reds = info.get_reds()
    ubl = []
    for r in reds:
        ubl.append(r[0])
    ind = np.where(uv_model.polarization_array==aipy.miriad.str2pol[p])[0][0]
    v1 = {}
    v1[p] = {}
    for ii in range(0,uv_model.Nbls):
        if (a1[ii],a2[ii]) in ubl: v1[p][(a1[ii],a2[ii])] = vis_data[:,ii,:,ind]
        elif (a2[ii],a1[ii]) in ubl: v1[p][(a2[ii],a1[ii])] = vis_data[:,ii,:,ind].conj()
    print '   Lincal-ing'
    m2,g2,v2 = omni.redcal(d, info, gains=g1, vis=v1, uselogcal=False, removedegen=False)
    wgts = {}
    wgts[p] = {}
    for bl in f:
        i,j = bl
        wgts[p][(j,i)] = wgts[p][(i,j)] = np.logical_not(f[bl][p]).astype(np.int)
    xtalk = omni.compute_xtalk(m2['res'], wgts) #xtalk is time-average of residual
    m2['history'] = 'FHD_lin: '+''.join(sys.argv) + '\n'
    m2['jds'] = timeinfo['times']
    m2['freqs'] = freqs
    npzname = opts.omnipath+obsid+'.'+ p +'.FHDlin.npz'
    print '   Saving %s'%npzname
    omni.to_npz(npzname, m2, g2, v2, xtalk)
    return npzname

#******************************** Load Data and calibrate ************************************
npzlist = []
info_dict = []
print "   Loading Data"
Data = wyl.uv_read_omni([obsid+'.uvfits'],filetype='uvfits',antstr='cross',p_list=pols)
for p in pols:
    info_dict.append(Data[p])
par = Pool(2)
npzlist = par.map(lincal, info_dict)


