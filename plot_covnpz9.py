#! /usr/bin/env python
import numpy as n, pylab as p, sys, aipy as a
import optparse

colors = 'kbrgcmy'

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, pol=True)
opts, args = o.parse_args(sys.argv[1:])

p.rcParams['legend.fontsize'] = 6

filegroups = {}
for filename in args:
    basefile = filename.split('__')[0]
    filegroups[basefile] = filegroups.get(basefile, []) + [filename]
srcest_bm, srcest_ant, srcest_bl = {}, {}, {}
for basefile in filegroups.keys():
    for filename in filegroups[basefile]:
        fwords = filename[:-len('.npz')].split('__')
        try: f = n.load(filename)
        except(IOError): continue
        if fwords[1] == 'times': times = f['times']
        elif fwords[1] == 'afreqs': afreqs= f['freqs']
        elif fwords[1] == 'srcest_bm':
            for k in f.files: srcest_bm[k] = f[k]
        elif fwords[1] == 'srcest_ant':
            k = fwords[2]
            srcest_ant[k] = {}
            for i in f.files: srcest_ant[k][int(i)] = f[i]
        elif fwords[1] == 'srcest_bl':
            k = fwords[2]
            srcest_bl[k] = {}
            for bl in f.files: srcest_bl[k][int(bl)] = f[bl]

srcs = srcest_bm.keys(); srcs.sort()
if opts.cal != None:
    srclist = []
    for src in srcs:
        radec = src.split('_')
        if len(radec) == 2:
            src = a.phs.RadioFixedBody(ra=radec[0], dec=radec[1], name=src)
        srclist.append(src)
    cat = a.cal.get_catalog(opts.cal, srclist)
    aa = a.cal.get_aa(opts.cal, afreqs)
else: cat = {}

#srcs = ['cyg'] + srcs

norm=1
for cnt, k in enumerate(srcs):
    #d,t = srcest_bm[k], times
    #order = n.argsort(t)
    #d,t = d.take(order, axis=0), t.take(order)
    #I = 1
    #shape = (int(t.shape[0]/I), I)
    #ints = shape[0] * shape[1]
    #d,t = d[:ints], t[:ints]
    #d.shape,t.shape = shape + d.shape[1:], shape
    #d,t = n.average(d, axis=1), n.average(t, axis=1)
    #d *= norm

    ## Calculate beam response
    #bm = []
    #for jd in t:
    #    aa.set_jultime(jd)
    #    cat[src].compute(aa)
    #    bm.append(aa[0].bm_response(cat[src].get_crds('top'), pol=opts.pol)**2)
    #bm = n.array(bm).squeeze()
    #spec = n.sum(d*bm, axis=0)/n.sum(bm**2, axis=0)
    #if i == 0 and src == 'cyg':
    #    norm = cat['cyg'].jys / spec
    #    norm.shape = (1,norm.size)
    #    continue
    #ind, flx = n.polyfit(n.log10(afreqs/.150), n.log10(spec), deg=1)
    #
    #q = n.average((d-n.average(d))*(bm - n.average(bm))) / n.std(d) / n.std(bm)
    #print '%25s: FLX=%6.1f IND=%+4.2f Q=%+4.2f' % (src, 10**flx, n.round(ind,2), n.round(q,2))
    ##d /= bm
    ##_f = flx[1:-1] - .5 * (flx[2:] + flx[:-2])
    ##q = n.sum(n.abs(_f)) / n.sum(n.abs(flx[1:-1]))
    ##if q < .5: continue
    color = colors[cnt%len(colors)]
    #p.subplot(211)
    #p.semilogy(times, n.average(n.abs(srcest_bm[k]), axis=1), color+':', label=k)
    #for i in srcest_ant[k]:
    #    p.semilogy(t-.0002, n.average(n.abs(n.sqrt(srcest_bm[k])+srcest_ant[k][i])**2, axis=1), color+',', label=k)
    d,w = 0.,0.
    bmsqrt = n.sqrt(srcest_bm.get(k,0.))
    dw = n.ones_like(bmsqrt)
    for bl in srcest_bl[k]:
        i,j = a.miriad.bl2ij(bl)
        ai = srcest_ant.get(k,{}).get(i,0.)
        aj = srcest_ant.get(k,{}).get(j,0.)
        d += (bmsqrt + ai) * n.conj(bmsqrt + aj) + srcest_bl[k][bl]
        w += dw
    d /= w
    p.semilogy(times, n.average(n.abs(d), axis=1), color+'-', label=k)
    p.ylim(.1,1e5)

    #p.subplot(212)
    #p.loglog(afreqs, spec, color+',', label=k)
    #p.loglog(afreqs, 10**n.polyval([ind,flx], n.log10(afreqs/.150)), color+'-', label=k)
    #p.xlim(afreqs[0], afreqs[-1])
    #p.ylim(10,1e5)

#p.subplot(211)
#p.legend(loc='best')
p.show()

