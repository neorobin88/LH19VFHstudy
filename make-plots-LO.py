import yoda
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict
import re
import os

def extractNNJETHisto(fname, ynorm=1.0) :
        with open(fname) as f :
                lines = f.readlines()
        matching = False      
        histo = []
        for line in lines :
                if(not ("#" in line)):
                        # collapse blanks
                        line=' '.join(line.split())
                        line=re.sub("D","E",line)
                        try :
                                historow = map(float,line.split(' '))
                        except ValueError as e :
                                print "\t{}".format(e)
                                print "\t{}".format(line)
                        histo.append([historow[0], historow[1], historow[2], ynorm*historow[3], ynorm*historow[4]])
        histo=np.array(histo)
        return histo
        
os.system("mkdir -p plots-LO/MC_HJETSVBF")
NNLOJETprefix=['NNLOJETdataplot/1stpaperrun/VFH/Data/VFHPDF4LHC15R04LO.']
suffix='.dat'
yodaFiles=['yodafiles-v2/VBF-LO' ]

yodaData =[ yoda.read(yodafile+'.yoda') for yodafile in yodaFiles]

color=['blue', 'limegreen', 'violet', 'orange', 'darkgreen']
for key in yodaData[0].keys():
        if ('RAW' in key) or (not 'MC_HJETSVBF' in key) or ('cross' in key) or ('dr10' in key) or ('dr07' in key) or( 'incl' in key) or ('fb' in key) or ('dy01' in key) or ('vbfvh' in key) or ('x' in key) or ('pt2_pt1' in key) or ('njets' in key) or ('pt3_pt1' in key) or ('log' in key):
                continue
        else:
                print key

        

                
        print key[13:]
        NNLOJETkey=key[13:]
        if("rstudy" in NNLOJETkey):
                NNLOJETkey=NNLOJETkey[12:]

        if((NNLOJETkey=="pth200_dy10_pth") or (NNLOJETkey=="pth500_dy10_pth")):
                continue

        
        print NNLOJETkey
        NNLOJETfname=NNLOJETprefix[0]+NNLOJETkey+'.dat'
        nnj = extractNNJETHisto(NNLOJETfname)
        
        
        
                
        plt.suptitle(NNLOJETkey)
        ax1 = plt.subplot(211)
        ax2 = plt.subplot(212)

        
        ref=np.array(nnj[:,3])
        ax1.errorbar(nnj[:,1], nnj[:,3], yerr=nnj[:,4], xerr=nnj[:,1]-nnj[:,0], color='r', ls='-', label='NNLOJET')
        ax2.errorbar(nnj[:,1], nnj[:,3]/ref, yerr=nnj[:,4]/ref, xerr=nnj[:,1]-nnj[:,0], color='r', ls='-', label='')
        # nplot[-1][0].set_linestyle('-')
        # nplot[-1][-1].set_linestyle('-')
        print ref[0]
        
        for yd in range(len(yodaData)):
                histo=yodaData[yd][key]
                xmin = np.array([bin.xMin for bin in histo.bins])
                xmax = np.array([bin.xMax for bin in histo.bins])
                xav  = np.array([bin.xMid for bin in histo.bins])
                binsize = xmax - xmin 
                yval = np.array([bin.area for bin in histo.bins])/binsize
                yerr = np.array([bin.relErr for bin in histo.bins])*yval
                if ("phi" in key):
                        #temporary fix to plot also nnlojet data
                        print "PHI IN KEY", key
                        xav=nnj[:,1]
                        yval=yval/np.pi
                        yerr=yerr/np.pi
                ax1.errorbar(xav, yval, yerr=yerr, xerr=binsize/2, color=color[yd], ls='-', label='LO')
                ax2.errorbar(xav, yval/ref, yerr=yerr/yval, xerr=binsize/2, color=color[yd], ls='-', label='')
        ylims=ax2.get_ylim()
        ylims=[max(ylims[0],0.9), min(ylims[1],1.1)]
        ax2.set_ylim(ylims)
        ax1.legend()
        plt.savefig('plots-LO/'+key+'.pdf')
        plt.close()
