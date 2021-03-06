import yoda
import numpy as np
import pylab as plt
from collections import OrderedDict
import re
import os

# routine that extracts all data from a top file
def extractPWHGData(fname,ynorm=1.0) :
        with open(fname) as f :
                lines = f.readlines()
        matching = False
        observable = None
        data = OrderedDict()
        for line in lines :
                if matching and re.match("^$",line) != None :
                        matching = False
                        data[observable] = np.array(data[observable])
                if matching :
                        # collapse blanks
                        line=' '.join(line.split())
                        line=re.sub("D","E",line)
                        try :
                                datarow = map(float,line.split(' '))
                        except ValueError as e :
                                print "\t{}".format(e)
                                print "\t{}".format(line)
                        data[observable].append([datarow[0], datarow[1], ynorm*datarow[2], ynorm*datarow[3]])
                match = re.match("^# ([^ ]*).*$",line)
                if match != None :
                        matching = True
                        observable = match.group(1)
                        data[observable]=[]
        return data


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
        
os.system("mkdir -p plots-NLO/MC_HJETSVBF")
NNLOJETprefix='NNLOJETdataplot/LH19VFHNLOdata/R04NLO.'
suffix='.dat'
yodaFiles=['yodafiles/fNLO-HW7-new' ]

yodaData =[ yoda.read(yodafile+'.yoda') for yodafile in yodaFiles]

color=['blue', 'limegreen', 'violet', 'orange', 'darkgreen']
for key in yodaData[0].keys():
        if ('RAW' in key) or (not 'MC_HJETSVBF' in key) or ('cross' in key) or ('deltaphi_jj_ATLAS' in key) or ('[W1]' in key):
                continue
        NNLOJETkey=key[13:]
        #for nnlojet ht -> htpt, not yet done at NLO        
        if('HT' in NNLOJETkey):
                NNLOJETkey='HTpt'+NNLOJETkey[2:] 
        NNLOJETfname=NNLOJETprefix+NNLOJETkey+'.dat'
        os.system("ls "+NNLOJETfname)
        nnj = extractNNJETHisto(NNLOJETfname)

        plt.suptitle(NNLOJETkey)
        ax1 = plt.subplot(211)
        ax2 = plt.subplot(212)
        
        ax1.errorbar(nnj[:,1], nnj[:,3], yerr=nnj[:,4], xerr=nnj[:,1]-nnj[:,0], color='r', ls='-', label='NNLOJET')
        ax2.errorbar(nnj[:,1], nnj[:,3]/nnj[:,3], yerr=nnj[:,4]/nnj[:,3], xerr=nnj[:,1]-nnj[:,0], color='r', ls='-', label='')
        # nplot[-1][0].set_linestyle('-')
        # nplot[-1][-1].set_linestyle('-')

        
        for yd in range(len(yodaData)):
                histo=yodaData[yd][key]
                xmin = np.array([bin.xMin for bin in histo.bins])
                xmax = np.array([bin.xMax for bin in histo.bins])
                xav  = np.array([bin.xMid for bin in histo.bins])
                binsize = xmax - xmin 
                yval = np.array([bin.area for bin in histo.bins])/binsize
                yerr = np.array([bin.relErr for bin in histo.bins])*yval
                ax1.errorbar(xav, yval, yerr=yerr, xerr=binsize/2, color=color[yd], ls='-', label=yodaFiles[yd][11:])
                ax2.errorbar(xav, yval/nnj[:,3], yerr=yerr/nnj[:,3], xerr=binsize/2, color=color[yd], ls='-', label='')

        ylims=ax2.get_ylim()
        ylims=[max(ylims[0],0.9), min(ylims[1],1.1)]
        ax2.set_ylim(ylims)
        ax1.legend()
        plt.savefig('plots-NLO/'+key+'.pdf')
        plt.close()
