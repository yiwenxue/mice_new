import time
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
from scipy.signal import argrelextrema

print ("Import finished!")

def sizeof_fmt(num, suffix='B'):
    '''output data sizes in a nice way, i.e. B, KiB, MiB, GiB, ...
    by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254,
    modified'''
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)

def readMiceData(batch, name):
    ''' Read the preprocessed data from the dir
    ''' 
    if debug:
        print('{}/batch{}/{}.Activity.txt'.format(path, batch, name))
        print('{}/batch{}/{}.Temperature.txt'.format(path, batch, name))
    a = np.genfromtxt('{}/batch{}/{}.Activity.txt'.format(path, batch, name))
    b = np.genfromtxt('{}/batch{}/{}.Temperature.txt'.format(path, batch, name))

    if a[0,0] != b[0,0]:
        start = max(a[0,0], b[0,0])
    else:
        start = a[0,0]
    if a[-1,0] != b[-1,0]:
        end = min(a[-1,0], b[-1,0])
    else:
        end = a[-1,0]

    a = a[(a[:,0]<=end)]
    a = a[(a[:,0]>=start)]
    b = b[(b[:,0]<=end)]
    b = b[(b[:,0]>=start)]
    
    return a, b

def getMidpeek(data):
    density = np.histogram(data, bins=100 ,density=True)
    
    weight = density[0]
    tempx = density[1][0:-1]

    temp_offset = np.mean(tempx)
    tempx = tempx - temp_offset
    
    # fit with a smooth line
    z = np.polyfit (tempx, weight, 20)
    fitfunc = np.poly1d(z)
    
    x = np.arange(tempx[0],tempx[-1],(tempx[-1]-tempx[0])/100)
    y = fitfunc(x)
    
    plt.plot(x,y)
    plt.plot(tempx,weight)
    plt.savefig("fit_vs_ori.pdf")
    
    maxid = argrelextrema(y, np.greater)
    maxes = y[maxid]
    a = np.sort(maxes, kind="quicksort")
    maxima1 = a[-1]
    maxima2 = a[-2]
    maxid1 = np.where(y == maxima1)[0]
    maxid2 = np.where(y == maxima2)[0]
    # maxid1, maxid2
    maxtemp1 = x[maxid1]
    maxtemp2 = x[maxid2]
    midpeek = (maxtemp1[0] + maxtemp2[0]) / 2.0 + temp_offset
    
    return midpeek


winSize = 6
winNum = len(temp) // 6

debug = 0
path = '.'
# The mouse to test
batch = 3
mouse = '12Otx2'

if __name__ == 'main':

    # Readin the meta data
    metaData = pd.read_csv("./mice", index_col=0)
    # Read test data
    act, temp = readMiceData(batch, mouse)

    # data overview
    plt.plot(temp[:, 0], temp[:, 1])
    plt.savefig("overview_temp.pdf")
    plt.plot(act[:, 0], act[:, 1])
    plt.savefig("overview_act.pdf")

    # calculate temp average

    tmp = temp[:,1].reshape(winNum,winSize)
    tempave = np.mean(tmp, axis=1)
    sb.distplot(tempave)
    plt.savefig("density_plot.pdf")

    midpeek = getMidpeek(tempave)
    temp_tmp = temp[:,1]
    awakeidx, sleepidx = index_by_value(tmep[:,1], midpeek)
    
    
    
def index_by_value(temp, midpeek):
    sleepidx = np.nonzero(temp < midpeek)[0]
    awakeidx = np.nonzero(temp >= midpeek)[0]
    return awakeidx, sleepidx

def idx_to_stage(idx):
    stages = []
    last = idx[0]
    i=1
    while i < len(idx):
        current = idx[i]
        if (current - idx[i-1] > 1):
            stages.append([last, idx[i-1]])
            last = current
        i+=1 
    return stages