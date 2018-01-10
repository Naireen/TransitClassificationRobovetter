#!/usr/bin/env python
import numpy as np
import pandas as pds
import matplotlib
from matplotlib import pyplot as plt
from dataio import readcolumn
import os

def cal_distance(p1, p2):
    #print p1, p2
    paarr = (np.arange(5)+1)*p1
    pbarr = (np.arange(5)+1)*p2
    xv, yv = np.meshgrid(1./paarr, 1./pbarr)
    diff = np.abs(xv-yv)
    #print np.min(diff)
    #print np.where(diff == np.min(diff))
    xarr, yarr = np.where(diff == np.min(diff))
    pa = p1*(xarr[0]+1)
    pb = p2*(yarr[0]+1)
    err = 2*np.abs(1./pa-1./pb)/(1./pa+1./pb)
    #print err
    return err 

def main():
    inlist = "kplr_dr25_inj1_shortperiod.txt"
    kepid = []; readcolumn(kepid, 1, inlist); kepid = np.array(kepid)
    injperiod = []; readcolumn(injperiod, 3, inlist); injperiod = np.array(injperiod)
    depth = []; readcolumn(depth, 6, inlist); depth = np.array(depth)
    epoch = []; readcolumn(epoch, 4, inlist); epoch = np.array(epoch)
    expectedMES = []; readcolumn(expectedMES, 14, inlist); expectedMES = np.array(expectedMES)
    ntransits = []; readcolumn(ntransits, 5, inlist); ntransits = np.array(ntransits)
    # plotting everything and check
    print len(expectedMES[expectedMES>6])
    plt.hist(expectedMES)
    plt.xlabel("Expected MES")
    plt.show()

    plt.hist(injperiod)
    plt.xlabel("Injected Period")
    plt.show()
    
    print max(depth), min(depth)
    plt.hist(depth)
    plt.xlabel("Injected Depth")
    plt.xlim([10, 10000])
    plt.xscale('log')
    plt.show()

    scaled_MES = expectedMES*(27/injperiod)**0.5/ntransits**0.5
    print len(scaled_MES[scaled_MES>6])
    plt.hist(scaled_MES)
    plt.xlabel("scaled_MES")
    plt.show()
    #exit()

    for i in xrange(len(kepid)):
        for j in xrange(17):

            infile = "../simulation/primaryinj1/kplr%.9d-%d_prim_ltf.blsanal" % (int(kepid[i]), j) 
            if os.path.exists(infile):
                blsanal_info = np.loadtxt(infile)
                #print blsanal_info.shape
                detect_period = blsanal_info[:,1]
                SNR = blsanal_info[:, -1]
                white_noise = blsanal_info[:, -2]
                detect_snr = 0
                detect_index = -1
                detect_distance = -1
                #print infile
                for k in xrange(len(detect_period)):
                    distance = cal_distance(detect_period[k], injperiod[i])
                    if distance<0.01:
                        detect_index = k
                        detect_distance = distance
                        
                        if SNR[k] > detect_snr:
                            detect_snr = SNR[k]
                if detect_index == -1:
                    print int(kepid[i]),",",j,",",injperiod[i],",",depth[i],",",epoch[i],",",expectedMES[i],",",scaled_MES[i],",",ntransits[i],",",
                    for k in xrange(blsanal_info.shape[1]-1):
                        print blsanal_info[0, k+1], ",",
                    print detect_distance
                else :
                    print int(kepid[i]),",",j,",",injperiod[i],",",depth[i],",",epoch[i],",",expectedMES[i],",",scaled_MES[i],",",ntransits[i],",",
                    for k in xrange(blsanal_info.shape[1]-1):
                        print blsanal_info[detect_index, k+1], ",",
                    print detect_distance

                #print infile
                #print blsanal_info
                #break 
        #break

    return

if __name__ == '__main__':
    main()
