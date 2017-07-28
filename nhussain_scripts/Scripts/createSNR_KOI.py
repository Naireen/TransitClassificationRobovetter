import numpy as np
import pandas as pd
import os

#file_name = "../nhussain_scripts/InvModShiftVals/10129857-10_feats.txt"
file_name = os.listdir("KOI_ModShiftVals/.")
data  = np.loadtxt(file_name, dtype = str, delimiter=",")#.reshape(data.shape[0], 5)
#print data.shape
#print data[0]
#data = data.reshape((data.shape[0],5))
#data = data[:, :-1]
#id, segment, depth, period
#print data.shape
#concate the first two, create list of file names
#values = np.core.defchararray.add(data[:,0], "-")
#values = np.core.defchararray.add(values, data[:, 1])
#data = data.astype("float64")#print data[0,:]
#print type(data[1,3])
#print values



SNR = np.zeros(data.shape[0])
#error = np.empty(data.shape[0])
#concatenate the first two values
for i, val in enumerate(data):
    print val
    lc_data = np.loadtxt("../Training/simulation/primary/"+val[:-9]+"ltf.lc")
    break
    depth =float(data[i,2])**2 #depth
    meidan = np.median(lc_data[:,1]) # flux
    error = np.median(np.abs(lc_data[:,1]-meidan)**2)**0.5
    transits = float(26)//float(data[i,3]) +1
    SNR[i] = depth*np.sqrt(transits)/(error*np.sqrt(2))
    if i %999 ==1:
        print i

np.savetxt("SNR_KOIS_510.txt", SNR)
