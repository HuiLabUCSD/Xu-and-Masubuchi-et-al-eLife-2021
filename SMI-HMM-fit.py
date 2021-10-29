# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 18:39:31 2020

@author: MT

102220 add Count State Length section
091020 create
"""


"""
import data from excel
-> analyze with HMM
-> output 

***caution***
tp_array and dwell_array go backwards
"""

import numpy as np
from hmmlearn import hmm
import math
import pandas as pd
import matplotlib.pyplot as plt
import glob

# data information
target_protein = "shp2"
folder = glob.glob("%s/*.xlsx" %(target_protein))

path = "Trace.xlsx"


# import data
raw_data = pd.read_excel(path, dtype="float32")
raw_np = raw_data.values
del raw_data
    
# count the number of sample
plot_N = raw_np.shape[0]    # number of timepoints
sample_N = raw_np.shape[1]  # number of spots

# fit by HMM
model = hmm.GaussianHMM(n_components=2, covariance_type="full")

# prepare arrays to store calculated numbers
dwell_array = np.empty(0)   # storing dwell times
tp_array = np.empty(0)      # stroing transmat 
change = np.zeros(plot_N)        # state change points
state_len = np.zeros(plot_N)     # length of each states
on_len = np.empty(0)        # summary of the length of on states
off_len = np.empty(0)       # summary of the length of off states

for i in range(0,sample_N):
    
    target = raw_np[:,i].reshape(plot_N,1)
    
    model.fit(target)
    predict = model.predict(target)
    data_normalize = (target-min(model.means_))/(abs(model.means_[1]-model.means_[0]))
    
    if model.means_[0] > model.means_[1]:
        predict = 1 - predict
        model.transmat_ = np.flip(model.transmat_)
        
    # tp
    tp_array = np.append(model.transmat_, tp_array)
    
    #dwell time
    dwell = math.log(2) / model.transmat_[1,0]
    dwell_array = np.append(dwell, dwell_array)
    
    # on/off length
    for k in range(plot_N):
        if k == plot_N:
            change[k] = 1
        else:
            if predict[k] == predict[k-1]:
                change[k] = 0
            else:
                change[k] = 1
        if i == 0:
            state_len[k] = 1
        else:
            if change[k-1] == 0:
                state_len[k] = state_len[k-1]+1
            else:
                state_len[k] = 1
        if change[k] == 1:
            if predict[k] == 1:
                on_len = np.append(on_len, state_len[k])
            else:
                off_len = np.append(off_len, state_len[k])
    
    #plot data
    fig = plt.figure()
    plt.plot(data_normalize, marker="", label="obs", color="green")
    plt.plot(predict, marker="", label="fit", color="gray")
    plt.legend(loc = "upper right")
    fig.savefig("%d.png" %(i+1))
    
fig = plt.figure()
plt.hist(dwell_array*0.05, range=[0,20])
fig.savefig("%s_dwell.png" %(path))
"""
np.savetxt("%s_tp.csv" %(path), tp_array, delimiter=",")
np.savetxt("%s_dwell.csv" %(path), dwell_array, delimiter=",")
np.savetxt("%s_on_length.csv" %(path), on_len, delimiter=",")
np.savetxt("%s_off_length.csv" %(path), off_len, delimiter=",")
"""
np.savetxt("%s_predict.csv" %(path), predict, delimiter=",")
np.savetxt("%s_data_normalize.csv" %(path), data_normalize, delimiter=",")


