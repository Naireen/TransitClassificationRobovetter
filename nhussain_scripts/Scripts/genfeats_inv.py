import numpy as np
import pandas as pd
import sklearn
import matplotlib 
from matplotlib import pyplot as plt
from scipy.special import erfcinv

data_names = np.loadtxt("../Training/simulation/invert/invert.ls", delimiter=",",dtype=str)

#print data_names[1000:1010]

#data_names =  np.core.defchararray.add(data_names[:-4] , "_inv.txt")
#print data_names[1000:1010]
#data_names = np.core.defchararray.add(data_names, data[:,1])
#data_names = np.core.defchararray.add(data_names, "_inv.blsanal")


for name in data_names[20000:]:
    
    #print name[:-3] + "blsanal", name[:-4]+".txt"
    try:
        first_bls = np.loadtxt("../Training/simulation/invert/" + name[:-3]+"blsanal")[0, 1:]
        first_lc = np.loadtxt("../Training/simulation/invert/"+name[:-4]+".txt")
    except:
        print "Faulty file"
        continue 
    #print name[: -3]+"blsanal", name[:-4]+".txt" 
    #fold it in phase space
    phase = ((first_lc[:,0] -first_bls[1])/first_bls[0])%1
    depth = first_bls[5]
    transit_duration = first_bls[6]
    ingress_duration = first_bls[7]*first_bls[6]*1.5


    #find white noise, and remove primary signal
    white_noise = np.std(first_lc[:,1][(phase>transit_duration/2) & (phase<1-transit_duration/2)])
    mean_out_transit_flux = np.mean(first_lc[:,1][(phase>transit_duration/2) & (phase<1-transit_duration/2)])
    #remove in transit portion

    cleaned_vals = np.empty(shape=[0,3])
    #remove transit signal
    slope = depth /ingress_duration

    
    for i, vals in enumerate(first_lc):
        vals =vals.reshape(1,3)
        if phase[i] <transit_duration-ingress_duration or phase[i]>1-transit_duration+ingress_duration:
             vals[0][1] = vals[0][1]-depth
        elif ((phase[i] > transit_duration-ingress_duration  and phase[i]<transit_duration)) : 
             vals[0][1] -= (ingress_duration - phase[i]) *slope

        elif ((phase[i]>1-transit_duration) and (phase[i]<1-transit_duration+ingress_duration)):
             vals[0][1] -= (ingress_duration - (1 - phase[i])) *slope

        cleaned_vals = np.append(cleaned_vals, vals, axis = 0)



    #using transit width, search through phase space to find secondary and tertiary signals
    comb = zip(phase, cleaned_vals[:,1]) 
    comb.sort()
    phase_time, total_flux = zip(*comb)
    total_flux = np.asarray(total_flux).reshape((len(total_flux),1))

    mod_sig_sec_dv = 0
    mod_sig_pos_dv = 0
    sec_loc = 1
    phase_time = np.asarray(phase_time)

    cut_phase =  phase_time[(phase_time > transit_duration) & (phase_time<(1-transit_duration))]
    cut_flux = total_flux[(phase_time > transit_duration) & (phase_time<(1-transit_duration))]

    for i, vals in enumerate(cut_phase):
        section = cut_flux[(cut_phase>cut_phase[i]) & (cut_phase<cut_phase[i]+transit_duration*2)]
        depth_two = np.mean(section) - mean_out_transit_flux
        if depth_two == max( mod_sig_sec_dv, depth_two):
            sec_loc = phase_time[i]
        mod_sig_sec_dv = max(mod_sig_sec_dv, depth_two)
        mod_sig_pos_dv = min(mod_sig_pos_dv, depth_two) 




    for i, vals in enumerate(cut_phase):
        if (vals< sec_loc+transit_duration) and (vals>sec_loc-transit_duration):
            cut_flux[i] -= mod_sig_sec_dv

    mod_sig_ter_dv = 0
    for i, vals in enumerate(cut_phase):
        section = cut_flux[(cut_phase>cut_phase[i]) & (cut_phase<cut_phase[i]+transit_duration*2)]
        sig = np.mean(section) - mean_out_transit_flux
        mod_sig_ter_dv = max(mod_sig_sec_dv, sig)
                                        
    mod_sig_sec_dv = mod_sig_sec_dv - mean_out_transit_flux


    #red noise
    #red noise^2 = white noise /npoiint + red noise
    #inverse relationship
    #try with three timescales, one, two, and threee t durations
    #look at portions outside of transit?
    red_noise_1 = []
    red_noise_2 = []
    red_noise_3 = []

    red_noise_1_p = []
    red_noise_2_p = []
    red_noise_3_p = []

    for i, vals in enumerate(cut_phase[:-1]):
        section1 = cut_flux[(cut_phase>cut_phase[i]) & (cut_phase<cut_phase[i]+transit_duration)]
        section2 = cut_flux[(cut_phase>cut_phase[i]) & (cut_phase<cut_phase[i]+transit_duration*2)]
        section3 = cut_flux[(cut_phase>cut_phase[i]) & (cut_phase<cut_phase[i]+transit_duration*3)]
        red_noise_1_p.append(section1.shape[0])
        red_noise_2_p.append(section2.shape[0])
        red_noise_3_p.append(section3.shape[0])
        red_noise_1.append(np.var(section1))
        red_noise_2.append(np.var(section2))
        red_noise_3.append(np.var(section3))
    
    N_p_1 = np.median(red_noise_1_p)
    N_p_2 = np.median(red_noise_2_p)
    N_p_3 = np.median(red_noise_3_p)

    white_noise_2 = (np.mean(red_noise_1) - np.mean(red_noise_2))* (N_p_1 * N_p_2) / (N_p_1+N_p_2)

    red_noise = np.sqrt(np.mean(red_noise_1) - white_noise_2/N_p_1)

    #from scipy.special import erfcinv

    period = first_bls[0]
    n_trans = first_bls[14]


    mod_sig_pri_dv = depth
    mod_fred_dv = red_noise / white_noise # ratio between the two
    mod_fa1_dv = np.sqrt(2) * erfcinv(transit_duration/(period* n_trans))
    mod_fa2_dv = np.sqrt(2) * erfcinv(transit_duration/(period))

    modshiftval1_dv = mod_sig_pri_dv/mod_fred_dv - mod_fa1_dv
    modshiftval2_dv = mod_sig_pri_dv - mod_sig_ter_dv-mod_fa2_dv
    modshiftval3_dv = mod_sig_pri_dv - mod_sig_pos_dv - mod_fa2_dv
    modshiftval4_dv = mod_sig_sec_dv / mod_fred_dv - mod_fa1_dv
    modshiftval5_dv = mod_sig_sec_dv - mod_sig_ter_dv - mod_fa2_dv
    modshiftval6_dv = mod_sig_sec_dv-mod_sig_pos_dv - mod_fa2_dv


    vals = [mod_sig_pri_dv, mod_sig_sec_dv, mod_sig_ter_dv, mod_sig_pos_dv, mod_fa1_dv, mod_fa2_dv, mod_fred_dv, modshiftval1_dv, modshiftval2_dv, modshiftval3_dv, modshiftval4_dv, modshiftval5_dv, modshiftval6_dv, white_noise, red_noise]

    total_Vals = np.concatenate([vals, first_bls])

    head_txt = """mod_sig_pri_dv mod_sig_sec_dv mod_sig_ter_dv mod_sig_pos_dv mod_fa1_dv mod_fa2_dv mod_fred_dv modshiftval1_dv modshiftval2_dv modshiftval3_dv modshiftval4_dv modshiftval5_dv modshiftval6_dv white_noise red_noise BLS_Period_1_0 BLS_Tc_1_0 BLS_SN_1_0 BLS_SR_1_0 BLS_SDE_1_0 BLS_Depth_1_0 BLS_Qtran_1_0 BLS_Qingress_1_0 BLS_OOTmag_1_0 BLS_i1_1_0 BLS_i2_1_0 BLS_deltaChi2_1_0 BLS_fraconenight_1_0 BLS_Npointsintransit_1_0 BLS_Ntransits_1_0 BLS_Npointsbeforetransit_1_0 BLS_Npointsaftertransit_1_0 BLS_Rednoise_1_0 BLS_Whitenoise_1_0 BLS_SignaltoPinknoise_1_0 """
    np.savetxt("InvModShiftVals/" + name[:-8] +"_feats.txt", total_Vals, header=head_txt)
    #break
