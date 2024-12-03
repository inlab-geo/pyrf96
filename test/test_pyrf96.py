import numpy as np
import pyrf96

##################################################################################

vtype = 2
if(vtype==0):
    velmod = np.zeros([7,3]) # Set up a velocity reference model in Voronoi cell format
    velmod[0] = [8.370596, 3.249075, 1.7]
    velmod[1] = [17.23163, 3.001270, 1.7]
    velmod[2] = [1.9126695E-02, 2.509443, 1.7]
    velmod[3] = [19.78145, 3.562691, 1.7]
    velmod[4] = [41.73066, 4.225965, 1.7]
    velmod[5] = [14.35261, 2.963322, 1.7]
    velmod[6] = [49.92358, 4.586726, 1.7]
elif(vtype==1):
    velmod = np.zeros([13,3]) # Set up a velocity reference model in Layer Thickness  format
    velmod[0] = [20.0, 3.4600, 1.73]
    velmod[1] = [15.0, 3.8500, 1.73]
    velmod[2] = [42.5, 4.4850, 1.73]
    velmod[3] = [42.5, 4.4950, 1.73]
    velmod[4] = [45.0, 4.5045, 1.73]
    velmod[5] = [45.0, 4.5145, 1.73]
    velmod[6] = [50.0, 4.5660, 1.73]
    velmod[7] = [50.0, 4.6525, 1.73]
    velmod[8] = [50.0, 4.7395, 1.73]
    velmod[9] = [50.0, 4.8265, 1.73]
    velmod[10] = [50.0, 5.1330, 1.73]
    velmod[11] = [50.0, 5.2390, 1.73]
    velmod[12] = [0.0, 5.3450, 1.73]
elif(vtype==2):
    velmod = np.zeros([13,3]) # Set up a velocity reference model in Interface depth format
    velmod[0] = [20.0, 3.4600, 1.73]
    velmod[1] = [35.0, 3.8500, 1.73]
    velmod[2] = [77.5, 4.4850, 1.73]
    velmod[3] = [120.0, 4.4950, 1.73]
    velmod[4] = [165.0, 4.5045, 1.73]
    velmod[5] = [210.0, 4.5145, 1.73]
    velmod[6] = [260.0, 4.5660, 1.73]
    velmod[7] = [310.0, 4.6525, 1.73]
    velmod[8] = [360.0, 4.7395, 1.73]
    velmod[9] = [410.0, 4.8265, 1.73]
    velmod[10] = [460.0, 5.1330, 1.73]
    velmod[11] = [510.0, 5.2390, 1.73]
    velmod[12] = [510.0, 5.3450, 1.73]

# calculate and plot waveforms

#time1, RFo = np.loadtxt("RF_obs.dat", unpack=True) # read in observed Receiver function
time1, RFo = pyrf96.rfcalc(velmod,sn=0.25,mtype=vtype,seed=61254557)
#print('time1',time1)
#print(' RFo',RFo)
time2, RFp = pyrf96.rfcalc(velmod,mtype=vtype)

#print('time2',time2)
#print(' RFp',RFp)

pyrf96.plot_RFs(time1,RFo,time2,RFp) # plot a pair of RFs in a single frame

