# -----------------------------------------------------------------------
# This script sets up a 1-D seismic velocity model and calls a Fortran
# subroutine to calculate a receiver function. It then performs the following tasks
#       1.) plots predicted and observed receiver function at a reference model together with velocity model
#       2.) Optionally plots predicted misfit along a chosen axis through reference model
#
# Most parameters relating to the Receiver Function are set within the Fortran routine.
#
# The fortran source was interfaced with ctypes
#
# Fortran routines for Receiver function calculation by T. Shibutani
#
# M. Sambridge and Juerg Hauser,
# RSES, Nov., 2024.
#
# -----------------------------------------------------------------------
import pyrf96
import numpy as np

##################################################################################

vtype = 2
if vtype == 0:
    # Set up a velocity reference model in Voronoi cell format
    velmod = np.zeros([7, 3])
    velmod[0] = [8.370596, 3.249075, 1.7]
    velmod[1] = [17.23163, 3.001270, 1.7]
    velmod[2] = [1.9126695e-02, 2.509443, 1.7]
    velmod[3] = [19.78145, 3.562691, 1.7]
    velmod[4] = [41.73066, 4.225965, 1.7]
    velmod[5] = [14.35261, 2.963322, 1.7]
    velmod[6] = [49.92358, 4.586726, 1.7]
elif vtype == 1:
    # Set up a velocity reference model in Layer Thickness  format
    velmod = np.zeros([13, 3])
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
elif vtype == 2:
    # Set up a velocity reference model in Interface depth format
    velmod = np.zeros([13, 3])
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

# time1, RFo = np.loadtxt("RF_obs.dat", unpack=True) # read in observed Receiver function
time1, RFo = pyrf96.rfcalc(velmod, sn=0.25, mtype=vtype, seed=61254557)
time2, RFp = pyrf96.rfcalc(velmod, mtype=vtype)

pyrf96.plot_RFs(time1, RFo, time2, RFp)  # plot a pair of RFs in a single frame

pv, pd = pyrf96.plotRFm(
    velmod, time1, RFo, time2, RFp, mtype=vtype, dmax=60.0
)  # plot a pair of RFs with velocity model

# Set up data Covariance matrix
# 2.5, 0.01 = correlation half-width and amplitude of Gaussian noise
ndata = len(RFp)
Cdinv = pyrf96.InvDataCov(2.5, 0.01, ndata)  # Calculate Data covariance matrix

# Calculate waveform misfit for reference model

res = RFo - RFp
mref = np.dot(res, np.transpose(np.dot(Cdinv, res))) / 2.0
print(" Waveform misfit of reference model", mref)

# Calculate waveform misfit profile along chosen model parameter axis

iparam = [1, 0]  # choose model parameter axis
# 1,0 = depth layer 2
# 4,0 = depth layer 5
# 1,1 = velocity layer 2
# 3,1 = velocity layer 4

# Set up axis limits
limits0 = [
    [0.0, 2.0],
    [0.0, 2.0],
    [0.0, 2.0],
    [0.0, 3.0],
    [0.0, 3.0],
    [0.0, 2.0],
    [0.0, 2.0],
]
limits1 = [
    [60.0, 4.0],
    [60.0, 4.0],
    [60.0, 4.0],
    [60.0, 5.0],
    [60.0, 5.0],
    [60.0, 4.0],
    [60.0, 5.0],
]
x0 = limits0[iparam[0]][iparam[1]]  # axis lower limit
x1 = limits1[iparam[0]][iparam[1]]  # axis upper limit
nints = 512  # Number of points along axis
x = np.linspace(x0, x1, nints)  # set axis values
xtrue = velmod[iparam[0], iparam[1]]  # reference parameter value along axis

profile = True  # switch to calculate and plot misfit profile
# Calculate profile
if profile:
    misfit = np.zeros(nints)
    for i in range(nints):
        velmod[iparam[0], iparam[1]] = x[i]
        t, RF = pyrf96.rfcalc(velmod)  # calculate receiver function for velocity model
        res = RFo - RF
        misfit[i] = np.dot(res, np.transpose(np.dot(Cdinv, res))) / 2.0 - mref

    # Plot misfit profile
    pyrf96.plot_misfit_profile(x[1:-2], misfit[1:-2], xtrue, iparam)
