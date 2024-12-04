#-----------------------------------------------------------------------
# This script is a c-types interface between Receiver function routines and python
# 
# Most parameters relating to the Receiver Function are set within the Fortran routine.
#
# Fortran routines for Receiver function calculation by T. Shibutani
#
# M. Sambridge and Juerg Hauser, 
# RSES, Dec., 2024.
#
#-----------------------------------------------------------------------
from typing import Literal;

import os
import glob
import platform
import numpy
import numpy.ctypeslib
import ctypes
import numpy as np
from matplotlib import pyplot as plt

librf96=ctypes.cdll.LoadLibrary(glob.glob(os.path.dirname(__file__)+"/rf96*.so")[0])

def rfcalc(model,sn=0.0,mtype=0,fs=25.0,gauss_a=2.5,water_c=0.0001,angle=35.0,time_shift=5.0,ndatar=626,v60=8.043,seed=1): # Calculate Data covariance matrix for evaluation of waveform fit
    """
	Calculate synthetic surface wave dispersion curves for a given earth model, with optional addition of noise added in frequency domain.

	This is a slim Fortran wrapper around RF.F90 from Shibutani et al. (1996), which uses the Thomson-Haskell matrix formulation. 
	Further details may be found.


	Args:
    	model (np.ndarray)               : Triplet defining layered model. meaning depends on mytpe, with shape (npts,3).
    	sn (float,optional)              : Signal to noise ratio used to add correlated Gaussian noise to output.
    	mtype (int, optional)            : Indicator for format of velocity model (default=0) 
                                       	   model(1,i) is Vs velocity of layer i; 
                                           model(2,i) is vpvs ratio of layer i;
                                           mtype = 0 -> model(0,i) is the depth of Voronoi nuclei defining layer i;
                                           mtype = 1 -> model(0,i) is the thickness of layer i;
                                           mtype = 2 -> model(0,i) is depth of lower interface of layer i;
        fs (float, optional)             : Sampling frequency (default=25 samples/s)
        gauss_a (float, optional)        : Number 'a' defining the width of the gaussian filter in the deconvolution (default=2.5) 
        water_c (float, optional)        : Water level used in deconvolution (default=0.0001) 
        angle (float, optional)          : Angle in degrees of incoming teleseismic plane wave from vertical (default=35 degrees)
        time_shift (float, optional)     : Time shift before the first p pusle (default=5s)   
        ndatar (int,optional)            : Number of time time steps of output waveform
        v60 (float,optional)             : P-wave velocity (km/s) needed to compute the ray parameter from angle (default=8.043 km/s)
        seed (int,optional)              : Random set for noise gneration

    Returns:
        time (np.array, size ndatar )  : Time series time in seconds.
        wdata (np.array, size ndatar)  : The Receiver function amplitude.
    """

    npt = np.shape(model)[0]
    model_f = numpy.asfortranarray(model,dtype=numpy.float32)
    model_c = model_f.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    mtype_f = ctypes.c_int(mtype)
    mtype_c = ctypes.byref(mtype_f) 
    sn_f = ctypes.c_float(sn)
    sn_c = ctypes.byref(sn_f) 
    fs_f = ctypes.c_float(fs)
    fs_c = ctypes.byref(fs_f) 
    gauss_af = ctypes.c_float(gauss_a)
    gauss_ac = ctypes.byref( gauss_af) 
    water_cf = ctypes.c_float(water_c)
    water_cc = ctypes.byref( water_cf) 
    angle_f = ctypes.c_float(angle)
    angle_c = ctypes.byref( angle_f) 
    time_shift_f = ctypes.c_float(time_shift)
    time_shift_c = ctypes.byref(time_shift_f) 
    ndatar_f = ctypes.c_int(ndatar)
    ndatar_c = ctypes.byref(ndatar_f) 
    v60_f = ctypes.c_float(v60)
    v60_c = ctypes.byref(v60_f) 
    seed_f = ctypes.c_int(seed) 
    seed_c = ctypes.byref(seed_f)
    npt_f = ctypes.c_int(npt) 
    npt_c = ctypes.byref(npt_f)
    # returned variables
    time_f = numpy.asfortranarray(numpy.zeros(ndatar),dtype=numpy.float32)
    time_c = time_f.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    wdata_f = numpy.asfortranarray(numpy.zeros(ndatar),dtype=numpy.float32)
    wdata_c = wdata_f.ctypes.data_as(ctypes.POINTER(ctypes.c_float))

    if(sn==0.0):
        #a,b = rfm.rfcalc_nonoise(model,mtype,fs,gauss_a,water_c,angle,time_shift,ndatar,v60)
        librf96.rfcalc_nonoise(model_c,mtype_c,fs_c,gauss_ac,water_cc,angle_c,time_shift_c,ndatar_c,v60_c,npt_c,time_c,wdata_c)
    else:
        #a,b = rfm.rfcalc_noise(model,mtype,sn,fs,gauss_a,water_c,angle,time_shift,ndatar,v60,seed)
        librf96.rfcalc_noise(model_c,mtype_c,sn_c,fs_c,gauss_ac,water_cc,angle_c,time_shift_c,ndatar_c,v60_c,seed_c,npt_c,time_c,wdata_c)
        #print(time_f)
        #print(wdata_f)
    a = time_f
    b = wdata_f
    return a,b

def v2mod(model,vmin=2.4,vmax=4.7,dmin=0.0,dmax=60.0): # Transform Voronoi nucleus representation to (depth vel) plot format
    #a,b,c,d,e = rfm.voro2mod(model)
    # map variables for ctypes
    model_f = numpy.asfortranarray(model,dtype=numpy.float32)
    model_c = model_f.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    npt = np.shape(model)[0]
    npt_f = ctypes.c_int(npt)
    npt_c = ctypes.byref(npt_f)
    h_f = numpy.asfortranarray(numpy.zeros(npt),dtype=numpy.float32)
    h_c = h_f.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    beta_f = numpy.asfortranarray(numpy.zeros(npt),dtype=numpy.float32)
    beta_c = beta_f.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    vpvs_f = numpy.asfortranarray(numpy.zeros(npt),dtype=numpy.float32)
    vpvs_c = vpvs_f.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    qa_f = numpy.asfortranarray(numpy.zeros(npt),dtype=numpy.float32)
    qa_c = qa_f.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    qb_f = numpy.asfortranarray(numpy.zeros(npt),dtype=numpy.float32)
    qb_c = qb_f.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
    librf96.voro2mod(model_c,npt_c,h_c,beta_c,vpvs_c,qa_c,qb_c)

    px = np.zeros([2*len(qa_f)])
    py = np.zeros([2*len(qa_f)])
    #a,b,c,d,e = rfm.voro2mod(model)
    py[1::2],py[2::2],px[0::2],px[1::2] = list(np.cumsum(h_f)),list(np.cumsum(h_f))[:-1],beta_f[:],beta_f[:]
    py[-1] = dmax
    a,b,c,d,e = None,None,None,None,None
    return np.cumsum(h_f),beta_f,vpvs_f,qa_f,qb_f,px,py

##################################################################################
def InvDataCov(width,Ar,ndata): # Calculate Data covariance matrix for evaluation of waveform fit
    sigsq = width**2
    Arsq = Ar*Ar # waveform noise variance    
    Cd = np.zeros((ndata,ndata))
    for k1 in range(ndata):
        for k2 in range(ndata):
            Cd[k1][k2] = np.exp(-0.5*(k1-k2)**2/sigsq)
    U, s, V = np.linalg.svd(Cd)
    pmax = [ n for n,a in enumerate(s) if a<0.000001 ][0]
    Sinv = np.diag(1.0/s[:pmax])
    Cdinv = np.dot(V.T[:,:pmax], np.dot(Sinv, U.T[:pmax,:]))
    return Cdinv/Arsq

##################################################################################
def InvDataCovSub(width,Ar,ndata,ind): # Calculate Data sub-covariance matrix for evaluation of waveform fit
    sigsq = width**2
    Arsq = Ar*Ar # waveform noise variance    
    Cd = np.zeros((ndata,ndata))
    for k1 in range(ndata):
        for k2 in range(ndata):
            Cd[k1][k2] = np.exp(-0.5*(k1-k2)**2/sigsq)
    Cdused = Cd[ind,:][:,ind]
    U, s, V = np.linalg.svd(Cdused)
    pmax = [ n for n,a in enumerate(s) if a>0.000001 ][-1]+1
    Sinv = np.diag(1.0/s[:pmax])
    Cdinv = np.dot(V.T[:,:pmax], np.dot(Sinv, U.T[:pmax,:]))
    return Cdinv/Arsq        
##################################################################################
def plot_misfit_profile(x,misfit,xtrue,iparam): # Plot calculated and observed RF waveform 
    fig, ax = plt.subplots()
    
    if(iparam[1]==1): strx = "Layer "+repr(iparam[0]+1)+" velocity (km/s)" # plot axis labels
    if(iparam[1]==0): strx = "Layer "+repr(iparam[0]+1)+" node depth (km)"

    mismin = np.min(misfit)
    mismax = np.max(misfit)
    xt = [xtrue,xtrue]
    yt = [mismin,mismin+0.1*(mismax-mismin)]        
    ax.plot(x, misfit, 'k-')
    ax.plot(xt,yt,'r-',label='True')
    if(iparam[1] == 0):
        ax.set_xlabel('Node Depth (km)') 
    else:
        ax.set_xlabel('Velocity (km)/s')
    ax.set_xlabel(strx)
    #ax.set_ylabel('-Log L(m)')
    #ax.set_title("-Log-Likelihood through reference model")
    ax.set_ylabel('Misfit')
    ax.set_title("Waveform misfit profile")
    ax.grid(True)
    xt2 = [x[np.argmin(misfit)],x[np.argmin(misfit)]] 
    ax.plot(xt2,yt,'c-',label='Minimum') 
    plt.title("Misfit profile along single axis")
    plt.legend()
    plt.show()
##################################################################################
def l2mod(model,vmin=2.4,vmax=4.7,dmin=0.0,dmax=60.0): # Transform layer thickness representation to (depth vel) plot format
    px = np.zeros([2*len(model)])
    py = np.zeros([2*len(model)])
    py[1::2],py[2::2],px[0::2],px[1::2] = list(np.cumsum(model.T[0])),list(np.cumsum(model.T[0]))[:-1],model.T[1],model.T[1]
    py[-1] = dmax
    return px,py
##################################################################################
def d2mod(model,vmin=2.4,vmax=4.7,dmin=0.0,dmax=60.0): # Transform depth representation to (depth vel) plot format
    px = np.zeros([2*len(model)])
    py = np.zeros([2*len(model)])
    py[1::2],py[2::2],px[0::2],px[1::2] = list(model.T[0]),list(model.T[0])[:-1],model.T[1],model.T[1]
    py[-1] = dmax
    return px,py
##################################################################################
def plot_RFs(time1,RFo,time2,RFp,string="Observed and predicted receiver functions",filename=None):
    fig, ax = plt.subplots()
    plt.title(string)
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.grid(True)
    plt.plot(time1, RFo, 'k-', label='Observed')
    plt.plot(time2, RFp, 'r-', label='Predicted')
    plt.legend()
    if(filename is not None): plt.savefig(filename)
    plt.show()
##################################################################################
def plotRFm(velmod,time1,RFo,time2,RFp,mtype = 0,vmin=2.4,vmax=4.7,dmin=0.0,dmax=60.0,
            filename=None,string="Observed and predicted receiver functions"): # plot velocity model and receiver function
    f, (a0, a1) = plt.subplots(1,2, figsize=(12,4), gridspec_kw = {'width_ratios':[1, 3]})

    a1.set_title(string)
    a1.set_xlabel("Time (s)")
    a1.set_ylabel("Amplitude")
    #if((k==3) & (j==0)):a1.set_ylim(-0.2,0.6)
    #if((k==3) & (j==1)):a1.set_ylim(-0.12,0.4)
    a1.grid(True)
    a1.plot(time1, RFo, '-', color="grey",label='Observed')
    #a1.plot(time2[:626], RFp[:626], 'b-', label='Predicted')
    a1.plot(time2, RFp, 'b-', label='Predicted')
    a1.legend()

    if(mtype ==0):d, beta, vpvs, qa, qb, pv, pd = v2mod(velmod)  # Convert velocity model from Voronoi format to plot format
    if(mtype ==1):pv, pd = l2mod(velmod)  # Convert velocity model from Layer format to plot format
    if(mtype ==2):pv, pd = d2mod(velmod)  # Convert velocity model from Layer format to plot format
    
    a0.set_title(" Velocity model")                   # Plot velocity model with Receiver function
    a0.set_xlabel('Vs (km/s)')
    a0.set_ylabel('Depth (km)')
    a0.plot(pv,pd,'g-')
    a0.set_xlim(vmin,vmax)
    a0.set_ylim(dmin,dmax)
    a0.invert_yaxis()
    if(mtype==0): a0.plot(velmod[:,1],velmod[:,0],ls="",marker="o",markerfacecolor='none', markeredgecolor='k')

    #plt.tight_layout()
    #plt.savefig('RF_mod_plot.pdf',format='PDF')
    if(filename is not None): plt.savefig(filename)
    plt.show()
    return pv,pd
##################################################################################

__all__ = ["rfcalc","plotRFm","plot_RFs"]

