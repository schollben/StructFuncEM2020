# Python script used to simulate spine inputs from EM data collected in 
# Scholl et al. 2020. A 'data' structure contains anatomical information 
# (spine head volume, PSD area, neck length, and somatic distance) which 
# is used to describe each spine. Then a single presynaptic spike is simulated 
# and the resulting depolarizations are recorded.
#
# Uses the NEURON-Python module and compiling of mechanisms
# Make sure to cd to directory where .mod files are located
# NMDA can be included (see commented lines below)
#
# B. Scholl 2020


import os

#cd to working directory with .mod files
os.chdir('[insert path here]')

import numpy
from matplotlib import pyplot
from neuron import h, gui
from pylab import *
import scipy.io as sio

# helpful functions:
# shape_window = h.PlotShape() - see shape of everything
den_l = 400 #this will make a 400 micron long dendrite so no effects of the distal tip
soma = h.Section(name='soma')
dend = [h.Section(name='dend[%d]' % i) for i in xrange(den_l)]

for i in xrange(den_l):
	dend[i].L = 1
	dend[i].diam = 1
	dend[i].insert('pas')
	dend[i].nseg = 1    
    
dend[0].connect(soma(0.5),0)

for i in xrange(den_l-1):
	dend[i+1].connect(dend[i](1),0)

# Geometry
soma.L = soma.diam = 13 # microns

h.define_shape() # Translate into 3D points.

# Insert active Hodgkin-Huxley current in the soma
soma.insert('hh')
soma.nseg = 21

# Insert active Hodgkin-Huxley current in the soma
for seg in soma:
	seg.hh.gnabar = 0
	seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
	seg.hh.gl = 0.0003    # Leak conductance in S/cm2
	seg.hh.el = -50     # Reversal potential in mV

# Biophysics
for sec in h.allsec():
	sec.Ra = 105      # Axial resistance in Ohm * cm
	sec.cm = 1       # Membrane capacitance in micro Farads / cm^2
    
    
#specify spine anatomical parameters or load data 
#below is an example of parameters for 2 different input2
 
data = [[1,1,1,100], [1,1,4,100]]

spineV = numpy.zeros(len(data))
dendV = numpy.zeros(len(data))
somaV = numpy.zeros(len(data))

for i in xrange(len(data)):
    
    print(i)
    
    spine_head_volume = data[i][0] # volume in um^3, going to assume a cylindar with length 1
    spine_PSD = data[i][1] # area in um^2
    spine_neck_length = data[i][3] # length in um
    spine_soma_distance = data[i][4] # length in um
    
    if spine_neck_length < 0.001:
        spine_neck_length = 0.001
        
    if ~isnan(data[i][0]):
        
        #synapse assumptions: 
        #AMPA grows linearly w PSD length (Takumi 1999, SCC): 
        #0.87 particles per 0.1 micron PSD
        #take lower limit: 1 particle ~ 1 receptor 
        #assume measured PSD is circular
        #single AMPA-R ~10pS
        #AMPA/NMDA ratio linearly grows with radius (Takumi 1999)
        #NMDA particles (Takumi 1999, SCC): 2.3 particles / 100nm PSD
        #NMDA g ~ scaled AMPA
        
        mininumArea = (0.09 * 0.09 * pi) #180 nm diameter minimum (Takumi 1999)- assume circular
        
        mSpineDiameter = 2 * sqrt(spine_PSD / pi) # um
		
        AMPAgmax = mSpineDiameter * (0.87 / 0.1) * 50 #convert um^2 to receptor pS
        
        AMPAgmax = ((AMPAgmax / 1000) / 1000) #convert pS to uS
        
#	     #uncomment for including of NMDA channels
#        #only issue here is that NMDA-Rs are more random
#        NMDAgmax = (2.3 / 0.87) * AMPAgmax # scale AMAP gmax based on density
        
        spine_head = h.Section(name='spine_head')
        spine_head.diam = 2 * sqrt(spine_head_volume /pi)
        spine_head.L = 1
        spine_head.Ra = 250
        spine_head.cm = 1
        spine_head.insert('pas') #active conductances?
        
        spine_neck = h.Section(name='spine_neck') 
        spine_neck.L = spine_neck_length
        spine_neck.diam = 0.200
        spine_neck.Ra = 250
        spine_neck.cm = 1
        spine_neck.insert('pas')
        
        spineLoc= int(spine_soma_distance) #den_l/2
        spine_neck.connect(dend[spineLoc](0.5),0)
        spine_head.connect(spine_neck(1),0)
        
        h.define_shape() # Translate into 3D points.
        
        syn = h.AlphaSynapse(spine_head(0.5))
        syn.gmax = AMPAgmax #uS
        spike_onset = 20 #ms
        syn.onset = spike_onset #ms
        syn.tau = 1.5 #ms off kinetic
        syn.e = 0
        
#	     #uncomment for including of NMDA channels
#        #insert additional channels to synpase (NMDA for now, add Ca?)
#        NMDA_syn = h.Exp2SynNMDA(spine_head(0.5)) # conductance is set in uS
#        NMDA_syn.tau1 = 1
#        NMDA_syn.tau2 = 25
#        NC = h.NetCon(None, NMDA_syn, 0, 0, NMDAgmax)	# connection; from, to, threshold, delay, weight
#        def synStim():
#            NC.event(spike_onset)
        
        #record
        v_vecSoma = h.Vector()            
        v_vecDend = h.Vector()
        v_vecSpine = h.Vector()                 
        t_vec = h.Vector()             
        v_vecSoma.record(soma(0.5)._ref_v)
        v_vecDend.record(dend[spineLoc](0.5)._ref_v)
        v_vecSpine.record(spine_head(0.5)._ref_v)
        t_vec.record(h._ref_t)
        
        #init and run
        h.v_init = -67.5
        h.tstop = 100
#	    #uncomment for including of NMDA channels
#       fih = h.FInitializeHandler(synStim)
        h.run()
        
        t = numpy.array(t_vec)
        v_Soma = numpy.array(v_vecSoma)
        v_Dend = numpy.array(v_vecDend)
        v_Spine = numpy.array(v_vecSpine)
        
        
        spineV[i] = numpy.amax(v_Spine) - v_Spine[400]
        dendV[i] = numpy.amax(v_Dend) - v_Dend[400]
        somaV[i] = numpy.amax(v_Soma) - v_Soma[400]

    else:
        
        spineV[i] = NaN
        dendV[i] = NaN
        somaV[i] = NaN      

pyplot.figure(figsize=(8,4))
pyplot.plot(t,v_Spine, color='red')
pyplot.plot(t,v_Dend, color='blue')
pyplot.plot(t,v_Soma, color='black')
pyplot.xlabel('time (ms)')
pyplot.ylabel('mV')
pyplot.ylim(-70,0)
pyplot.show(block=False)

pyplot.figure(figsize=(8,4))
pyplot.plot(spineV,somaV,'.')



