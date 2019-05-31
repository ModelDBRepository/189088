from __future__ import division

import xppy
import os.path
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

path_old = os.path.expanduser("./booth_bose.ode")
path_cont = os.path.expanduser("./booth_bose_cont.ode")

def plot_firing2(Isom0, gnmda, gc, time,gs1,gs2,title, VNa=60):

    xppy.changeOde([["par", "Isom0", Isom0],
                    ["par", "gnmda", gnmda],
                    ["par", "gc", gc],
                    ["par", "VNa", VNa],
                    ["@", "total", time]],
                    ode_file=path_old)
                
    xppy.changeOde([["par", "Isom0", Isom0],
                    ["par", "gnmda", gnmda],
                    ["par", "gc", gc],
                    ["par", "VNa", VNa],
                    ["@", "total", time]],
                    ode_file=path_cont)

    
    results_old = xppy.run(ode_file=path_old, set_file=os.path.abspath("./Data/data_old.dat"))
    results_cont = xppy.run(ode_file=path_cont, set_file=os.path.abspath("./Data/data_cont.dat"))
    
    t_old = results_old['t']
    Vs_old = results_old['Vs']
    Ca_old = results_old['Cad']
    
    t_cont = results_cont['t']
    Vs_cont = results_cont['Vs']
    Ca_cont = results_cont['Cad']
    
    ax1 = plt.subplot(gs1)
    ax2 = plt.subplot(gs2)

    
    ax1.plot(t_old, Vs_old)
    axCa1 = ax1.twinx()
    axCa1.plot(t_old,Ca_old, 'r')
    ax1.set_xlim([0,time])
    ax1.set_title(title,fontsize=14)
    ax1.set_xticks([])
    ax1.set_yticks([-60,20])
    axCa1.set_yticks([0,round((max(Ca_old)/50))*50])
    
    ax2.plot(t_cont,Vs_cont)
    axCa2 = ax2.twinx()
    axCa2.plot(t_cont,Ca_cont, 'r')
    ax2.set_xlim([0,time])
    ax2.set_xticks([])
    ax2.set_yticks([-60,20])
    axCa2.set_yticks([0,round((max(Ca_cont)/50))*50])
    
    return ax1,ax2,axCa1,axCa2

def plot_fi(I,f1,f2,p1,p2,title,gs):
    
    axfi = plt.subplot(gs)
    axfi.semilogy(I, f1, '--b')
    axfi.semilogy(I, f2, '--r')
    axfi.set_title(title,fontsize=14)

    markers = ['^','o','s']

    for i in [1,2,3]:
        m1 = np.extract(p1==i,f1)
        I1 = np.extract(p1==i,I)
        axfi.semilogy(I1,m1,'b'+markers[i-1])
    
        m2 = np.extract(p2==i,f2)
        I2 = np.extract(p2==i,I)
        axfi.semilogy(I2,m2,'r'+markers[i-1])
    
    axfi.set_xlim([-0.5,3.25])
    axfi.set_xticks([])
    
    return axfi

def classify_trace(trace,vthresh=-50,durthresh=30):
    pattern = np.extract(trace>vthresh,np.arange(0,len(trace),1))
    pos = np.where(np.diff(pattern)>1)[0] + 1
    Iter = iter(np.split(pattern,pos))
    durs = []
    for line in Iter:
        durs.append(len(line))
    
    num_bursts = np.sum(np.array(durs)>=durthresh)
    num_spikes = np.sum(np.array(durs)<durthresh)
    pbursts = num_bursts/(num_spikes+num_bursts)
    
    if not durs:
        c=0
    elif pbursts <= 0.1:
        c = 1 #spiking
    elif pbursts >= 0.9:
        c = 2 #bursting
    else:
        c = 3 #mixed
            
    return c,num_spikes,num_bursts

xppy.changeOde([["par","VNa", 60],
                ["par","gnmda",0]],
                ode_file=path_old)
xppy.changeOde([["par","VNa", 60],
                ["par","gnmda",0]],
                ode_file=path_cont)

Isom0 = np.arange(-0.5,3.0,0.2)
Iden0 = np.arange(-0.5,3.0,0.2)
gc = [0,2.1,2.1,1e6,1.85]

fI_som_old = np.zeros([len(gc),len(Isom0)])
fI_som_cont = np.zeros([len(gc),len(Isom0)])

fI_den_old = np.zeros([len(gc),len(Iden0)])
fI_den_cont = np.zeros([len(gc),len(Iden0)])

pat_som_old = np.zeros([len(gc),len(Isom0)])
pat_som_cont = np.zeros([len(gc),len(Isom0)])

pat_den_old = np.zeros([len(gc),len(Iden0)])
pat_den_cont = np.zeros([len(gc),len(Iden0)])

time = 10000

for i in xrange(len(gc)):
    for j in xrange(len(Isom0)):
        isomj = Isom0[j]
        if gc[i] == 0:
            idenj = Iden0[j]
        elif i == 2:
            idenj = Iden0[j]
            isomj = 0
        else:
            idenj = 0
    
        xppy.changeOde([["par", "Isom0", isomj],
                        ["par", "Iden0", idenj],
                        ["par", "gc", gc[i]],
                        ["@", "total", time]],
                        ode_file=path_old)
                
        xppy.changeOde([["par", "Isom0", isomj],
                        ["par", "Iden0", idenj],
                        ["par", "gc", gc[i]],
                        ["@", "total", time]],
                        ode_file=path_cont)
        
        results_old = xppy.run(ode_file=path_old, set_file=os.path.abspath("./Data/data_old.dat"))
        results_cont = xppy.run(ode_file=path_cont, set_file=os.path.abspath("./Data/data_cont.dat"))
    
        t_old = results_old['t']
        Vs_old = results_old['Vs']
        Vd_old = results_old['Vd']
    
        t_cont = results_cont['t']
        Vs_cont = results_cont['Vs']
        Vd_cont = results_cont['Vd']
    
        sptimes_som_old = []
        sptimes_som_cont = []
        sptimes_den_old = []
        sptimes_den_cont = []
        dt = 0.3
    
        for k in xrange(len(t_old)):
            if not(sptimes_som_old):
                cond_som_old = not(sptimes_som_old)
            else:
                cond_som_old = (k*dt - sptimes_som_old[-1]) > 1.5
            
            if not(sptimes_som_cont):
                cond_som_cont = not(sptimes_som_cont)
            else:
                cond_som_cont = (k*dt - sptimes_som_cont[-1]) > 1.5
                
            
            if cond_som_old & (Vs_old[k] > 10):
                sptimes_som_old.append(k*dt)
            
            if cond_som_cont & (Vs_cont[k] > 10):
                sptimes_som_cont.append(k*dt)
                
            if not(sptimes_den_old):
                cond_den_old = not(sptimes_den_old)
            else:
                cond_den_old = (k*dt - sptimes_den_old[-1]) > 10

            if not(sptimes_den_cont):
                cond_den_cont = not(sptimes_den_cont)
            else:
                cond_den_cont = (k*dt - sptimes_den_cont[-1]) > 10
            
            
            if cond_den_old & (Vd_old[k] > -20):
                sptimes_den_old.append(k*dt)
                
            if cond_den_cont & (Vd_cont[k] > -20):
                sptimes_den_cont.append(k*dt)
                    
    
        fI_som_old[i,j] = len(sptimes_som_old)*1000.0/time
        fI_som_cont[i,j] = len(sptimes_som_cont)*1000.0/time
        fI_den_old[i,j] = len(sptimes_den_old)*1000.0/time
        fI_den_cont[i,j] = len(sptimes_den_cont)*1000.0/time

        pat_som_old[i,j],s,b = classify_trace(Vs_old[500/0.3:],vthresh=-50,durthresh=30)
        pat_som_cont[i,j],s,b = classify_trace(Vs_cont[500/0.3:],vthresh=-50,durthresh=30)
        pat_den_old[i,j],s,b = classify_trace(Vd_old[500/0.3:],vthresh=-50,durthresh=35)
        pat_den_cont[i,j],s,b = classify_trace(Vd_cont[500/0.3:],vthresh=-50,durthresh=35)


traub_path = os.path.expanduser("./Data/traub.csv")
traub = np.genfromtxt(traub_path,delimiter=',')
traub_I = traub[:,0]
traub_Hz = traub[:,1]

bfreq = np.extract(traub_Hz <10,traub_Hz)
bI = np.extract(traub_Hz <10,traub_I)
sfreq = np.extract(traub_Hz >10,traub_Hz)
sI = np.extract(traub_Hz >10,traub_I)
        
        
fig = plt.figure()
outer = gridspec.GridSpec(nrows=1,ncols =2, width_ratios = [3,2],wspace=0.3)
gs_firing = gridspec.GridSpecFromSubplotSpec(nrows=5,ncols=1,subplot_spec=outer[0],hspace=0.3)

gs_firing1 = gridspec.GridSpecFromSubplotSpec(nrows=2,ncols=1,subplot_spec=gs_firing[0])
gs_firing2 = gridspec.GridSpecFromSubplotSpec(nrows=2,ncols=1,subplot_spec=gs_firing[1])
gs_firing3 = gridspec.GridSpecFromSubplotSpec(nrows=2,ncols=1,subplot_spec=gs_firing[2])
gs_firing4 = gridspec.GridSpecFromSubplotSpec(nrows=2,ncols=1,subplot_spec=gs_firing[3])
gs_firing5 = gridspec.GridSpecFromSubplotSpec(nrows=2,ncols=1,subplot_spec=gs_firing[4])

time = 1000
ax1,ax2,axCa1,axCa2= plot_firing2(0.75, 0, 2.1, time,gs_firing1[0],gs_firing1[1],'Somatic Input')
ax3,ax4,axCa3,axCa4= plot_firing2(-0.5, 1.25, 2.1, time,gs_firing2[0],gs_firing2[1], 'Dendritic Input')
ax5,ax6,axCa5,axCa6= plot_firing2(2.5, 0, 2.1, time,gs_firing3[0],gs_firing3[1], 'Strong Somatic Input')
ax7,ax8,axCa7,axCa8= plot_firing2(2.5, 0, 10.5, time,gs_firing4[0],gs_firing4[1], 'Strong Somatic Input, tight coupling')
ax9,ax10,axCa9,axCa10= plot_firing2(-0.5, 1.75, 1.425, time,gs_firing5[0],gs_firing5[1], 'Dendritic Input, weak coupling',VNa=55)

ax10.set_xticks(np.linspace(0,time,5))

gs_fi = gridspec.GridSpecFromSubplotSpec(nrows=3,ncols=2, subplot_spec=outer[1],wspace=0.15,hspace=0.15)

ax11 = plot_fi(Isom0,fI_som_old[0,:],fI_som_cont[0,:],pat_som_old[0,:],pat_som_cont[0,:],'Soma Only',gs_fi[0,0])
ax12 = plot_fi(Iden0,fI_den_old[0,:],fI_den_cont[0,:],pat_den_old[0,:],pat_den_cont[0,:],'Dendrite Only',gs_fi[0,1])
ax13 = plot_fi(Isom0,fI_som_old[1,:],fI_som_cont[1,:],pat_som_old[1,:],pat_som_cont[1,:],'Somatic input',gs_fi[1,0])
ax14 = plot_fi(Iden0,fI_den_old[2,:],fI_den_cont[2,:],pat_den_old[2,:],pat_den_cont[2,:],'Dendritic input',gs_fi[1,1])
ax15 = plot_fi(Isom0,fI_som_old[3,:],fI_som_cont[3,:],pat_som_old[3,:],pat_som_cont[3,:],'Infinite coupling',gs_fi[2,0])
ax16 = plot_fi(Isom0,fI_som_old[4,:],fI_som_cont[4,:],pat_som_old[4,:],pat_som_cont[4,:],'Traub',gs_fi[2,1])

ax16.semilogy(traub_I,traub_Hz,'--k', label='Traub')
ax16.semilogy(bI,bfreq,'ko', label='bursting')
ax16.semilogy(sI,sfreq,'k^', label='spiking')

ax12.yaxis.tick_right()
ax14.yaxis.tick_right()
ax16.yaxis.tick_right()

ax15.set_xticks([0,1,2,3])
ax16.set_xticks([0,1,2,3])

fig.set_size_inches([9,7.5])

fig.text(0.3,0.01, 'time (ms)',ha='center',va='center')
fig.text(0.8,0.01, 'I $\mu$A/cm$^2$',ha='center',va='center')
fig.text(0.01,0.5, 'Membrane Potential (mV)',ha='center',va='center',rotation='vertical')
fig.text(0.5775,0.5, 'Calcium (au)',ha='center',va='center',rotation=270)
fig.text(0.99,0.5,'Frequency (Hz)',ha='center',va='center',rotation=270)
outer.tight_layout(fig)

fig.savefig('./fig2.png',format='png')

fig.show()
