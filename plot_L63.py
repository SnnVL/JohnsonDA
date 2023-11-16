
import numpy as np 
import matplotlib.pyplot as plt
import pickle
from scipy import stats
from cmcrameri import cm

meths = [ \
    'Gaussian', \
    'Johnson ' \
]
n_meths = len(meths)
n_plot = n_meths

# State to compare (either 'a' for analysis or 'b' for background)
comp_state = 'a'

# File name
tf='bus'
with open('./data/Johnson_L63_'+tf+'.pkl','rb') as f:
    info = pickle.load(f)
    MAE_A = pickle.load(f)
    MAE_B = pickle.load(f)
    FAILS = pickle.load(f)
    MAE_A = np.flip(MAE_A,axis=-1)
    MAE_B = np.flip(MAE_B,axis=-1)
print(info)
print("Total number of fails = "+str(np.sum(FAILS)))

p_vec = info['period_obs']
o_vec = np.arange(info['var_obs'].size)

r_A = np.empty((p_vec.size,o_vec.size,n_meths))
r_B = np.empty((p_vec.size,o_vec.size,n_meths))
min_A = np.empty((p_vec.size,o_vec.size,n_meths))
max_A = np.empty((p_vec.size,o_vec.size,n_meths))

# Savefiles
fls = np.empty(n_plot,dtype='object')
for n_c in range(n_plot):
    fls[n_c] = open("./values/L63_"+meths[n_c].strip()+"_"+tf+".txt",'w')

# Get mean values of the mean absolute errors
for ii in range(p_vec.size):
    for jj in range(o_vec.size):

        r_A[ii,jj,:] = np.mean(MAE_A[ii,jj,:,:], axis = 0)
        min_A[ii,jj,:] = np.min(MAE_A[ii,jj,:,:], axis = 0)
        max_A[ii,jj,:] = np.max(MAE_A[ii,jj,:,:], axis = 0)
        r_B[ii,jj,:] = np.mean(MAE_B[ii,jj,:,:], axis = 0)
        
        for n_c in range(n_plot):
            if n_c==0:
                fls[n_c].write(str(o_vec[jj])+" "+str(p_vec[ii])+" "+str(r_A[ii,jj,n_c])+"\n")
            else:
                fls[n_c].write(str(o_vec[jj])+" "+str(p_vec[ii])+" "+str(r_A[ii,jj,n_c]/r_A[ii,jj,0])+"\n")
            
    for n_c in range(n_plot):
        fls[n_c].write("\n")

if comp_state == 'a':
    r = r_A
elif comp_state == 'b':
    r = r_B

# Significance testing
alpha = 1-0.99**(1/p_vec.size/o_vec.size)
alpha=1e-4
dof = 250*50*30-1
f_min = stats.f.ppf(alpha/2  ,dof,dof)
f_max = stats.f.ppf(1-alpha/2,dof,dof)
print("[f_min, f_max] = "+str(f_min)+","+str(f_max))
print("alpha   = "+str(alpha))
print("alpha^n = "+str((1-alpha)**(p_vec.size*o_vec.size)))

# Create figure
fig, ax = plt.subplots(1,n_plot,figsize=(n_plot*(o_vec.size-2),(p_vec.size)))
ii = 0
for n_c in range(0,n_plot):
    if n_c == 0:
        print(np.min(r[:,:,n_c]),np.max(r[:,:,n_c]))
        cq = ax[ii].pcolormesh(o_vec,p_vec,r[:,:,n_c], cmap=cm.batlow)
        cb = plt.colorbar(cq,ax=ax[ii], orientation = "horizontal")
    else:
        print(np.min(r[:,:,n_c]/r[:,:,0]),np.max(r[:,:,n_c]/r[:,:,0]))
        vmin, vmax = 0.8, 1.2
        cq = ax[ii].pcolormesh(o_vec,p_vec,r[:,:,n_c]/r[:,:,0], cmap=cm.vik,vmin = vmin,vmax=vmax)
        cb = plt.colorbar(cq,ax=ax[ii], orientation = "horizontal", \
            label = "MAE/MAE$_{\\rm g}$", extend='both')

        fls[n_c].write("\n\n\n")
        fls[n_c].write("v_max,p_max\n")
        for i_p in range(p_vec.size):
            for i_v in range(o_vec.size):
                if r[i_p,i_v,n_c]/r[i_p,i_v,0]>f_max:
                    ax[ii].scatter(o_vec[i_v],p_vec[i_p],s=100,marker='x',facecolors='k')
                    fls[n_c].write(str(o_vec[i_v])+" "+str(p_vec[i_p])+"\n")
        fls[n_c].write("\n\n\n")
        fls[n_c].write("v_min,p_min\n")
        for i_p in range(p_vec.size):
            for i_v in range(o_vec.size):
                if r[i_p,i_v,n_c]/r[i_p,i_v,0]<f_min:
                    ax[ii].scatter(o_vec[i_v],p_vec[i_p],s=100,marker='o',facecolors='none',edgecolors='k')
                    fls[n_c].write(str(o_vec[i_v])+" "+str(p_vec[i_p])+"\n")

    ax[ii].set_xticks(o_vec,np.round((info['var_obs']),2))
    ax[ii].set_title(meths[n_c])
    ax[ii].set_xlabel("Observation variance")
    ax[ii].set_ylabel("Observation period")
    ax[ii].set_yticks(p_vec)
    ii += 1

for n_c in range(n_plot):
    fls[n_c].close()
fig.savefig("plot_L63_"+tf+".png")


# Histograms
print("")
print("")
# ie=0
n_p = 3
fig, ax = plt.subplots((n_p),1,figsize=(7, 5*(n_p)))
for ie in range(n_p):
    for n_c in range(n_plot):
        ax[ie].plot(o_vec,r[ie,:,n_c],label=meths[n_c])
        ax[ie].plot(o_vec,min_A[ie,:,n_c],'--')
        ax[ie].plot(o_vec,max_A[ie,:,n_c],'--')
        print("")
        print("p = "+str(p_vec[ie])+", "+meths[n_c]+", "+tf)
        for iv in o_vec:
            # print(iv,r[ie,iv,n_c],min_A[ie,iv,n_c],max_A[ie,iv,n_c])
            print(iv,r[ie,iv,n_c])
    ax[ie].set_xticks(o_vec,np.round(info['var_obs'],2))
    ax[ie].legend()
    ax[ie].set_title("$period = $"+str(p_vec[ie]))
fig.savefig("histograms_"+tf+".png")
