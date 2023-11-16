""" 
Run a nongaussian 3DVAR twin experiment with Lorenz-63 
"""

# Basic modules
import numpy as np 
import pickle

# Lorenz-63 and DA methods
import mod_JohnsonDA as da

# Options for DA run
n_runs = 50                           # Number of DA runs
n_wind = 250                          # Number of DA windows for each run
seed = None                           # Seed for random number generator, can be set to None
SV_init = np.array([-5.0,-6.0,22.0])  # Initial values of the Lorenz-63 run

# Model definition
N = 3
def model(t, x0):
    y = da.sol_L63(t, x0, np.array([10.0,28.0,8.0/3.0]))
    return y

# Observation operator
def obs_h(x):
    return np.copy(x)
def phi(t):
    return np.copy(t)
def dphi(t):
    return np.eye(N)

# mean absolute error
def mae(SV_true,SV_DA,period_DA = 1):
    return np.nanmean(np.abs(SV_true[:,::period_DA] - SV_DA))

# Make sure initial conditions are consistent with model
dt = 0.01
t_prep = np.arange(0.0,1000*dt,dt)
SV_prep= model(t_prep,SV_init)
SV_init_clim = SV_prep[:,-1]

# Climatology
t = np.arange(0.0,100000*dt,dt)
SV_clim = model(t,SV_init_clim)
print(SV_init_clim)

min_SV = np.min(SV_clim, axis=1)-2.
max_SV = np.max(SV_clim, axis=1)+2.
mean_SV = np.mean(SV_clim, axis=1)
std_SV = np.std(SV_clim, axis=1)

# Define necessary functions
# str_transform = ['bounded','bounded','bounded']
# xi  = [min_SV[0],min_SV[1],min_SV[2]]
# lam = [max_SV[0]-min_SV[0],max_SV[1]-min_SV[1],max_SV[2]-min_SV[2]]
# str_transform = ['bounded','unbounded','bounded']
# xi  = [min_SV[0],mean_SV[1],min_SV[2]]
# lam = [max_SV[0]-min_SV[0],std_SV[1],max_SV[2]-min_SV[2]]
# str_transform = ['bounded','unbounded','unbounded']
# xi  = [min_SV[0],mean_SV[1],mean_SV[2]]
# lam = [max_SV[0]-min_SV[0],std_SV[1],std_SV[2]]
# str_transform = ['unbounded','unbounded','unbounded']
# xi  = [mean_SV[0],mean_SV[1],mean_SV[2]]
# lam = [std_SV[0],std_SV[1],std_SV[2]]
str_transform = ['bounded','unbounded','semibounded']
xi  = [min_SV[0],mean_SV[1],min_SV[2]]
lam = [max_SV[0]-min_SV[0],std_SV[1],std_SV[2]]

# print(str_transform)
# print(xi)
# print(lam)

t_fun_SV = lambda x: da.t_fun_L63(x,xi,lam,str_transform,N)
inv_t_fun_SV = lambda t: da.inv_t_fun_L63(t,xi,lam,str_transform,N)
log_s_fun_SV = lambda t: da.log_s_fun_L63(t,str_transform,N)
f_fun_SV = lambda t: da.f_fun_L63(t,str_transform,N)

t_fun_obs = t_fun_SV
log_s_fun_obs = log_s_fun_SV
f_fun_obs = f_fun_SV

t_SV_clim = t_fun_SV(SV_clim)

# Define necessary functions
str_transform_gauss = ['normal','normal','normal']
xi_gauss = [0.,0.,0.]
lam_gauss = [1.,1.,1.]

t_fun_SV_gauss = lambda x: da.t_fun_L63(x,xi_gauss,lam_gauss,str_transform_gauss,N)
inv_t_fun_SV_gauss = lambda t: da.inv_t_fun_L63(t,xi_gauss,lam_gauss,str_transform_gauss,N)
log_s_fun_SV_gauss = lambda t: da.log_s_fun_L63(t,str_transform_gauss,N)
f_fun_SV_gauss = lambda t: da.f_fun_L63(t,str_transform_gauss,N)

t_fun_obs_gauss = t_fun_SV_gauss
log_s_fun_obs_gauss = log_s_fun_SV_gauss
f_fun_obs_gauss = f_fun_SV_gauss


#################################################################################
#################################################################################
#################################################################################

def one_run(jj):
    global period_obs
    global var_obs

    # Randomize initial conditions
    rng = np.random.default_rng(seed)
    init_guess = SV_init_clim + rng.standard_normal(N)

    # Make sure initial conditions are consistent with model
    t_prep = np.arange(0.0,1000*dt,dt)
    SV_prep= model(t_prep,init_guess)
    SV_init = SV_prep[:,-1]

    # Nature run
    t_max = n_wind*period_obs*dt        # End time of the Lorenz-63 run
    t = np.arange(0.0, t_max, dt)  # Evaluation time
    SV = model(t, SV_init)

    # Generate observations
    t_obs, y, R = da.gen_obs(t, SV, period_obs, obs_h, \
                            var_obs, t_fun_SV, inv_t_fun_SV, seed = seed)
    
    # Generate B matrix
    B = np.cov(t_SV_clim[:,::period_obs])
    B_gauss = np.cov(SV_clim[:,::period_obs])

    # Johnson DA
    init_guess_da = SV_init + rng.standard_normal(N)
    x_a, x_b, _ = da.johnson_var3d(\
        init_guess_da, t_obs, period_obs, y, phi, dphi, B, R, model, \
        t_fun_SV, t_fun_obs, inv_t_fun_SV, \
        log_s_fun_SV, log_s_fun_obs, \
        f_fun_SV, f_fun_obs)

    # Gaussian DA
    x_a_gauss, x_b_gauss, _ = da.johnson_var3d(\
        init_guess_da, t_obs, period_obs, y, phi, dphi, B_gauss, R, model, \
        t_fun_SV_gauss, t_fun_obs_gauss, inv_t_fun_SV_gauss, \
        log_s_fun_SV_gauss, log_s_fun_obs_gauss, \
        f_fun_SV_gauss, f_fun_obs_gauss)

    # Mean absolute errors
    mae_a = mae(SV, x_a[:,:], period_DA=period_obs)
    mae_b = mae(SV, x_b[:,:-1])
    mae_a_gauss = mae(SV, x_a_gauss[:,:], period_DA=period_obs)
    mae_b_gauss = mae(SV, x_b_gauss[:,:-1])

    n_fails = int(np.sum(np.isnan(x_a))/N)
    n_fails_gauss = int(np.sum(np.isnan(x_a_gauss))/N)

    return mae_a, mae_b, mae_a_gauss, mae_b_gauss, n_fails, n_fails_gauss


# Main program
meths=['Johnson ','Gaussian']
if __name__ == '__main__':

    # Values of period and error to test
    p_vec = np.arange(20,200,30)
    o_vec = np.array([0.01,0.02,0.05,0.1,0.2,0.5])

    MAE_A = np.empty((p_vec.size,o_vec.size,n_runs,2))
    MAE_B = np.empty((p_vec.size,o_vec.size,n_runs,2))
    FAILS = np.empty((p_vec.size,o_vec.size,n_runs,2))

    # Temporary savefile in case program is interrupted
    tmpFile = open("./data/tmp_cL63.pkl",'wb')
    for ii, period_obs in enumerate(p_vec):
        for jj, var_obs in enumerate(o_vec):
            print(ii,jj)

            for kk in range(n_runs):
                MAE_A[ii,jj,kk,0], MAE_B[ii,jj,kk,0], \
                MAE_A[ii,jj,kk,1], MAE_B[ii,jj,kk,1], \
                FAILS[ii,jj,kk,0], FAILS[ii,jj,kk,1] = one_run(kk)

            print("Finished p = "+str(period_obs)+", s = "+str(var_obs))

            for iM in range(2):
                r_a = np.mean(MAE_A[ii,jj,:,iM])
                r_b = np.mean(MAE_B[ii,jj,:,iM])

                print(meths[iM]+ \
                    ': mean(r_a) = '+format(np.round(r_a,3),".3f")+ \
                    ', mean(r_b) = '+format(np.round(r_b,3),".3f")+ \
                    ', num fails = '+str(np.sum(FAILS[ii,jj,:,iM])))

    info = {
        "n_runs": n_runs,
        "n_wind": n_wind,
        "period_obs": p_vec,
        "var_obs": o_vec,
        "meths": meths,
        "transforms":str_transform,
        "xi":xi,
        "lam":lam,
        "seed": seed,
        "SV_init": SV_init_clim
    }
    with open('./data/Johnson_L63' \
        +"_"+str_transform[0][0] \
            +str_transform[1][0] \
            +str_transform[2][0] \
        +'.pkl','wb') as f:
        pickle.dump(info, f)
        pickle.dump(MAE_A, f)
        pickle.dump(MAE_B, f)
        pickle.dump(FAILS, f)
    print("Finished")