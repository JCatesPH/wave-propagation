#%%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#%%
nprocs = 8
nx = 2400
ny = 101
idx1 = 320
idx2 = 2080


#%% Read in output
df = pd.read_csv('data/X.csv', header=None)
X = df.values[:,0]

df = pd.read_csv('data/Y.csv', header=None)
Y = df.values[:,0]

#%%
def filterprocs(arr, nprocs):
    procs = ['Process[{}]'.format(n) for n in range(nprocs)]
    for p in procs:
        idx = np.where(arr == p)
        arr = np.delete(arr, idx[0])
    return arr


# %%
df = pd.read_csv('data/epszz.csv', header=None, skiprows=2, delimiter=',')
epszz = df.values[:,0].astype(str)
epszz = np.char.replace(epszz,'i','j')
epszz = np.char.replace(epszz,' ','')
epszz = filterprocs(epszz, nprocs)
epszz = np.reshape(epszz.astype(complex), (ny,nx))

#%%
plt.clf()
plt.contourf(X, Y, np.abs(epszz), 50)
plt.title(r'Magnitude of relative permittivity')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$y$ [m]')
plt.colorbar()
#plt.show()
plt.savefig('figs/epszz.png')

# %%
df = pd.read_csv('data/q.csv', header=None, skiprows=2, delimiter=',')
Q = df.values[:,0].astype(str)
Q = np.char.replace(Q,'i','j')
Q = np.char.replace(Q,' ','')
Q = filterprocs(Q, nprocs)
Q = np.reshape(Q.astype(complex), (ny,nx))

#%%
plt.clf()
plt.contourf(X, Y, np.abs(Q))
plt.title(r'Masking matrix, $Q$')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$y$ [m]')
#plt.show()
plt.savefig('figs/Q.png')

# %%
df = pd.read_csv('data/src.csv', header=None, skiprows=2, delimiter=',')
src = df.values[:,0].astype(str)
src = np.char.replace(src,'i','j')
src = np.char.replace(src,' ','')
src = filterprocs(src, nprocs)
src = np.reshape(src.astype(complex), (ny,nx))

#%%
plt.clf()
plt.contourf(X, Y, np.real(src), 50)
plt.title(r'Source, Re[$f_{src}(x,y)$]')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$y$ [m]')
cbar = plt.colorbar()
cbar.formatter.set_powerlimits((0, 0))
#plt.show()
plt.savefig('figs/fsrc.png')

# %%
df = pd.read_csv('data/b.csv', header=None, skiprows=2, delimiter=',')
b = df.values[:,0].astype(str)
b = np.char.replace(b,'i','j')
b = np.char.replace(b,' ','')
b = filterprocs(b, nprocs)
b = np.reshape(b.astype(complex), (ny,nx))

#%%
plt.clf()
plt.contourf(X, Y, np.abs(b), 50)
plt.title(r'Reshaped source vector, $|\mathbf{b}|$')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$y$ [m]')
cbar = plt.colorbar()
cbar.formatter.set_powerlimits((0, 0))
#plt.show()
plt.savefig('figs/b.png')


#%%
df = pd.read_csv('data/Ez_om.csv', header=None, skiprows=2, delimiter=',')
Ez1 = df.values[:,0].astype(str)
Ez1 = np.char.replace(Ez1,'i','j')
Ez1 = np.char.replace(Ez1,' ','')
Ez1 = filterprocs(Ez1, nprocs)
Ez1 = np.reshape(Ez1.astype(complex), (ny,nx))

#%%
plt.clf()
plt.contourf(X, Y, np.real(Ez1), 50)
plt.title(r'Solution, $E_z(\omega)$')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$y$ [m]')

cbar = plt.colorbar()
cbar.formatter.set_powerlimits((0, 0))

#plt.show()
plt.savefig('figs/Ez_om.png')

#%%
df = pd.read_csv('data/Ez_2om.csv', header=None, skiprows=2, delimiter=',')
Ez2 = df.values[:,0].astype(str)
Ez2 = np.char.replace(Ez2,'i','j')
Ez2 = np.char.replace(Ez2,' ','')
Ez2 = filterprocs(Ez2, nprocs)
Ez2 = np.reshape(Ez2.astype(complex), (ny,nx))

#%%
plt.clf()
plt.contourf(X, Y, np.real(Ez2), 50)
plt.title(r'Solution, $E_z(2\omega)$')
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$y$ [m]')
cbar = plt.colorbar()
cbar.formatter.set_powerlimits((0, 0))
#plt.show()
plt.savefig('figs/Ez_2om.png')



# %% Calculate conversion efficiency
n_om = np.abs(Ez1) / (2.156*np.abs(src))
n_2om = np.abs(Ez2) / (2.156*np.abs(src))

#%%
plt.clf()
fig, ax = plt.subplots(2)
#fig.suptitle('Conversion efficiency')
ax[0].plot(X, n_om[int(n_om.shape[0]/2),:])
ax[0].set_title('Conversion efficiency for first harmonic, $\eta_{\omega}$')

ax[1].plot(X, n_2om[int(n_om.shape[0]/2),:])
ax[1].set_title('Conversion efficiency for second harmonic, $\eta_{2\omega}$')

ax[0].ticklabel_format(axis='both', style='sci', scilimits=(0,0))
ax[1].ticklabel_format(axis='both', style='sci', scilimits=(0,0))

ax[1].set_xlabel(r'Sample length [$\mu$m]')

plt.tight_layout()
#plt.show()
plt.savefig('figs/conv-eff.png')

# %%
plt.clf()
fig, ax = plt.subplots(2)
#fig.suptitle('Conversion efficiency')
ax[0].plot(X[idx1:idx2]*1e6, n_om[int(n_om.shape[0]/2), idx1:idx2])
ax[0].set_title('Conversion efficiency for first harmonic, $\eta_{\omega}$')

ax[1].plot(X[idx1:idx2]*1e6, n_2om[int(n_om.shape[0]/2), idx1:idx2])
ax[1].set_title('Conversion efficiency for second harmonic, $\eta_{2\omega}$')

ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0,0))

ax[1].set_xlabel(r'Sample length [$\mu$m]')

plt.tight_layout()
#plt.show()
plt.savefig('figs/conv-eff_inSHG.png')

# %%
