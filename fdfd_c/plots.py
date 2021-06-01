#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%
nprocs = 8
nx = 1000
ny = 201

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
plt.contourf(X, Y, np.abs(epszz))
plt.title(r'Magnitude of relative permittivity, $|\varepsilon_{zz}|$')
plt.colorbar()
plt.show()
plt.savefig('figs/epszz.png')

# %%
df = pd.read_csv('data/q.csv', header=None, skiprows=2, delimiter=',')
Q = df.values[:,0].astype(str)
Q = np.char.replace(Q,'i','j')
Q = np.char.replace(Q,' ','')
Q = filterprocs(Q, nprocs)
Q = np.reshape(Q.astype(complex), (ny,nx))

#%%
plt.contourf(X, Y, np.abs(Q))
plt.title(r'Masking matrix, $Q$')
plt.colorbar()
plt.show()
plt.savefig('figs/Q.png')

# %%
df = pd.read_csv('data/src.csv', header=None, skiprows=2, delimiter=',')
src = df.values[:,0].astype(str)
src = np.char.replace(src,'i','j')
src = np.char.replace(src,' ','')
src = filterprocs(src, nprocs)
src = np.reshape(src.astype(complex), (ny,nx))

#%%
plt.contourf(X, Y, np.real(src))
plt.title(r'Source, Re[$f_{src}(x,y)$]')
plt.colorbar()
plt.show()
plt.savefig('figs/fsrc.png')

# %%
df = pd.read_csv('data/b.csv', header=None, skiprows=2, delimiter=',')
b = df.values[:,0].astype(str)
b = np.char.replace(b,'i','j')
b = np.char.replace(b,' ','')
b = filterprocs(b, nprocs)
b = np.reshape(b.astype(complex), (ny,nx))

#%%
plt.contourf(X, Y, np.abs(b))
plt.title(r'Reshaped source vector, $|\mathbf{b}|$')
plt.colorbar()
plt.show()
plt.savefig('figs/b.png')


#%%
df = pd.read_csv('data/Ez_om.csv', header=None, skiprows=2, delimiter=',')
Ez = df.values[:,0].astype(str)
Ez = np.char.replace(Ez,'i','j')
Ez = np.char.replace(Ez,' ','')
Ez = filterprocs(Ez, nprocs)
Ez = np.reshape(Ez.astype(complex), (ny,nx))

#%%
plt.contourf(X, Y, np.real(Ez), 50)
plt.title(r'Solution, $E_z(\omega)$')
plt.colorbar()
plt.show()
plt.savefig('figs/Ez_om.png')

#%%
df = pd.read_csv('data/Ez_2om.csv', header=None, skiprows=2, delimiter=',')
Ez = df.values[:,0].astype(str)
Ez = np.char.replace(Ez,'i','j')
Ez = np.char.replace(Ez,' ','')
Ez = filterprocs(Ez, nprocs)
Ez = np.reshape(Ez.astype(complex), (ny,nx))

#%%
plt.contourf(X, Y, np.real(Ez), 50)
plt.title(r'Solution, $E_z(2\omega)$')
plt.colorbar()
plt.show()
plt.savefig('figs/Ez_2om.png')
