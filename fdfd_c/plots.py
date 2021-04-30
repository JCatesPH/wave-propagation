#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%% Read in output
df = pd.read_csv('data/X.csv', header=None)
X = df.values

df = pd.read_csv('data/Y.csv', header=None)
Y = df.values


# %%
df = pd.read_csv('data/epszz.csv', header=None, delimiter=' ')
epszz = df.values[:,0] + 1j*df.values[:,1]
epszz = np.reshape(epszz, (1000,200))

#%%
plt.contourf(X[:,0], Y[:,0], np.transpose(np.abs(epszz)))
plt.colorbar()

# %%
df = pd.read_csv('data/src.csv', header=None, delimiter=' ')
src = df.values[:,0] + 1j*df.values[:,1]
src = np.reshape(src, (1000,200))

#%%
plt.contourf(X[:,0], Y[:,0], np.transpose(np.real(src)))
plt.colorbar()