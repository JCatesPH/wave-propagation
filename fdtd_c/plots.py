#%%
import os, sys
import numpy as np
import glob
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import subprocess


#%%
def makemp4(framerate=2, outname=None, rm=False):
    os.chdir('figs')
    bash1 = './makemp4.sh -f {0}'.format(framerate)
    if outname is not None:
        bash1 += ' -o {0}'.format(outname)
    if rm is True:
        bash1 += ' -r'
    subprocess.run(bash1.split(),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=True)

    os.chdir('..')

#%%
if sys.argv[2] == 0:
    rmopt = False
else:
    rmopt = True

#%% Read in output
print('Reading in grid.')
df = pd.read_csv('data/X.csv', header=None)
X = np.transpose(df.values)

df = pd.read_csv('data/Y.csv', header=None)
Y = np.transpose(df.values)

#%%
print('Reading in field arrays.')
Ez_li = []
for fi in glob.glob('data/Ez*'):
    df = pd.read_csv(fi, header=None, delimiter=',')
    Ez = df.values
    Ez_li.append(Ez)



# %%
print('Printing slices.')
plt.clf()
fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(X[:,0], Ez_li[2][250,:], 'b-')
ax.set_ylim(-100.0, 100.0)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$E_z$')

for n in range(1,15):
    T = (n-1)*10

    line1.set_ydata(Ez_li[n-1][250,:])
    ax.set_title(r"$E_z(x,y={y0}$ at $T={t}$".format(y0=Y[250,0], t=T))
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.savefig('figs/output{:}.png'.format(n))

print('Making slice animation.')
makemp4(framerate=sys.argv[1], outname='xslice', rm=rmopt)


#%%
print('Printing contours.')
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
cax = make_axes_locatable(ax).append_axes("right", size="5%", pad="2%")

for n in range(1,15):
    T = (n-1)*10

    cax.clear()
    cf = ax.contourf(X[:,0], Y[:,0], Ez_li[n-1])
    ax.set_title(r"$E_z$ at $T={t}$".format(t=T))
    fig.colorbar(cf, cax=cax)
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.savefig('figs/output{:}.png'.format(n))

print('Making contours animation.')
makemp4(framerate=sys.argv[1], outname='contours', rm=rmopt)

print('Exiting python script.')
