#%%
import os, sys, subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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
print(sys.argv)
if sys.argv[2] == 'true':
    rmopt = True
else:
    rmopt = False

DOUT = 200

#%% Read in output
print('Reading in grid.')
df = pd.read_csv('data/X.csv', header=None)
Z = df.values[:,0]


#%%
print('Reading in field array.')
df = pd.read_csv("data/Ex.csv", header=None, delimiter=',')
Ez = df.values


# %%
print('Plotting field.')
plt.clf()
fig = plt.figure(figsize=(6,6), dpi=200)
ax = fig.add_subplot(111)
line1, = ax.plot(Z, Ez[2], 'b-')
ax.set_ylim(min(0.9*Ez.min(),1.1*Ez.min()), 1.1*Ez.max())
ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$E_z$')
text1 = ax.text(0.75, 0.9, 
    r'max$(E_z)=${:.3e}'.format(Ez[2].max()),
    horizontalalignment='center',
    verticalalignment='center',
    transform=ax.transAxes)

#%%
for n in range(0, Ez.shape[0]):
    T = n*DOUT

    line1.set_ydata(Ez[n])
    ax.set_title(r"$E_x(z)$ at step {t}".format(t=T))
    text1.set_text(r'max$(E_z)=${:.3e}'.format(Ez[n].max()))
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.savefig('figs/output{:}.png'.format(n))


print('Making animation.')
makemp4(framerate=sys.argv[1], outname='test', rm=rmopt)

# %%
print("Reading in frequencies.")
df = pd.read_csv('data/omeg.csv', header=None)
omeg = df.values[:,0]

#%%
print('Reading in frequency domain array.')
df = pd.read_csv("data/Ew_re.csv", header=None, delimiter=',')
Ew_re = df.values
df = pd.read_csv("data/Ew_im.csv", header=None, delimiter=',')
Ew_im = df.values
Ew = np.power(Ew_re,2) + np.power(Ew_im,2)

#%%
print('Plotting frequency domain array.')
fig = plt.figure(figsize=(7,6), dpi=200)
ax = fig.add_subplot(111)
m = ax.contourf(omeg, Z, Ew)

ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
ax.set_xlabel(r'$\omega$')
ax.set_ylabel(r'$z$')
ax.set_title(r'Intensity, $|E_z(\omega)|^2$')
fig.colorbar(m)
plt.savefig('figs/Ew.png')

# %%
fig = plt.figure(figsize=(6,6), dpi=200)
ax = fig.add_subplot(111)
ax.plot(omeg, Ew[int(Ew.shape[0]/2),:])

ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
ax.set_xlabel(r'$\omega$')
ax.set_title(r'Intensity spectrum, $|E_z(\omega)|^2$, at $z={z0:.2g}$'.format(z0=Z[int(Ew.shape[0]/2)]))
plt.savefig('figs/Ew_midpoint.png')

# %%
iom = 1
fig = plt.figure(figsize=(6,6), dpi=200)
ax = fig.add_subplot(111)
ax.plot(Z, Ew[:,iom])

ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
ax.set_xlabel(r'$z$')
ax.set_title(r'Intensity spectrum, $|E_z(\omega)|^2$, at $\omega={omeg1:.2g}$'.format(omeg1=omeg[iom]))
plt.savefig('figs/Ew_omeg1.png')

# %%
