#%%
import numpy as np
import matplotlib.pyplot as plt
# %%
loaded = np.loadtxt('stdout_merged')
# %%
histo = np.zeros((65,2))
# %%
for l in loaded:
    histo[int(l[0])][0] += l[1]
    histo[int(l[0])][1] += l[2]
# %%
histo
# %%
plt.clf()
x = np.arange(33)
plt.yscale('linear')
plt.bar(x, histo[:33,0], width=0.7, color='r', label='full')
plt.savefig('full_linear.pdf')
# %%
plt.clf()
x = np.arange(33)
plt.yscale('log')
plt.bar(x, histo[:33,0], width=0.7, color='r', label='full')
plt.savefig('full_log.pdf')
# %%
plt.clf()
x = np.arange(33)
plt.yscale('linear')
plt.bar(x, histo[:33,1], width=0.7, color='b', label='simplified')
plt.savefig('simplified_linear.pdf')
# %%
plt.clf()
plt.yscale('log')
plt.bar(x, histo[:33,1], width=0.7, color='b', label='simplified')
plt.savefig('simplified_log.pdf')
