import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy.stats import norm
import random

u = MDAnalysis.Universe('/ExtDrive/maya/Data/Simulation/pepg_bilayer/step7_1/bilayer.gro',
                        '/ExtDrive/maya/Data/Simulation/pepg_bilayer/step7_1/bilayer_trunc.xtc' )#Simulation time-series data of a 2:1 POPE/POPG lipid bilayer.

#Randomly sample 101 POPE and POPG lipids and store them as a list of resid numbers.

random_sample_POPE = list(random.sample(range(0,224), 101))
random_sample_POPG = list(random.sample(range(225, 336), 101))

str1_POPE = '  '.join(str(i) for i in random_sample_POPE)
str2_POPE = str1_POPE.split()

str1_POPG = '  '.join(str(i) for i in random_sample_POPG)
str2_POPG = str1_POPG.split()

l1_POPE =['resid '+ i for i in str2_POPE]
l2_POPE = ' or '.join(str(i) for i in l1_POPE)

l1_POPG =['resid '+ i for i in str2_POPG]
l2_POPG = ' or '.join(str(i) for i in l1_POPG)

#Grab the POPG and POPE atoms in respective lists.

POPG_lipids = u.select_atoms(l2_POPG)
POPE_lipids = u.select_atoms(l2_POPE)

#Set two counters to zero.

i = 0
j = 0

#Initialize the mean (sigma) and standard deviation (mu) arrays for POPE and POPG.

POPE_sigma_array = []
POPE_mu_array = []

POPG_mu_array = []
POPG_sigma_array = []

#Calculate the radius of gyration for each randomly sampled POPE and POPG lipid.
#Calculate the mean and standard deviation and store in above arrays.

while i < len(random_sample_POPE) and j < len(random_sample_POPG):

    POPE_Rgyr = []
    POPG_Rgyr = []

    POPE_lipid_group = POPE_lipids.residues[i]
    POPG_lipid_group = POPG_lipids.residues[j]

    i = i + 1
    j = j + 1

    for ts in u.trajectory:
        POPE_Rgyr.append((u.trajectory.time, POPE_lipid_group.radius_of_gyration()))
        POPG_Rgyr.append((u.trajectory.time, POPG_lipid_group.radius_of_gyration()))

    POPE_Rgyr = np.array(POPE_Rgyr)
    POPG_Rgyr = np.array(POPG_Rgyr)

    (mu1, sigma1) = norm.fit(POPE_Rgyr[:,1])
    (mu2, sigma2) = norm.fit(POPG_Rgyr[:,1])

    POPE_mu_array.append(mu1)
    POPE_sigma_array.append(sigma1)

    POPG_mu_array.append(mu2)
    POPG_sigma_array.append(sigma2)
    
#Plot histograms.

plt.subplot(2,1,1)

n_mu, bins_mu, patches_mu = plt.hist([POPE_mu_array, POPG_mu_array], bins=20, normed=True, histtype='bar', color=['red', 'blue'], label=['POPE','POPG'] )
plt.legend()
plt.title(r"$R_G$ ($\AA$) Histogram of Randomly Sampled Lipids")
plt.ylabel('Frequency')
plt.xlabel('$\mu$')

plt.subplot(2,1,2)

n_sigma, bins_sigma, patches_sigma = plt.hist([POPE_sigma_array, POPG_sigma_array], bins=20, normed=True, histtype='bar', color=['red', 'blue'], label=['POPE','POPG'] )
plt.legend()
plt.title(r"$R_G$ ($\AA$) Histogram of Randomly Sampled Lipids")
plt.ylabel('Frequency')
plt.xlabel('$\sigma$')

plt.tight_layout()
plt.show()
plt.savefig('SD_mu_Histogram.png')
plt.close()









