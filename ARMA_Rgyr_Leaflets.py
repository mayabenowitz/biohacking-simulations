import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
from spectrum import arma_estimate, arma, arma2psd
from MDAnalysis.analysis.leaflet import LeafletFinder

u = MDAnalysis.Universe('/home/jarod/Documents/Simulation/pepg_bilayer/step7_1/bilayer.gro',
                        '/home/jarod/Documents/Simulation/pepg_bilayer/step7_1/bilayer_trunc.xtc' )

L = LeafletFinder(u, 'name P*')

leaflet0 = L.groups(0)
leaflet1 = L.groups(1)

##Rgyr = []
Rgyr_POPG =[]
Rgyr_POPE =[]

Rgyr0 = []
Rgyr1 = []

Rgyr0_Heads = []
Rgyr0_Tails = []

Rgyr1_Heads = []
Rgyr1_Tails = []

##lipids = u.select_atoms("resid 1:42224") # bilayer

lipids_POPG = u.select_atoms("segid POPG")
lipids_POPE = u.select_atoms("segid POPE")

lipids0 = leaflet0.residues.atoms  #upper leaflet
lipids1 = leaflet1.residues.atoms #lower leaflet

lipids0_Heads = leaflet0.select_atoms("name N or name HN1 or name HN2 or name HN3 or name C12 or name H12A or name H12B or name C11 or name H11A or name H11B or name P or name O13 or name O14 or name O11 or "
                                      "name O12 or name C1 or name HA or name HB or name C2 or name HS or name O21 or name C21 or name O22 or name C22 or "
                                      "name H2R or name H2S or name C3 or name HX or name HY or name O31 or name C31 or name O32") # upper leaflet headgroup

lipids0_Tails = leaflet0.select_atoms("not name N or name HN1 or name HN2 or name HN3 or name C12 or name H12A or name H12B or name C11 or name H11A or name H11B or name P or name O13 or name O14 or name O11 or "
                                      "name O12 or name C1 or name HA or name HB or name C2 or name HS or name O21 or name C21 or name O22 or name C22 or "
                                      "name H2R or name H2S or name C3 or name HX or name HY or name O31 or name C31 or name O32") # upper leaflet tails

lipids1_Heads = leaflet1.select_atoms("name N or name HN1 or name HN2 or name HN3 or name C12 or name H12A or name H12B or name C11 or name H11A or name H11B or name P or name O13 or name O14 or name O11 or "
                                      "name O12 or name C1 or name HA or name HB or name C2 or name HS or name O21 or name C21 or name O22 or name C22 or "
                                      "name H2R or name H2S or name C3 or name HX or name HY or name O31 or name C31 or name O32") # lower leaflet headgroup

lipids1_Tails = leaflet1.select_atoms("not name N or name HN1 or name HN2 or name HN3 or name C12 or name H12A or name H12B or name C11 or name H11A or name H11B or name P or name O13 or name O14 or name O11 or "
                                      "name O12 or name C1 or name HA or name HB or name C2 or name HS or name O21 or name C21 or name O22 or name C22 or "
                                      "name H2R or name H2S or name C3 or name HX or name HY or name O31 or name C31 or name O32") # lower leaflet tails

for ts in u.trajectory:

    ##Rgyr.append((u.trajectory.time, lipids.radius_of_gyration()))

    Rgyr0.append((u.trajectory.time, lipids0.radius_of_gyration()))
    Rgyr1.append((u.trajectory.time, lipids1.radius_of_gyration()))

    # Radius of Gyration for upper(0) and lower(1) leaflets

    Rgyr0_Heads.append((u.trajectory.time, lipids0_Heads.radius_of_gyration()))
    Rgyr0_Tails.append((u.trajectory.time, lipids0_Tails.radius_of_gyration()))

    # Radius of Gyration of headgroups and tails for upper leaflet

    Rgyr1_Heads.append((u.trajectory.time, lipids1_Heads.radius_of_gyration()))
    Rgyr1_Tails.append((u.trajectory.time, lipids1_Tails.radius_of_gyration()))

    # Radius of Gyration of headgroups and tails for lower leaflet

##Rgyr = np.array(Rgyr)

Rgyr0 = np.array(Rgyr0)
Rgyr1 = np.array(Rgyr1)

Rgyr0_Heads = np.array(Rgyr0_Heads)
Rgyr0_Tails = np.array(Rgyr0_Tails)

Rgyr1_Heads = np.array(Rgyr1_Heads)
Rgyr1_Tails = np.array(Rgyr1_Tails)

# 2D array of Radius of Gyration values

##a, b,  rho = arma_estimate(Rgyr.flatten('F'), 20, 20, 40)

a0,b0, rho0 = arma_estimate(Rgyr0.flatten('F'), 20, 20, 40)
a1,b1, rho1 = arma_estimate(Rgyr1.flatten('F'), 20, 20, 40)

a0_Heads,b0_Heads, rho0_Heads = arma_estimate(Rgyr0_Heads.flatten('F'), 20, 20, 40)
a0_Tails,b0_Tails, rho0_Tails = arma_estimate(Rgyr0_Tails.flatten('F'), 20, 20, 40)

a1_Heads,b1_Heads, rho1_Heads = arma_estimate(Rgyr1_Heads.flatten('F'), 20, 20, 40)
a1_Tails,b1_Tails, rho1_Tails = arma_estimate(Rgyr1_Tails.flatten('F'), 20, 20, 40)

# ARMA(20, 20) Model; Autoregressive and Moving Average parameters. Note: flatten('F') takes the 2D array "Gyration Tensor" and vectorizes it.

##psd  = arma2psd(A = a,B=b,  rho=rho,  sides ='centerdc', NFFT=4096)

psd0 = arma2psd(A=a0, B=b0, rho=rho0, sides ='centerdc', NFFT=4096)
psd1 = arma2psd(A=a1, B=b1, rho=rho1, sides ='centerdc', NFFT=4096)

psd0_Heads = arma2psd(A=a0_Heads, B=b0_Tails, rho=rho0_Heads, sides ='centerdc', NFFT=4096)
psd0_Tails = arma2psd(A=a0_Tails, B=b0_Tails, rho=rho0_Tails, sides ='centerdc', NFFT=4096)

psd1_Heads = arma2psd(A=a1_Heads, B=b1_Heads, rho=rho1_Heads, sides ='centerdc', NFFT=4096)
psd1_Tails = arma2psd(A=a1_Tails, B=b1_Tails, rho=rho1_Tails, sides ='centerdc', NFFT=4096)

# Power Spectral Density (PSD)


UL, = plt.plot(10*np.log10(psd0/max(psd0)),'b-', label= "Upper Leaflet")
LL, = plt.plot(10*np.log10(psd1/max(psd1)), 'r-', label= "Lower Leaflet")
##plt.plot(10*np.log(psd/max(psd)), 'g-')

first_legend = plt.legend(handles=[UL,LL], loc=1)
ax = plt.gca().add_artist(first_legend)


plt.xlim([0,4000])
plt.xlabel("$f$")
plt.ylabel(r"$P_{ARMA}$ ($f$)")
plt.title("Power Spectrum of $R_G$ ($t$) for a POPG bilayer")

plt.show()
plt.savefig('ARMA_Rgyr_leafletsANDbilayer.png', facecolor='lightgrey', edgecolor='w')
plt.close()