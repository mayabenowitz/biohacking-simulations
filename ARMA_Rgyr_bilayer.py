import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
from spectrum import arma_estimate, arma, arma2psd
from MDAnalysis.analysis.leaflet import LeafletFinder

u = MDAnalysis.Universe('/home/jarod/Documents/Simulation/pepg_bilayer/step7_1/bilayer.gro',
                        '/home/jarod/Documents/Simulation/pepg_bilayer/step7_1/bilayer.xtc' )#Simulation time-series of a 2:1 POPE/POPG lipid bilayer.

L = LeafletFinder(u, 'name P*')#Parses .gro files for Phosphorus groups P*, i.e., the lipid heads.

leaflet0 = L.groups(0)#Upper leaflet.
leaflet1 = L.groups(1)#Lower leaflet.

#Initialzes POPE and POPG (types of lipids) arrays.

Rgyr_POPE = []
Rgyr_POPG = []

#Initializes upper and lower leaflet arrays.

Rgyr_leaflet0 = []
Rgyr_leaflet1 = []

#Initializes POPE and POPG upper leaflet arrays.

Rgyr_leaflet0_POPE = []
Rgyr_leaflet0_POPG = []

#Initializes POPE and POPG lower leaflet arrays.

Rgyr_leaflet1_POPE = []
Rgyr_leaflet1_POPG = []

#Instance 1: grabs atoms of the POPE and POPG lipids from the entire membrane.

lipids_POPE = u.select_atoms("resname POPE")
lipids_POPG = u.select_atoms("resname POPG")

#Instance 2: grabs all atoms in upper and lower leaflets.

lipids0 = leaflet0.residues.atoms
lipids1 = leaflet1.residues.atoms

#Instance 3: grabs atoms of the POPE and POPG lipids in the upper leaflet.

lipids0_POPE = leaflet0.select_atoms("resname POPE")
lipids0_POPG = leaflet0.select_atoms("resname POPG")

#Instance 4: grabs atoms of the POPE and POPG lipids in the lower leaflet.

lipids1_POPE = leaflet1.select_atoms("resname POPE")
lipids1_POPG = leaflet1.select_atoms("resname POPG")

#Calculates the radius of gyration (i.e., the shape) of the above instances for every
#10 nanosecond (ns) interval of the simulation and stores the numerical values (in angstroms) in the above arrays.

for ts in u.trajectory:

    Rgyr_POPE.append((u.trajectory.time, lipids_POPE.radius_of_gyration()))
    Rgyr_POPG.append((u.trajectory.time, lipids_POPG.radius_of_gyration()))

    Rgyr_leaflet0.append((u.trajectory.time, lipids0.radius_of_gyration()))
    Rgyr_leaflet1.append((u.trajectory.time, lipids1.radius_of_gyration()))

    Rgyr_leaflet0_POPE.append((u.trajectory.time, lipids0_POPE.radius_of_gyration()))
    Rgyr_leaflet0_POPG.append((u.trajectory.time, lipids0_POPG.radius_of_gyration()))

    Rgyr_leaflet1_POPE.append((u.trajectory.time, lipids1_POPE.radius_of_gyration()))
    Rgyr_leaflet1_POPG.append((u.trajectory.time, lipids1_POPG.radius_of_gyration()))

Rgyr_POPE = np.array(Rgyr_POPE)
Rgyr_POPG = np.array(Rgyr_POPG)

Rgyr_leaflet0 = np.array(Rgyr_leaflet0)
Rgyr_leaflet1 = np.array(Rgyr_leaflet1)

Rgyr_leaflet0_POPE = np.array(Rgyr_leaflet0_POPE)
Rgyr_leaflet0_POPG = np.array(Rgyr_leaflet0_POPG)

Rgyr_leaflet1_POPE = np.array(Rgyr_leaflet1_POPE)
Rgyr_leaflet1_POPG = np.array(Rgyr_leaflet1_POPG)

#ARMA model of above instances. (See documentation linked in the description for discussion of parameter choice.)

ar_POPE, ma_POPE, rho_POPE = arma_estimate(Rgyr_POPE.flatten('C'), 300, 300, 600)
ar_POPG, ma_POPG, rho_POPG = arma_estimate(Rgyr_POPG.flatten('C'), 300, 300, 600)

ar_leaflet0, ma_leaflet0, rho_leaflet0 = arma_estimate(Rgyr_leaflet0.flatten('C'), 300, 300, 600)
ar_leaflet1, ma_leaflet1, rho_leaflet1 = arma_estimate(Rgyr_leaflet1.flatten('C'), 300, 300, 600)

ar_leaflet0_POPE, ma_leaflet0_POPE, rho_leaflet0_POPE = arma_estimate(Rgyr_leaflet0_POPE.flatten('C'), 300, 300, 600)
ar_leaflet0_POPG, ma_leaflet0_POPG, rho_leaflet0_POPG = arma_estimate(Rgyr_leaflet0_POPG.flatten('C'), 300, 300, 600)

ar_leaflet1_POPE, ma_leaflet1_POPE, rho_leaflet1_POPE = arma_estimate(Rgyr_leaflet1_POPE.flatten('C'), 300, 300, 600)
ar_leaflet1_POPG, ma_leaflet1_POPG, rho_leaflet1_POPG = arma_estimate(Rgyr_leaflet1_POPG.flatten('C'), 300, 300, 600)

#Power spectrum analysis of above ARMA model. (See documentation linked in the description for discussion of parameter choice.)

psd_POPE = arma2psd(A=ar_POPE, B=ma_POPE, rho=rho_POPE, sides ='centerdc', NFFT=4096)
psd_POPG = arma2psd(A=ar_POPG, B=ma_POPG, rho=rho_POPG, sides='centerdc', NFFT=4096)

psd_leaflet0 = arma2psd(A=ar_leaflet0, B=ma_leaflet0, rho=rho_leaflet0, sides ='centerdc', NFFT=4096)
psd_leaflet1 = arma2psd(A=ar_leaflet1, B=ma_leaflet1, rho=rho_leaflet1, sides='centerdc', NFFT=4096)

psd_leaflet0_POPE = arma2psd(A=ar_leaflet0_POPE, B=ma_leaflet0_POPE, rho=rho_leaflet0_POPE, sides ='centerdc', NFFT=4096)
psd_leaflet0_POPG = arma2psd(A=ar_leaflet0_POPG, B=ma_leaflet0_POPG, rho=rho_leaflet0_POPG, sides ='centerdc', NFFT=4096)

psd_leaflet1_POPE = arma2psd(A=ar_leaflet1_POPE, B=ma_leaflet1_POPE, rho=rho_leaflet1_POPE, sides ='centerdc', NFFT=4096)
psd_leaflet1_POPG = arma2psd(A=ar_leaflet1_POPG, B=ma_leaflet1_POPG, rho=rho_leaflet1_POPG, sides ='centerdc', NFFT=4096)

#Plots results of above analyses.

plt.subplot(2,1,1)
POPE, = plt.loglog(psd_POPE, 'r-', label="POPE")
POPG, = plt.loglog(psd_POPG, 'b-', label="POPG")
plt.xlim([0, 4000])
plt.xlabel("frequency (Hz)")
plt.ylabel(r"$P_{ARMA}$ ($f$)")
plt.title("Power Spectrum")

first_legend = plt.legend(handles=[POPE,POPG], loc=1)
ax1 = plt.gca().add_artist(first_legend)

plt.subplot(2,1,2)
plt.plot(Rgyr_POPE[:,0], Rgyr_POPE[:,1], 'r-',lw=2, label=r'"$R_G$')
plt.plot(Rgyr_POPG[:,0], Rgyr_POPG[:,1], 'b-', lw=2, label=r'"$R_G$')
plt.xlabel("time (ps)")
plt.ylabel(r"$R_G$ ($\AA$)")
plt.title("Radius of Gyration Time Series of a POPE/POPG 2:1 Bilayer")

plt.tight_layout()
plt.show()
plt.savefig('ARMA_Rgyr_POPE-POPG.png', facecolor='lightgrey', edgecolor='w')
plt.close()

plt.subplot(2,1,1)
Upper_Leaflet, = plt.loglog(psd_leaflet0, 'r-', label="Upper Leaflet")
Lower_Leaflet, = plt.loglog(psd_leaflet1, 'b-', label="Lower Leaflet")
plt.xlim([0, 4000])
plt.xlabel("frequency (Hz)")
plt.ylabel(r"$P_{ARMA}$ ($f$)")
plt.title("Power Spectrum")

second_legend = plt.legend(handles = [Upper_Leaflet, Lower_Leaflet], loc=1)
ax2 = plt.gca().add_artist(second_legend)

plt.subplot(2,1,2)
plt.plot(Rgyr_leaflet0[:,0], Rgyr_leaflet0[:,1], 'r-',lw=2, label=r'"$R_G$')
plt.plot(Rgyr_leaflet1[:,0], Rgyr_leaflet1[:,1], 'b-', lw=2, label=r'"$R_G$')
plt.xlabel("time (ps)")
plt.ylabel(r"$R_G$ ($\AA$)")
plt.title("Radius of Gyration Time Series of a POPE/POPG 2:1 Bilayer")

plt.tight_layout()
plt.show()
plt.savefig('ARMA_Rgyr_Upper-Lower_Leaflets.png', facecolor='lightgrey', edgecolor='w')
plt.close()

plt.subplot(2,1,1)
Upper_Leaflet_POPE, = plt.loglog(psd_leaflet0_POPE, 'r-', label="Upper Leaflet (POPE)")
Lower_Leaflet_POPE, = plt.loglog(psd_leaflet1_POPE, 'b-', label="Lower Leaflet (POPE)")
plt.xlim([0, 4000])
plt.xlabel("frequency (Hz)")
plt.ylabel(r"$P_{ARMA}$ ($f$)")
plt.title("Power Spectrum")

third_legend = plt.legend(handles = [Upper_Leaflet_POPE, Lower_Leaflet_POPE], loc=1)
ax3 = plt.gca().add_artist(third_legend)

plt.subplot(2,1,2)
plt.plot(Rgyr_leaflet0_POPE[:,0], Rgyr_leaflet0_POPE[:,1], 'r-',lw=2, label=r'"$R_G$')
plt.plot(Rgyr_leaflet1_POPE[:,0], Rgyr_leaflet1_POPE[:,1], 'b-', lw=2, label=r'"$R_G$')
plt.xlabel("time (ps)")
plt.ylabel(r"$R_G$ ($\AA$)")
plt.title("Radius of Gyration Time Series of POPE")

plt.tight_layout()
plt.show()
plt.savefig('ARMA_Rgyr_Upper-Lower_Leaflets-POPE.png', facecolor='lightgrey', edgecolor='w')
plt.close()

plt.subplot(2,1,1)
Upper_Leaflet_POPG, = plt.loglog(psd_leaflet0_POPG, 'r-', label="Upper Leaflet (POPG)")
Lower_Leaflet_POPG, = plt.loglog(psd_leaflet1_POPG, 'b-', label="Lower Leaflet (POPG)")
plt.xlim([0, 4000])
plt.xlabel("frequency (Hz)")
plt.ylabel(r"$P_{ARMA}$ ($f$)")
plt.title("Power Spectrum")

fourth_legend = plt.legend(handles = [Upper_Leaflet_POPG, Lower_Leaflet_POPG], loc=1)
ax4 = plt.gca().add_artist(fourth_legend)

plt.subplot(2,1,2)
plt.plot(Rgyr_leaflet0_POPG[:,0], Rgyr_leaflet0_POPG[:,1], 'r-',lw=2, label=r'"$R_G$')
plt.plot(Rgyr_leaflet1_POPG[:,0], Rgyr_leaflet1_POPG[:,1], 'b-', lw=2, label=r'"$R_G$')
plt.xlabel("time (ps)")
plt.ylabel(r"$R_G$ ($\AA$)")
plt.title("Radius of Gyration Time Series of POPG")

plt.tight_layout()
plt.show()
plt.savefig('ARMA_Rgyr_Upper-Lower_Leaflets-POPG.png', facecolor='lightgrey', edgecolor='w')
plt.close()
