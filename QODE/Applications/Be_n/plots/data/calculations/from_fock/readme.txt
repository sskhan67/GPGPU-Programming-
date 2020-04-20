	from_fock.tar.gz

This is a relatively complete export of the directory in which I was running under yhliu@fock.  The only things missing that would be necessary to reproduce
the results would be the psi4 code itself, the python3 scipy stack, and the Hamiltonian files, which are the same as what are being used here (exported as
a tarball ... only the directory /morescratch/adutoi/programs/Qode/Applications/Be_n/dimer_H/6-31g/run/molecule_CI/H/16-115-550/load=H:16-115-550:thresh=1e-6:4.5:u.pickle/
was copied, and installed on fock under .../Qode/Applications/Be_n/dimer_H/6-31g/run/molecule_CI/H/16-115-550/load=H:16-115-550:thresh=1e-6:4.5:u.pickle/).

The following things were extracted from this tarball (directory names changed):

	psi4/

The psi4 timings are from runs that were forced to use a single core.
It was observed that, at CCSD convergence, energies were stable out to the 13th decimal place.
The timings given are for the individual CCSD and (T) modules (ie the total CCSD(T) time would be the sum of these plus the HF and integral transformation times).

	calculations/

The log files in curves/ should be exact matches to the calculations here (save timing and numerical noise), and have been brought over just as a veracity check.
Also the calculations in timings/bf_thresh_changed/ should match the ones done here.

The log files in timings/ are for calculations in which the X-CCSD convergence was lowered to 1e-12 (relative change of T amplitude vector norm), because this
was observed to converge around the same time that energies became stable in the 13th decimal place, matching the quality of the psi4 answer.
To mitigate any potential complaints about "cheating," I also lowered the threshold on Hamiltonian matrix elements to 1e-16, which interestingly has
 only a very slight effect on timings (indicating that those elements that are thresholded out, which does make a big efficiency difference, are indeed very small).
