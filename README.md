# SILC

Scale-discretised wavelet Internal Linear Combination

- deconvolve_beam.py: convolve input maps to same effective beam (also correct to thermodynamic temperature as necessary). Set nprocess = # input frequency channels. Input: map FITS files; beam transfer function Numpy array binary files. Output: map FITS files.
- s2let_analysis_hybrid.py: wavelet analysis of input maps. Set nprocess = # input frequency channels. Input: map FITS files. Output: Numpy array binary file per wavelet coefficient map.
- s2let_ilc_hybrid.py: internal linear combination (ILC) for each set of wavelet coefficient maps. Will loop over wavelet scale & orientation [from jmin_real: jmax_real+1; ndir_min: ndir_max+1]. Set nprocess = # directions per wavelet scale [defaults to 1]; nprocess2 = # input frequency channels; nprocess3 = # available cores per direction. Input: wavelet coefficient map Numpy array binary files. Output: wavelet coefficient map Numpy array binary file per wavelet scale & orientation.
- s2let_synthesis_hybrid.py: wavelet synthesis of ILC maps. Input: ILC wavelet coefficient map Numpy array binary files. Output: ILC map FITS file; ILC map angular power spectrum FITS file.


# Spin-SILC

Spin Scale-discretised wavelet Internal Linear Combination

- spin_silc_deconvolve.py
- spin_silc_spin_analysis.py
- spin_silc_ilc.py
- spin_silc_scalar_synthesis.py


Please cite 'Rogers et al. (2016a,b)' (https://ui.adsabs.harvard.edu/abs/2016MNRAS.460.3014R/abstract; https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.2310R/abstract).