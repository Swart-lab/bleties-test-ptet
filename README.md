Benchmarking BleTIES with simulated reads from Paramecium
=========================================================

Source data and setup
---------------------

This repository supplements our
[preprint](https://www.biorxiv.org/content/10.1101/2021.05.18.444610v1)
describing the [BleTIES](https://github.com/Swart-lab/bleties) software.

MAC and MAC+IES assemblies (both v1) from Paramecium tetraurelia strain 51,
obtained from [ParameciumDB](https://paramecium.i2bc.paris-saclay.fr)

Download and install the following software in the `opt/` folder:

 * ART read simulator https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
 * PBSIM read simulator (v2): `git clone git@github.com:yukiteruono/pbsim2.git`
 * minimap2 mapper: `git clone https://github.com/lh3/minimap2`
 * seqtk: `git clone git@github.com:lh3/seqtk.git`
 * `faFilter` and `faSomeRecords` from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

Run read simulation and mapping, documented in script `do_bleties.sh`


Analysis
--------

Evaluations of BleTIES results and comparison with original sequences are in
the Jupyter notebooks:
 * `Get_input_IES_coordinates.ipynb`
 * `Evaluate_BleTIES_MILRAA_with_simulated_ONT_R9.5_reads.ipynb`
 * `Evaluate_BleTIES_MILRAA_with_simulated_PB_reads.ipynb`
