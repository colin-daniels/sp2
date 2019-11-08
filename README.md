# sp<sup>2</sup>
![](https://github.com/colin-daniels/sp2/workflows/C%2FC%2B%2B%20CI/badge.svg)

Computational physics research application dealing mainly with carbon-based structures 
(e.g. graphene, carbon nanotubes) developed at the ICMP research group by Colin Daniels and Michael Lamparski. A web interface to some of the Raman functionality is available [here](https://icmp.phys.rpi.edu/sp2).

> **Disclaimer**: sp<sup>2</sup> is very much research code, and has hardly followed best practices in terms of development. As such, some features may not work and it is unfortunately very possible that the exact code used for a particular publication may not even be present in the git commit history. The code itself is reproduced here for reference, and in particular the Raman functionality has been superceded by [rsp2](https://github.com/ExpHP/rsp2). 

## Installation
There are only a few dependencies, the primary ones being `MPI`, `boost`, `python 3`, and `phonopy`.
Look at the [GitHub workflow file](.github/workflows/ccpp.yml) if totally lost, since it does build fine in that environment.

## Related Publications
Below is a list of publications in which this software was utilized at least in part.

- J. Overbeck, G. B. Barin, C. Daniels, M. L. Perrin, O. Braun, Q. Sun, R. Darawish, M. De Luca, X.-Y. Wang, T. Dumslaff, A. Narita, K. Müllen, P. Ruffieux, V. Meunier, R. Fasel, and M. Calame, _A Universal Length-Dependent Vibrational Mode in Graphene Nanoribbons_, ACS Nano [10.1021/acsnano.9b05817](https://doi.org/10.1021/acsnano.9b05817) (2019).
- O. Gröning, S. Wang, X. Yao, C. A. Pignedoli, G. B. Barin, C. Daniels, A. Cupo, V. Meunier, X. Feng, A. Narita, K. Müllen, P. Ruffieux, and R. Fasel, _Engineering of Robust Topological Quantum Phases in Graphene Nanoribbons_, Nature **560**, 209 [10.1038/s41586-018-0375-9](https://doi.org/10.1038/s41586-018-0375-9) (2018).
- J. R. Owens, C. Daniels, A. Nicolaï, H. Terrones, and V. Meunier, _Structural, Energetic, and Electronic Properties of Gyroidal Graphene Nanostructures_, Carbon **96**, 998 [10.1016/j.carbon.2015.10.042](https://doi.org/10.1016/j.carbon.2015.10.042) (2016).
- Z. J. Qi, C. Daniels, S. J. Hong, Y. W. Park, V. Meunier, M. Drndić, and A. T. C. Johnson, _Electronic Transport of Recrystallized Freestanding Graphene Nanoribbons_, ACS Nano **9**, 3510 [10.1021/nn507452g](https://doi.org/10.1021/nn507452g) (2015).
- A. Nicolaï, J. Monti, C. Daniels, and V. Meunier, _Electrolyte Diffusion in Gyroidal Nanoporous Carbon_, J. Phys. Chem. C **119**, 2896 [10.1021/jp511919d](https://doi.org/10.1021/jp511919d) (2015).
- C. Daniels, A. Horning, A. Phillips, D. V. P. Massote, L. Liang, Z. Bullard, B. G. Sumpter, and V. Meunier, _Elastic, Plastic, and Fracture Mechanisms in Graphene Materials_, J. Phys.: Condens. Matter **27**, 373002 [10.1088/0953-8984/27/37/373002](https://doi.org/10.1088/0953-8984/27/37/373002) (2015).
- C. Daniels, Z. Bullard, E. C. Girão, and V. Meunier, _Emergent Magnetism in Irradiated Graphene Nanostructures_, Carbon **78**, 196 [10.1016/j.carbon.2014.06.072](https://doi.org/10.1016/j.carbon.2014.06.072) (2014).

## References
sp<sup>2</sup> would not be possible without building upon the following works:

- A. N. Kolmogorov and V. H. Crespi, _Registry-Dependent Interlayer Potential for Graphitic Systems_, Phys. Rev. B **71**, 235415 [10.1103/PhysRevB.71.235415](https://doi.org/10.1103/PhysRevB.71.235415) (2005).
- A. Togo and I. Tanaka, _First Principles Phonon Calculations in Materials Science_, Scripta Materialia **108**, 1 [10.1016/j.scriptamat.2015.07.021](https://doi.org/10.1016/j.scriptamat.2015.07.021) (2015).
- D. W. Brenner, O. A. Shenderova, J. A. Harrison, S. J. Stuart, B. Ni, and S. B. Sinnott, _A Second-Generation Reactive Empirical Bond Order (REBO) Potential Energy Expression for Hydrocarbons_, J. Phys.: Condens. Matter **14**, 783 [10.1088/0953-8984/14/4/312](https://doi.org/10.1088/0953-8984/14/4/312) (2002).
- E. Bitzek, P. Koskinen, F. Gähler, M. Moseler, and P. Gumbsch, _Structural Relaxation Made Simple_, Phys. Rev. Lett. **97**, 170201 [10.1103/PhysRevLett.97.170201](https://doi.org/10.1103/PhysRevLett.97.170201) (2006).
- N. Andrei, _Conjugate gradient algorithms for molecular formation under pairwise potential minimization_, in Mathematical Modeling of Environmental and Life Sciences Problems [camo.ici.ro/neculai/potential.pdf](https://camo.ici.ro/neculai/potential.pdf) (Constanta, Romania, 2008).
- R. Saito, M. Furukawa, G. Dresselhaus, and M. S. Dresselhaus, _Raman Spectra of Graphene Ribbons_, J. Phys.: Condens. Matter **22**, 334203 [10.1088/0953-8984/22/33/334203](https://doi.org/10.1088/0953-8984/22/33/334203) (2010).
- S. Guha, J. Menéndez, J. B. Page, and G. B. Adams, _Empirical Bond Polarizability Model for Fullerenes_, Phys. Rev. B **53**, 13106 [10.1103/PhysRevB.53.13106](https://doi.org/10.1103/PhysRevB.53.13106) (1996).
- S. J. Stuart, A. B. Tutein, and J. A. Harrison, _A Reactive Potential for Hydrocarbons with Intermolecular Interactions_, The Journal of Chemical Physics **112**, 6472 [10.1063/1.481208](https://doi.org/10.1063/1.481208) (2000).
- S. Plimpton, _Fast Parallel Algorithms for Short-Range Molecular Dynamics_, Journal of Computational Physics **117**, 1 [10.1006/jcph.1995.1039](https://doi.org/10.1006/jcph.1995.1039) (1995).
- Z. Bullard and V. Meunier, _Dynamical Properties of Carbon Nanotube Welding into X Junctions_, Phys. Rev. B **88**, 035422 [10.1103/PhysRevB.88.035422](https://doi.org/10.1103/PhysRevB.88.035422) (2013).
- Z.-H. Zhan, J. Zhang, Y. Li, and H. S.-H. Chung, _Adaptive Particle Swarm Optimization_, IEEE Transactions on Systems, Man, and Cybernetics, Part B (Cybernetics) **39**, 1362 [10.1109/TSMCB.2009.2015956](https://doi.org/10.1109/TSMCB.2009.2015956) (2009).



