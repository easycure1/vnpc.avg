# vnpc.avg

Obtain samples of the posterior of the averge multivariate likelihood in conjunction with an Hpd AGamma process prior on the spectral density matrix.

This work is a follow-up of [Yixuan Liu et al '23] and [Jianan Liu et al '24].



[Yixuan Liu et al '23]: http://dx.doi.org/10.1016/j.csda.2024.108010
[Jianan Liu et al '24]: https://arxiv.org/abs/2409.13224

```
@article{YixuanLiu_2024,
   title={A nonparametrically corrected likelihood for Bayesian spectral analysis of multivariate time series},
   volume={199},
   ISSN={0167-9473},
   url={http://dx.doi.org/10.1016/j.csda.2024.108010},
   DOI={10.1016/j.csda.2024.108010},
   journal={Computational Statistics &amp; Data Analysis},
   publisher={Elsevier BV},
   author={Liu, Yixuan and Kirch, Claudia and Lee, Jeong Eun and Meyer, Renate},
   year={2024},
   month=nov, pages={108010} }
   
@article{JiananLiu_2024,
      title={Variational inference for correlated gravitational wave detector network noise}, 
      author={Jianan Liu and Avi Vajpeyi and Renate Meyer and Kamiel Janssens and Jeung Eun Lee and Patricio Maturana-Russel and Nelson Christensen and Yixuan Liu},
      year={2024},
      eprint={2409.13224},
      archivePrefix={arXiv},
      primaryClass={gr-qc},
      url={https://arxiv.org/abs/2409.13224}, 
}

```


## Developer notes


First install the dependencies

```
install.packages("devtools")
library(devtools)
devtools::install_deps()
```


Install the package
```
devtools::install(quick = TRUE)
```