# Optimal Control for Restarting Minimal Residual Krylov Solvers via Data-Driven Residual Dynamics

This repository contains the MATLAB implementation and numerical experiments presented in the paper: **"OPTIMAL CONTROL FOR RESTARTING MINIMAL RESIDUAL KRYLOV SOLVERS VIA DATA-DRIVEN RESIDUAL DYNAMICS"**. This framework introduces a data-driven approach to adaptively select the restart parameter $m$ in Krylov subspace solvers (specifically GMRES and LGMRES) by modeling the restart cycles as a discrete-time dynamical system. Restarted GMRES(m) is widely used for large, sparse, non-symmetric linear systems. However, choosing the optimal restart parameter $m$ is challenging; a static choice can lead to stagnation or inefficient convergence.

Our approach utilizes:
- **Dynamic Mode Decomposition with control (DMDc):** To identify a low-dimensional Internal Linear Surrogate Model (ILSM) that captures cycle-to-cycle residual dynamics.
- **Linear Quadratic Regulator (LQR):** An optimal control policy applied to the ILSM to compute an adaptive sequence of restart parameters $\{m_k\}$.

## MATLAB test scripts

The MATLAB scripts in this repository are designed to facilitate the testing and validation of the DMDc-LQR framework. They are organized into three main scripts:

- **`main_test_performance.m`**: Performs a detailed performance comparison for a specific matrix (e.g., `sherman5`). It executes four different solvers—standard GMRES, LGMRES, and their DMDc-LQR adaptive variants—to compare convergence rates and residual history. The output includes comparative plots and specific metrics on how the adaptive controller adjusts the restart parameter $m$ during the iteration process.

- **`main_test_collection.m`**: This script is a batch-processing tool designed for large-scale benchmarking across multiple matrices stored in a local folder. It automates the testing of the entire collection, records execution times and cycle counts, and exports the results into CSV files for further analysis. It also generates log-scale bar charts to visualize the relative speedup and efficiency of the DMDc-LQR strategy compared to fixed-restart baselines.

- **`main_test_fieldAlpha.m`**: This script focuses on sensitivity analysis regarding the LQR penalty parameter $\alpha$. By running simulations over a range of alpha values for a single problem, it helps determine the optimal balance between control effort and convergence speed. The results are visualized in a 3D plot that illustrates the relationship between the penalty value, the number of cycles, and the final residual.

## Reproducibility

We provide the following scripts to reproduce the numerical experiments and figures presented in the paper. The adaptive DMDc-LQR framework identifies an internal linear surrogate model (ILSM) to compute the optimal restart parameter sequence $\{m_k\}$.

| Figure / Section | Description | Reference Script |
| :--- | :--- | :--- |
| **Fig. 4.1 & 4.2** | Residual norm evolution vs. cycles for `sherman5` and `orsirr1` matrices. | `main_test_performance.m` |
| **Fig. 4.3** | Comparison of adaptive $m_k$ trajectories vs. fixed-restart baselines. | `main_test_performance.m` |
| **Fig. 4.6** | Sensitivity analysis of the LQR penalty parameter $\alpha$ (3D Surface plot). | `main_test_fieldAlpha.m` |
| **Table 2 & Fig. 4.5** | Large-scale benchmark results (Speedup/Efficiency) on SuiteSparse collection. | `main_test_collection.m` |
| **Appendix A** | Eigenvalue/Ritz value distribution and DMDc mode analysis for the ILSM. | `main_test_performance.m` |

> **Note:** The matrices used in these experiments are part of the SuiteSparse Matrix Collection. The scripts are configured to load them from a local `mat_collection_test_full/` folder or download them automatically where applicable.

## Citation

If you use this code or the DMDc-LQR framework in your research, please cite:

```
@article{benitez2025optimal,
  title={Optimal Control for Restarting Minimal Residual Krylov Solvers via Data-Driven Residual Dynamics},
  author={Carlos Benitez, Juan C. Cabral, Christian E. Schaerer},
  journal={Submitted to SIAM (Under Review)},
  year={2025}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

- Carlos Benítez - [carlos.benitez@fpuna.edu.py]
- Juan C. Cabral - [jcabral@fpuna.edu.py]
- Christian E. Schaerer - [cschaer@pol.una.py]

## References

[1] **B. Anderson and J. Moore**, *Optimal Control: Linear Quadratic Methods*, Prentice-Hall, Englewood Cliffs, NJ, 1990.

[2] **B. Bamieh**, *Dynamic Mode Decomposition with Control*, in Proceedings of the 51st IEEE Conference on Decision and Control, 2012.

[3] **S. T. Bathelt, J. P. Hummel, and R. D. Falgout**, *An adaptive restart strategy for GMRES*, Numer. Linear Algebra Appl., 28 (2021), e2346.

[4] **A. H. Baker, J. M. Dennis, and E. R. Jessup**, *On improving linear solver performance: A block variant of GMRES*, SIAM J. Sci. Comput., 27 (2006), pp. 1608–1626.

[5] **A. H. Baker, E. R. Jessup, and T. Manteuffel**, *A technique for accelerating the convergence of restarted GMRES*, SIAM J. Matrix Anal. Appl., 26 (2005), pp. 962–984.

[6] **A. H. Baker, E. R. Jessup, and T. V. Kolev**, *A simple strategy for varying the restart parameter in GMRES(m)*, J. Comput. Appl. Math., 230 (2009), pp. 751–761.

[7] **C. Benítez, J. C. Cabral, and C. E. Schaerer**, *Optimal Control for Restarting Minimal Residual Krylov Solvers via Data-Driven Residual Dynamics*, Submitted to SIAM Journal on Scientific Computing (SISC), 2025.

[8] **A. Bhaya and E. Kaszkurewicz**, *Control Perspectives on Numerical Algorithms and Matrix Problems*, SIAM, Philadelphia, PA, 2006.

[9] **A. Bhaya, P.-A. Bliman, and F. Pazos**, *Cooperative concurrent asynchronous computation of the solution of symmetric linear systems*, Numer. Algorithms, 75 (2017), pp. 587–617.

[10] **C. Brezinski and M. Redivo-Zaglia**, *Hybrid procedures for solving linear systems*, Numer. Math., 67 (1994), pp. 1–19.

[11] **J. C. Cabral, C. E. Schaerer, and A. Bhaya**, *Improving GMRES(m) using an adaptive switching controller*, Numer. Linear Algebra Appl., 27 (2020), e2305.

[12] **J. C. Cabral, P. Torres, and C. E. Schaerer**, *On control strategy for stagnated GMRES(m)*, J. Comp. Int. Sci., 14 (2017), pp. 109–118.

[13] **R. Cuevas Núñez, C. E. Schaerer, and A. Bhaya**, *A proportional-derivative control strategy for restarting the GMRES(m) algorithm*, J. Comput. Appl. Math., 337 (2018), pp. 209–224.

[14] **T. A. Davis and Y. Hu**, *The University of Florida Sparse Matrix Collection*, ACM Trans. Math. Software, 38 (2011), pp. 1–25.

[15] **J. Drkosova, A. Greenbaum, M. Rozlozník, and Z. Strakos**, *Numerical stability of GMRES*, BIT Numer. Math., 35 (1995), pp. 309–330.

[16] **K. Du, J. Duintjer Tebbens, and G. Meurant**, *Any admissible harmonic Ritz value set is possible for GMRES*, Electron. Trans. Numer. Anal., 47 (2017), pp. 37–56.

[17] **J. Ely**, *Restarted GMRES with adaptive restart parameter*, doctoral dissertation, 2008.

[18] **M. Embree**, *The Arnoldi eigenvalue algorithm with implicit restarts*, 2002.

[19] **V. Frayssé, L. Giraud, and S. Gratton**, *A set of adaptive restart strategies for GMRES*, Technical Report, CERFACS, 1998.

[20] **A. Gaul, M. H. Gutknecht, J. Liesen, and R. Nabben**, *A framework for deflated and augmented Krylov subspace methods*, SIAM J. Matrix Anal. Appl., 34 (2013), pp. 495–518.

[21] **L. Giraud, S. Gratton, and J. Langou**, *A reference implementation of some adaptive restart strategies for GMRES*, 2003.

[22] **G. H. Golub and C. F. Van Loan**, *Matrix Computations*, 4th ed., Johns Hopkins University Press, Baltimore, MD, 2013.

[23] **S. Goossens and D. Roose**, *Ritz and harmonic Ritz values and the convergence of FOM and GMRES*, Numer. Linear Algebra Appl., 6 (1999), pp. 281–293.

[24] **A. Greenbaum**, *Iterative Methods for Solving Linear Systems*, SIAM, Philadelphia, PA, 1997.

[25] **M. Habu and T. Nodera**, *GMRES (m) algorithm with changing the restart cycle adaptively*, Proc. Algorithmy 2000, (2000), pp. 254–263.

[26] **W. Joubert**, *On the convergence behavior of the restarted GMRES algorithm for solving nonsymmetric linear systems*, Numer. Linear Algebra Appl., 1 (1994), pp. 427–447.

[27] **S. A. Kharitonov**, *Adaptive selection of the restart parameter in the GMRES method*, 2005.

[28] **T. G. Kolda and J. R. Mayo**, *An Adaptive Shifted Power Method for Computing Generalized Tensor Eigenpairs*, SIAM J. Matrix Anal. Appl., 35 (2014), pp. 1563–1581.

[29] **S. P. Kolodziej, et al.**, *The SuiteSparse Matrix Collection website interface*, J. Open Source Softw., 4 (2019), 1244.

[30] **J. N. Kutz, S. L. Brunton, B. W. Brunton, and J. L. Proctor**, *Dynamic Mode Decomposition: Data-Driven Modeling of Complex Systems*, SIAM, 2016.

[31] **R. B. Morgan**, *A Restarted GMRES Method Augmented with Eigenvectors*, SIAM J. Matrix Anal. Appl., 16 (1995), pp. 1154–1171.

[32] **R. B. Morgan**, *GMRES with Deflated Restarting*, SIAM J. Sci. Comput., 24 (2002), pp. 20–37.

[33] **N. M. Nachtigal, S. C. Reddy, and L. N. Trefethen**, *How Fast are Nonsymmetric Matrix Iterations?*, SIAM J. Matrix Anal. Appl., 13 (1992), pp. 778–795.

[34] **J. L. Proctor, S. L. Brunton, and J. N. Kutz**, *Dynamic mode decomposition with control*, SIAM J. Appl. Dyn. Syst., 15 (2016), pp. 142–161.

[35] **Y. Saad**, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM, Philadelphia, PA, 2003.

[36] **Y. Saad and M. H. Schultz**, *GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems*, SIAM J. Sci. Stat. Comput., 7 (1986), pp. 856–869.

[37] **P. J. Schmid**, *Dynamic mode decomposition of numerical and experimental data*, J. Fluid Mech., 656 (2010), pp. 5–28.

[38] **V. Simoncini and D. B. Szyld**, *Recent computational developments in Krylov subspace methods for linear systems*, Numer. Linear Algebra Appl., 8 (2017), pp. 1–59.

[39] **G. L. G. Sleijpen and H. A. Van der Vorst**, *A Jacobi–Davidson Iteration Method for Linear Eigenvalue Problems*, SIAM Rev., 42 (2000), pp. 267–293.

[40] **J. C. Strikwerda and S. C. Stodder**, *Convergence results for GMRES (m)*, Tech. Report 1280, University of Wisconsin, 1995.

[41] **D. S. Watkins**, *The matrix eigenvalue problem: GR and Krylov subspace methods*, SIAM, Vol. 101, 2007.

[42] **L. Zhang and T. Nodera**, *A new adaptive restart for GMRES(m) method*, ANZIAM J., 46 (2005), pp. 409–425.
