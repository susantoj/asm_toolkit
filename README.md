Asynchronous Motor Parameter Estimation Toolkit
===============================================

A MATLAB-based toolkit for estimating induction motor equivalent circuit parameters from manufacturer data.
All algorithms solve for the 8 parameter double-cage model.

Files in the repository:
- dnr_solver.m:         Damped Newton-Raphson solver (with linear restrictions kr and kx)
- dnr_solver2.m:        Damped Newton-Raphson solver (with Rs and Xr2 known)
- hybrid_dnr.m:         Hybrid Damped Newton-Raphson & Genetic algorithm solver
- hybrid_nr.m:          Hybrid Newton-Raphson & Genetic algorithm solver
- hybrid_nr.m:          Hybrid Levenberg-Marquardt & Genetic algorithm solver 
- ga_solver.m:          Genetic algorithm solver
- lm_solver.m:          Levenberg-Marquardt solver (with error term adjustment and linear restrictions kr and kx)
- lm_solver1a.m:        Levenberg-Marquardt solver (with error term adjustment and Rs and Xr2 known)
- lm_solver2.m:         Levenberg-Marquardt solver (with gain ratio adjustment and linear restrictions kr and kx)
- nr_solver.m:          Newton-Raphson solver (with linear restrictions kr and kx)
- nr_solver2.m:         Newton-Raphson solver (with Rs and Xr2 known)
- calc_pqt.m:           Calculates mechanical power, reactive power, current and torque from equivalent circuit parameters
- get_torque:           Calculates stator current and torque from equivalent circuit parameters
- IEC_test_set.csv      Test set of 4,002 IEC motors
- NEMA_test_set.csv     Test set of 2,378 NEMA motors

Author
======

Julius Susanto
+ http://github.com/susantoj
