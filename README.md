# paper_GRL_Verger-Miralles_etal_2025
Codes to accompany "SWOT enhances small-scale eddy detection in the Mediterranean Sea"

figure1.ipynb to figure5.ipynb
→ Computations and generation of plots to reproduce the manuscript figures.

adcp_compute_rmsd_all.ipynb
→ Computes total RMSD between 6 ADCP transects' cross-section velocities and SWOT/DUACS-derived cross-section velocities. Includes the bootstrap estimation of statistical error associated with the RMSD computation.

drifters_compute_rmsd_all_duacs.ipynb
→ Computes RMSD between SVPB drifter-derived velocities module and DUACS-derived velocities. Includes the bootstrap estimation of statistical error associated with the RMSD computation.

drifters_compute_rmsd_all_swot.ipynb
→ Computes RMSD between SVPB drifter-derived velocities module and SWOT-derived velocities. Includes the bootstrap estimation of statistical error associated with the RMSD computation.

functions.py
→ Utility functions for data processing and RMSD computation used across notebooks.

