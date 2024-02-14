This repository holds code used to generate the figures and tables in our paper "Variational design of sensory feedback for powerstroke-recovery systems". Figures 2-8 are for the HCO model. Figures 10-16 are for the Markin model. The rest figures are not generated via simulation.

Fig. 2 shows a typical solution of the HCO-motor system. To generate figure, run Fig_2.m

Fig. 3 shows the trajectoreis on the same time scale and the iSRC solution. To generate figure, run Fig_3.m

Fig. 4 shows all performance-sensitivity patterns in the HCO model. The data of this figure is generated using sensitivity_Yu.m as the simulation code. For Fig. 4, to produce the plot from precomputed data (excitation-contralateral.mat and inhibition-contralateral.mat), run Fig_4.m

Fig. 5 compares trajectories for the constant inhibitory feedback system and constant excitatory feedback system. To generate figure, run Fig_5.m

Fig. 6 compares trajectories of two systems controlled by the inhibitory-decreasing mechanism with L_0=11 and L_{slope}=0.5 versus L_{slope}=2. To generate figure, run Fig_6.m

Fig. 7 shows the effects of the increased load on the period and on the progress of the trajectories for the excitatory-decreasing systems and inhibitory-decreasing systems with L_0 fixed at L_0=9 and L_{slope} varied over (0.6, 40). The data of this figure is generated using sensitivity_Yu.m as the simulation code. For Fig. 7, to produce the plot from precomputed data (excitation-contralateral.mat and inhibition-contralateral.mat), run Fig_7.m

Fig. 8 shows the detail from Fig. 4A, expanded to show the region of approximate linearity. To generate figure, run Fig_8.m

Fig. 10 shows a typical solution of the Markin-motor system. To generate figure, run Fig_10.m

Fig. 11 compares trajectories of the unperturbed and perturbed systems. To generate figure, run Fig_11.m

Fig. 12 shows all performance-sensitivity patterns in the Markin model for different feedback strengths. The data of this figure is generated using sensitivity_Markin.m as the simulation code. For Fig. 12, to produce the plot from precomputed data (FB_strength_corrected.mat), run Fig_12.m

Fig. 13 shows the escape of In triggering CPG transitions. To generate figure, run Fig_13.m

Fig. 14 shows the effects of reduced sensory feedback gain s_{Ia-F} to In cells near the collapse point. To generate figure, run Fig_14.m

Fig. 15 compares the trajectories of two systems with s_{Iaf}=1 versus s_{Iaf}=1.1. To generate figure, run Fig_15.m

Fig. 16 shows the effects of the Ib-E feedback strength on the timing and on the shape of the solution trajectories for the Markin model. The data of this figure is generated using sensitivity_Markin.m as the simulation code. For Fig. 16, to produce the plot from precomputed data (s_Ibe.mat), run Fig_16.m
