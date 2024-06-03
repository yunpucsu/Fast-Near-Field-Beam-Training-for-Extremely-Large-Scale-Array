This simulation code package is mainly used to reproduce the results of the following paper [1]:

[1] Y. Zhang, X. Wu, and C. You, “Fast near-field beam training for extremely large-scale array,” IEEE Wireless Commun. Lett., vol. 11, no. 12, pp. 2625–2629, Dec. 2022.

*********************************************************************************************************************************
If you use this simulation code package in any way, please cite the original paper [1] above. 
 
The author in charge of this simulation code package is Yunpu Zhang (email: yunpuzhangcsu@gmail.com).

Please note that the MATLAB R2023a is used for this simulation code package,  
and there may be some incompatibility problems among different MATLAB versions. 

*********************************************************************************************************************************
Abstract of the paper: 

In this letter, we study efficient near-field beam training design for the extremely large-scale array (XL-array)
communication systems. Compared with the conventional far-field beam training method that searches for the best beam
direction only, the near-field beam training is more challenging since it requires a beam search over both the angular and
distance domains due to the spherical wavefront propagation model. To reduce the near-field beam-training overhead based
on the two-dimensional exhaustive search, we propose in this letter a new two-phase beam training method that decomposes the
two-dimensional search into two sequential phases. Specifically, in the first phase, the candidate angles of the user are determined by
a new method based on the conventional far-field codebook and angle-domain beam sweeping. Then, a customized polar-domain
codebook is employed in the second phase to find the best effective distance of the user given the shortlisted candidate angles.
Numerical results show that our proposed two-phase beam training method significantly reduces the training overhead of the
exhaustive search and yet achieves comparable beamforming performance for data transmission. 

*********************************************************************************************************************************
How to use this simulation code package?

The simulation results in Figs. 3(b) and 3(d) can be obtained by running the files 'Main_Rate_vs_SNR.m' and 'Main_Rate_vs_Distance.m', respectively.

*********************************************************************************************************************************
Enjoy the reproducible research!