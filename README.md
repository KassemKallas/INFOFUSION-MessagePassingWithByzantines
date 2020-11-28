# INFOFUSION-MessagePassingWithByzantines
 Simulation files for the paper "A message passing approach for decision fusion in adversarial multi-sensor networks" published at Information Fusion - ElSevier

 Authors: Andrea Abrardo, Mauro Barni, Kassem Kallas, Benedetta Tondi

 Abstract: We consider a simple, yet widely studied, set-up in which a Fusion Center (FC) is asked to make a binary decision about a sequence of system states by relying on the possibly corrupted decisions provided by byzantine nodes, i.e. nodes which deliberately alter the result of the local decision to induce an error at the fusion center. When independent states are considered, the optimum fusion rule over a batch of observations has already been derived, however its complexity prevents its use in conjunction with large observation windows.

In this paper, we propose a near-optimal algorithm based on message passing that greatly reduces the computational burden of the optimum fusion rule. In addition, the proposed algorithm retains very good performance also in the case of dependent system states. By first focusing on the case of small observation windows, we use numerical simulations to show that the proposed scheme introduces a negligible increase of the decision error probability compared to the optimum fusion rule. We then analyse the performance of the new scheme when the FC makes its decision by relying on long observation windows. We do so by considering both the case of independent and Markovian system states and show that the obtained performance are superior to those obtained with prior suboptimal schemes. As an additional result, we confirm the previous finding that, in some cases, it is preferable for the byzantine nodes to minimise the mutual information between the sequence system states and the reports submitted to the FC, rather than always flipping the local decision.

The repository contains two folders: one for the main simulations and the other to the run the experiments of the Byzantines dual behavior as reported in the paper.

For the main, set the parameters in the file Main_sim_turbo_def3.m
The rest of the files in both directories are to compute the messages exchanged over the factor graph and to simulate the state of art results.


The paper can be found at: https://www.sciencedirect.com/science/article/abs/pii/S156625351630135X

For any inquiries please contact Dr. Kassem Kallas at k_kallas@hotmail.com
