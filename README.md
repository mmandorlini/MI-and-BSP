# MI-and-BSP
Medical Imaging and Biomedical Signal Processing assignments, 2020 @Polimi

## Medical Imaging - Filtered Back Projection
Project focus: implement FBP algorithm (parallel ray) and evaluate its reconstruction quality with respect to noise and filter methods.

*Details:*

Implement a simple version of the FBP algorithm over the Shepp-Logan phantom without using radon and iradon matlab functions. \
Discuss how reconstruction performances vary according to different filtering methods, as well as the number of employed projections and image dimension. \
Add noise to the starting image, describing its physical interpretation, how it affects the reconstruction performances and what solutions can be implemented as post-processing steps.

## Biomedical Signal Processing - The potentialities of HRV as a tool for assessing fetal wellbeing
The focus of this project is on Intrauterine Growth Restriction (IUGR), a condition that affects the fetus growth and result in a baby small for gestational age (SGA). \
Efficient diagnosis and management of IUGR leads to perinatal mortality reduction; a digitalized signal processing approach is a helpful tool for clinicians in distinguishing between IUGR and Normal fetuses.

*CTG analysis pipeline:*
1. Pre-processing
2. Signal variance analysis - Eye inspection analysis
3. Signal variance analysis - Computerized system introduction (numerical digital context for variability analysis)
4. Spectral analysis - Introduction of methods in the frequency domain, giving emphasis on correlation between variability observed in HR signal and physiological mechanism
5. Nonlinear chaotic time series analysis - Introduction of new methodologies such as Approximate Entropy, Lempel Ziv complexity
6. Machine learning analysis - Principal components analysis and implementation of classification algorithm

Dataset: 10 CTG recordings. Traces were acquired at 2 Hz and labelled as healthy or IUGR depending on the clinical diagnosis at birth (5 healthy and 5 IUGR fetuses).
