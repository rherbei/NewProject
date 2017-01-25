This directory contains all data, code, and output pertinent to the paper
"Bayesian Inference of Selection in the Wright-Fisher Diffusion Model".
Details regarding the specific files are provided below.

Data:

Data/fly2Ln.txt
 - Text file containing the ending allele frequencies for the entire 
    2L chromosome for the neutral fly data

Data/fly2Ls.txt
 - Text file containing the ending allele frequencies for the entire 2L 
    2L chromosome for the hypoxic fly data

Data/fly-n.txt
 - Text file containing the ending allele frequencies for the subset of the 
    2L chromosome identified by Ronen et al (2013) for the neutral fly data

Data/fly-s.txt
 - Text file containing the ending allele frequencies for the subset of the 
    2L chromosome identified by Ronen et al (2013) for the hypoxic fly data

Data/nlocs2L.txt
 - Text file containing SNP locations for the entire 2L chromosome for the 
    neutral fly data

Data/slocs2L.txt
 - Text file containing SNP locations for the entire 2L chromosome for the 
    hypoxic fly data


Code:

WrightFisherSampling.cu
 - CUDA code for simulating the Wright-Fisher diffusion using the GPU
 - Produces "qSamples/final_samples_cuda.txt"

WFwithSelection-ApproxAll.py
 - Python code for simulating a general Wright-Fisher diffusion 
 - Uses a Normal approximation to compute m within the exact algorithm of
    Jenkins and Spano (2015)
 - Calls "NeutralWF-NormApproxAll.py"

NeutralWF-NormApproxAll.py
 - Python code for simulating a neutral Wright-Fisher diffusion
 - Uses a Normal approximation to compute m within the exact algorithm of
    Jenkins and Spano (2015), even when t is not too small
 - Calls "NormalApprox.py"

NormalApprox.py
 - Python code to compute Normal approximation for m within the exact algorithm
    of Jenkins and Spano (2015)    

ComputeKDEs.py
 - Python code for computing kernel density estimates
 - Reads in "qSamples/final_samples_cuda.txt"
 - Produces "DensityEstimates/kdes_fly1217"

MCMC_simXX.py	XX=0,55,11
 - Python code for generating simulated data with s=0.0, s=5.5, or s=11.0 then
    sampling from the posterior distribution via MCMC
 - Uses "informative" prior for initial allele frequencies 
 - Calls "WFwithSelection-ApproxAll.py"
 - Reads in "DensityEstimates/kdes_fly1217" and "Data/fly-n.txt"
 - Produces "Chains/MCMC_simXX_ssX.txt"

MCMC_simXX_unif.py	XX=0,55,11
 - Python code for generating simulated data with s=0.0, s=5.5, or s=11.0 then
    sampling from the posterior distribution via MCMC
 - Uses Uniform prior for initial allele frequencies 
 - Calls "WFwithSelection-ApproxAll.py"
 - Reads in "DensityEstimates/kdes_fly1217" and "Data/fly-n.txt"
 - Produces "Chains/MCMC_simXX_unif_ssX.txt"

MCMC_fly_neutral.py
 - Python code for sampling from the posterior distribution for the neutral fly 
    data using MCMC
 - Calls "WFwithSelection-ApproxAll.py"
 - Reads in "DensityEstimates/kdes_fly1217" and "Data/fly-n.txt"
 - Produces "Chains/MCMC_fly1n_output.txt"

MCMC_fly_select.py
 - Python code for sampling from the posterior distribution for the hypoxic fly
    data using MCMC
 - Calls "WFwithSelection-ApproxAll.py"
 - Reads in "DensityEstimates/kdes_fly1217", "Data/fly-n.txt", and 
    "Data/fly-s.txt"
 - Produces "Chains/MCMC_fly1s_output.txt"

MCMC_fly_windows.py
 - Python code for sampling from the posterior distribution for sliding windows 
    along the 2L chromosome of the hypoxic fly data using MCMC
 - Calls "WFwithSelection-ApproxAll.py"
 - Reads in "DensityEstimates/kdes_fly1217", "Data/fly2Ln.txt", 
    "Data/nlocs2L.txt", "Data/fly2Ls.txt", and "Data/slocs2L.txt"
 - Produces "Chains/MCMC_fly2L_windows.txt"

figure1.R
 - R code for generating histograms of the posterior of s for the simulations
    that use an informative prior for the initial allele frequencies 
 - Reads in "Chains/MCMC_simXX_ssX.txt" for XX=0,55,11 and X=1,2,3,4
 - Produces "Figures/figure1.pdf"

figure2.R
 - R code for generating histograms for the posterior of s for the simulations
    that use an informative prior for the initial allele frequencies 
 - Reads in "Chains/MCMC_simXX_unif_ssX.txt" for XX=0,55,11 and X=1,2,3,4
 - Produces "Figures/figure2.pdf"

figure3.R
 - R code for generating histograms for the posterior of s for the fly data
 - Reads in "Chains/MCMC_fly1s_output.txt" and "Chains/MCMC_fly1n_output.txt"
 - Produces "Figures/figure3.pdf"

figure4.R
 - R code for computing and plotting 99% credible intervals for s for sliding
    windows along the 2L chromosome of the hypoxic fly population
 - Reads in "Data/slocs2L.txt" and "Chains/MCMC_fly2L_windows.txt"
 - Produces "Figures/figure4.pdf"


Output:

qSamples/final_samples_cuda.txt
 - Text file containing simulated values of q from "WrightFisherSampling.cu"

DensityEstimates/kdes_fly1217
 - Pickle file containing the kernel density estimates from "ComputeKDEs.py"

Chains/MCMC_simXX_ssX.txt	XX=0,55,11	X=1,2,3,4
 - Text file containing Markov chains for s for simulated data and an 
    informative prior for the initial allele frequencies 

Chains/MCMC_simXX_unif_ssX.txt	XX=0,55,11	X=1,2,3,4
 - Text file containing Markov chains for s for simulated data and a 
    Uniform prior for the initial allele frequencies 

Chains/MCMC_fly1n_output.txt
 - Text file containing the Markov chain for (s,q) for the neutral fly data

Chains/MCMC_fly1s_output.txt
 - Text file containing the Markov chain for (s,q) for the hypoxic fly data

Chains/MCMC_fly2L_windows.txt
 - Text file containing Markov chains for s for sliding windows along the 
    2L chromosome of the hypoxic fly data 

Figures/figureX.pdf	X=1,2,3,4
 - PDF file containing Figure X from the paper

