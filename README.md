# MetabolicFluxAnalysis

MetabolicFluxAnalysis.jl provides basic methods for metabolic flux analysis based on isotopomer labeling experiments. For this, it provides datastructures to represent the experimental setup, experimental results and atom transition networks, together with an implementation of efficient methods for computing the forward simulation of stationary isotopomer labeling experiments based on the EMU decomposition [1].

In addition to that it provides basic methods for constraint-based modeling and the parameterization of flux spaces.

This package should provide a "batteries-included" environment to do state-of-the-art computational analysis of (stationary) isotopomer labeling experiments. Several experimental datasets from publications are already included (uhm, coming soon..), and it would be great to start a Julia-based community effort for collecting experimental data and devising standardized computational workflows for a thorough and reproducible analysis of labeling experiments.

I feel like the Julia community is the right place to address this issue, and the Julia language and environment provides the optimal tools to perform this task.




Thomas Liphardt, June 2018
