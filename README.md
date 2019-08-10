# GmicR
Identifies biomarkers and genes of interest by combining WGCNA, Cell Signatures, and Bayesian network learning to generate a Gene-Module Immune-Cell network

Are you a fan of CD4 T cells? Do B cells really interest you? Are you curious about how fibroblasts or macrophages may be influenced by gene modules in the lungs? These are the types of relationships that can be learned and visualized by this package from bulk RNAseq data .

The goal behind this tool is to provide biomedical researchers, like myself, with a way to come up with 
new molecules of interests or biological functions to explore. WGCNA module eigengenes and xCell immune cell signatures are
merged, relationships can be visualized in the form of a Bayesian network learned using the bnlearn package. 

Inverse relationships can also be highlighted!

This is my first repository on GitHub! 

I thank the community in advanced for helping to refine this package into a tool that may help advance biomedical research.

    # Installation
    install.packages("devtools")
    library(devtools)
    install_github("Rvirgenslane/GmicR")

