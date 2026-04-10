# ExopherGeneExpressionProfiling
Driscoll Lab at Rutgers University
Molecular Biology, Biochemistry, Genetics, functional genomics and bioinformatics
For additional information and consultation regarding the code, please email me, Nelson Mejia, Mark Saba, Durga, Julissa, and James Saba who is our computer scientist.  
This repository contains the R scripts and Python code to analyze two gene expression microarray profiles. 
The first Dataset is titled "A Fasting-Responsive Signaling Pathway that Extends Life Span in C. elegans"  DOI: 10.1016/j.celrep.2012.12.018 NCBI GEO Accession Number GSE27677 
Harvald E, Sprenger R, Dall K..Multi-omics Analyses of Starvation Responses Reveal a Central Role for Lipoprotein Metabolism in Acute Starvation Survival in C. elegans
Cell Systems, 2017; 5, 38-52.e4 GEO accession number; GSE98919
The second Dataset is titled "The Mediator Subunit MDT-15 Confers Metabolic Adaptation to Ingested Material" DOI: 10.1371/journal.pgen.1000021 MNCBI GEO Accession Number GSE9720
Each dataset was analyzed slightly differently given the microarrays used in each study, so for each script, the process of normalization, transformation, and final processing varies, so please pay close attention to which script is being used and execute the line of code with caution. 
This repository is still under construction, so please check back for updates!



of note,
Please follow up with these crucial updates!

High-impact README improvements (bioinformatics / differential expression)
Even without seeing it yet, these are the most common upgrades that materially help:

Clear purpose + scope (top 10 lines)

1–2 sentence summary of what the folder/pipeline does (inputs → method → outputs).
Explicitly state organism, tissue/cell type, and comparison being tested (if known).
Reproducibility section

Exact toolchain: R/Python version, key packages + versions, OS.
A renv.lock / environment.yml / requirements.txt mention and how to restore.
Inputs/Outputs contract

Bullet list of required input files with examples (counts matrix format, metadata columns, reference genome annotation).
A table of outputs (DE table columns, normalized counts, plots) + where they are written.
One “Quickstart” command path

The minimal commands to run end-to-end, ideally copy/paste:
install deps
run
expected output location
Method details (brief but precise)

Which DE method (DESeq2 / edgeR / limma-voom / Seurat, etc.)
Filtering thresholds, normalization, multiple testing correction (FDR), and model formula.
Troubleshooting + FAQ

Common errors (missing sample IDs, non-integer counts for DESeq2, batch effect confusion).
How to validate the design matrix / contrasts.
Project structure

Short tree of directories and what each contains.
Citations

Cite the main method paper(s) and any relevant dataset references.
If you have a DOI / preprint / manuscript, link it.
Make it skimmable

Use consistent headings, add a TOC, keep paragraphs short, add code blocks.
Add a “Results sanity checks” section

MA plot / volcano expectations
sample clustering/PCA
how to interpret log2FC and adjusted p-values
