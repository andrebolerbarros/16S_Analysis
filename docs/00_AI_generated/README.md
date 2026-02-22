### Prompt

> Prepare the notebooks for the analysis of 16s amplicon sequencing data, using R tools.
> 
> * The structure should be GitHub like, with a README.md file. Then, three separate folders - data, docs and workflow.
> * The pipeline should use DADA2 and ideally phyloseq.
> * It should cover the ASV determination, negative & positive controls processing, diversity analysis (for alpha diversity, you should include calculations with bias correction), and differential abundance.
> * Differential abundance should be done, ideally, with ANCOM-BC2 or Dirichlet regression. Do not use DESeq2.
> * Be sure to organize the DADA2 steps considering continuous quality score vs binned quality scores.
> * The folder already contains the rFigaro repository, which you have created. Please, use it in the workflow as well
> 
> Sources that you can get inspiration from:
> 
> * https://benjjneb.github.io/dada2/tutorial.html
> * https://joey711.github.io/phyloseq/
> * https://f1000research.com/articles/5-1492
> * https://microbiome.github.io/OMA/docs/devel/
> * https://f1000research.com/articles/13-369
> 
> Feel free to run this process iteratively, up to 20 iterations.