
# Differential expression analysis pipeline

## Overview

This is a basic differential expression analysis pipeline based developed at the 
Life Sciences Computation Centre, University of Melbourne.

It is based around the Tuxedo protocol (Tophat, Cufflinks) and EdgeR/Voom with HTSeq-count.

To run the pipeline you will need the Ruffus library: [http://www.ruffus.org.uk/](http://www.ruffus.org.uk/).

Documentation for the pipeline can be found [here](https://bitbucket.org/jessicachung/rna_seq_pipeline/wiki/Home).

Usage: 
    
    python RNA-seq_pipeline.py
        [-h | --help]
        --opts=<pipeline options files>
        --style=<run | print | flowchart>
        --force=<force this task to run>
        --end=<final task>
        --rebuild=<fromtarget | fromstart>
        --verbose=<0 | 1 | 2>

Example:

    python RNA-seq_pipeline.py \
        --opts=pipeline_config,pipeline_stages_config \
        --style=run \
        --end=fastQC

If arguments are not specified, the default setting listed in pipeline_config.py will be used.

## Running on VLSCI's clusters (i.e. Merri or Barcoo)

Use Python 2.7.5, which you can load with

    module load python-gcc/2.7.5

To use the flowchart option you will need graphviz which you can load with

    module load graphviz
    
If you're using a terminal multiplexer (such as GNU Screen) on merri, some environment variables aren't inherited so you'll have to reload your modules.

    module purge
    module load vlsci/2014-09
