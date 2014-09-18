# Version

This software was last updated 7/15/2014. Questions, suggestions, comments, etc? 
 
Send emails to nick.zhangyuan@gmail.com  

# Installation

1. Clone the repository:  

  `git clone git@github.com:zhangy72/MetaDomain.git`

2. g++ compiler is required in your Unix system. To install component bin files of MetaDomain, run the Makeme file using the following command:  

  `make`

3. Make sure MetaDomain.sh is executable in your environment. If not you can use the command:  

  `chmod 755 MetaDomain.sh`


# Run MetaDomain

The installation will generate two bin files: hmmscore and generate_domain_expression. The former one is to align short reads to query Pfam domain and the latter is to classify the Pfam domain based on the alignment result. Make sure these two bin files are placed in the same folder with the pipeline of MetaDomain, MetaDomain.sh.

To run MetaDomain pipeline, use the following command:  

`./MetaDomain.sh -m <Pfam HMM file> -f <fasta file> -o <output file> [other options]`

Other possible options:  
```
-h:  show this message
-g:  specify gamma (alignment score threshold rate, in the range of [0,1], default: 0.6)
-n:  specify aligned read number threshold (default: 20)
-c:  specify domain coverage threshold (default: 0.3)
```

The hmm file should be in HMMER3.0's hmm file format. These files can be downloaded from Pfam (http://pfam.xfam.org/). The nucleotide sequence file should be in fasta format.
 

# Output

The output of MetaDomain shows the model length, aligned read number, domain coverage and whether it is transcribed. Here is an example:  

`PF00411 110 19 0.881818 transcribed`

From the output we can tell that the length of PF00411 is 110. There are 19 reads aligned to this domain. The domain coverage is 88.18%. Based on the user-specified thresholds, this domain is transcribed.


# Reference MetaDomain  

MetaDomain can be referenced as:   

<a href="http://psb.stanford.edu/psb-online/proceedings/psb12/zhang-y.pdf">Zhang Y, Sun Y (2012) MetaDomain: a profile HMM-based protein domain classification tool for short sequences. Pac Symp Biocomput: 271–282.</a>


# License

Copyright (C) 2014, Yuan Zhang and Yanni Sun. 

You may redistribute this software under the terms of GNU GENERAL PUBLIC LICENSE.

