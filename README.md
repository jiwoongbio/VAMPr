# VAMPr
VAriant Mapping and Prediction of antimicrobial resistance


## Requirements

1. Perl - https://www.perl.org
2. R - http://www.r-project.org
3. Perl module Bio::DB::Fasta - https://metacpan.org/pod/Bio::DB::Fasta
4. Perl module Statistics::R - https://metacpan.org/pod/Statistics::R
5. R library caret - https://cran.r-project.org/web/packages/caret/index.html
6. R library xgboost - https://cran.r-project.org/web/packages/xgboost/index.html
7. DIAMOND - https://github.com/bbuchfink/diamond
8. Linux commands: sort, wget - https://www.gnu.org/software/wget/


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git. It will take a few seconds.
```
git clone https://github.com/jiwoongbio/VAMPr.git
```


## Usages

* VAMP.pl
```
Usage:   perl VAMP.pl [options] genome.fasta [variant.vcf [sample [...]]] > VAMP.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [$TMPDIR or /tmp]
         -C STR   codon and translation e.g. ATG=M [NCBI genetic code 11 (Bacterial, Archaeal and Plant Plastid)]
         -S STR   comma-separated start codons [GTG,ATG,CTG,TTG,ATA,ATC,ATT]
         -T STR   comma-separated termination codons [TAG,TAA,TGA]
         -L INT   minimum translation length [10]
         -p INT   number of threads [1]
         -e FLOAT maximum e-value to report alignments [10]
         -c FLOAT minimum coverage [0.8]
         -s FILE  output SAM file
         -a FILE  output alignment file
         -A       all variants
```

* VAMP_protein.pl
```
Usage:   perl VAMP_protein.pl [options] protein.fasta > VAMP.txt

Options: -h       display this help message
         -t DIR   directory for temporary files [$TMPDIR or /tmp]
         -p INT   number of threads [1]
         -e FLOAT maximum e-value to report alignments [10]
         -c FLOAT minimum coverage [0.8]
         -s FILE  output SAM file
         -a FILE  output alignment file
         -A       all variants
```

* VAMP_database.pl
```
Usage:   perl VAMP_database.pl [options] KEGG_orthology.txt

Options: -h       display this help message
         -r       redownload data
         -c FLOAT sequence identity threshold [0.7]
         -d       decoy
         -a FILE  additional protein sequence file including orthology, protein and sequence columns
         -p INT   number of threads [1]
```

* VAMP_database.feature.pl
```
Usage:   perl VAMP_database.feature.pl [options]

Options: -h       display this help message
         -r       redownload data
```

* VAMP_feature.pl
```
Usage:   perl VAMP_feature.pl [options] VAMP.txt > VAMP_feature.txt

Options: -h       display this help message
```

* VAMP_fisher.pl
```
Usage:   perl VAMP_fisher.pl [options] phenotype=VAMP.txt[,...] [...] > VAMP_fisher.txt

Options: -h       display this help message
         -g FILE  genotype file
         -c       cluster only
         -o       orthology only
```

* VAMP_fisher.filter.pl
```
Usage:   perl VAMP_fisher.filter.pl VAMP_fisher.txt > VAMP_fisher.filter.txt

Options: -h       display this help message
         -a FLOAT alpha, p-value cutoff
         -D       do not consider drive count
         -O       do not consider odds ratio
```

* VAMP_model.pl
```
Usage:   perl VAMP_model.pl [options] VAMP_model.RData phenotype=VAMP.txt[,...] [...]

Options: -h       display this help message
         -g FILE  genotype file
         -c       cluster only
         -o       orthology only
         -s INT   seed [1]
         -f INT   fold [5]
         -G STR   grid
```

* VAMP_nestedCV.pl
```
Usage:   perl VAMP_nestedCV.pl [options] VAMP_nestedCV.accuracy.txt phenotype=VAMP.txt,... [...]

Options: -h       display this help message
         -g FILE  genotype file
         -c       cluster only
         -o       orthology only
         -s INT   seed [1]
         -f INT   fold [5]
         -F INT   outer fold [10]
         -G STR   grid
         -i FILE  image file
```

* VAMP_model.genotype.pl
```
Usage:   perl VAMP_model.genotype.pl [options] VAMP_model.RData > VAMP_model.genotype.txt

Options: -h       display this help message
```

* VAMP_predict.pl
```
Usage:   perl VAMP_predict.pl [options] VAMP_model.RData sample=VAMP.txt [...]

Options: -h       display this help message
         -c       cluster only
         -o       orthology only
         -p STR   phenotype
         -i FILE  image file
```
