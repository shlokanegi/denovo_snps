Container to validate denovos by generating a pileup of reads at candidate de novo site in both parents.

The [`dnvval.py`](dnvval.py) script is included (available at `/opt/scripts/dnvval.py` in the container).
It will run samtools pileup for each variant with both mother's and father's BAM file, parse matches and mismatches and annotate VCF records with a true-positive flag if none of the the parents' reads support the variant.

The inputs of `dnvval.py` are:

```
root@0a5dd60eee5f:/home# python3 /opt/scripts/dnvval.py --help
usage: dnvval.py [-h] -mbam MBAM -dbam DBAM -r R -v V [-o O]

optional arguments:
  -h, --help  show this help message and exit
  -mbam MBAM  mother's BAM file (indexed)
  -dbam DBAM  father's BAM file (indexed)
  -r R        reference FASTA file (indexed)
  -v V        VCF file (can be bgzipped)
  -o O        output (annotated) VCF. Will not be bgzipped
```

One new INFO field is added in the output VCF:

- `flag` with the 'TP' label, if the variant was not round in any other the reads in both parents.

To build locally and upload to [quay.io](https://quay.io/shnegi/dvval).
