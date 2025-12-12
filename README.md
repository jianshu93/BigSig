# Large-scale Sequence Search with BItsliced Genomic Signature Index (BIGSIG)

This crate implements **BI**tsliced **G**enomic **SIG**nature Index (BIGSIG, also called BIGSI) for fast sequence search based on kmers or minimizers.

Some key functionalities:
1. Use [xxh3](https://crates.io/crates/xxh3) to suport aarch64 and x86-64 platforms;
2. Use [needletail](https://crates.io/crates/needletail) for fast and compressed fasta/fastq file processing;
3. 2-bit nucleitide sequence representation via NtHash; 
4. use simd-minimizers for SIMD support minimizer computation

## Install
```bash
git clone https://gitlab.com/Jianshu_Zhao/bigsig
cd bigsig
### NEON and AVX2 can be used for fast minimizer computation
RUSTFLAGS="-C target-cpu=native" cargo build --release

```

## Usage
```bash

************** initializing logger *****************

Large-scale Sequence Search with BItsliced Genomic Signature Index (BIGSIG)

Usage: bigsig [COMMAND]

Commands:
  construct  Construct a BIGSIG
  query      Query a BIGSIG on one or more fasta/fastq.gz files
  identify   Identify reads based on probability
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version



```
An example to build and query BigSig database
```bash
bigsig construct -r ref_file_example.txt -b test -k 31 -mv 21 -s 10000000 -n 4 -t 24
bigsig query -b ./test.mxi  -q ./test_data/test.fastq.gz 
bigsig identify -b test.mxi -q ./test_data/test.fastq.gz -n output -t 24 --high_mem_load

```
## Results
With the default settings BigSiq will report reference sequences that share >35% of their k-mers with the query. Here is the output of a query with SRA accession SRR4098796 (L. monocytogenes lineage I) as query:
```
SRR4098796_1.fastq.gz	3076072	Listeria_monocytogenes_F2365	0.87	134.25	126	475266
SRR4098796_1.fastq.gz	3076072	Listeria_monocytogenes_SRR2167842	0.40	128.25	122	7831
```
In the first column is the query, the second column shows the number of k-mers/minimizers in the query, the third column displays the reference sequence, the fourth column the proportion of kmers/minimizers in the reference shared with the query, the fifth column displays the average coverage based on k-mers/minimizers that were uniquely matched with this reference, the sixth column is the modus of the coverage based on uniquely matched k-mers and the last column the number of uniquely matched k-mers.



This crate is largely inspired the following papers:

## Reference
1. Bradley, Phelim, et al. "Ultrafast search of all deposited bacterial and viral genomic data." Nature biotechnology 37.2 (2019): 152-159.
2. Bingmann, Timo, et al. "COBS: a compact bit-sliced signature index." String Processing and Information Retrieval: 26th International Symposium, SPIRE 2019, Segovia, Spain, October 7–9, 2019, Proceedings 26. Springer International Publishing, 2019.