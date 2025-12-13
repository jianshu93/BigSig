# Large-scale Sequence Search with BItsliced Genomic Signature Index (BIGSIG)

This crate implements **BI**tsliced **G**enomic **SIG**nature Index (BIGSI) for fast sequence search based on kmers or minimizers. We name this software BIGSIG.

Some key functionalities:
1. Use [xxh3](https://crates.io/crates/xxh3) to suport aarch64 and x86-64 platforms;
2. Use [needletail](https://crates.io/crates/needletail) for fast and compressed fasta/fastq file processing;
3. 2-bit nucleitide sequence representation via [NtHash](https://crates.io/crates/nthash); 
4. use [simd-minimizers](https://crates.io/crates/simd-minimizers) for SIMD support minimizer computation
5. BIGSI idea from [here](https://crates.io/crates/bigsi_rs)

## Install
### Pre-compliled binary

```bash
wget https://github.com/jianshu93/BigSig/releases/download/v0.3.1/bigsig_Linux_x86-64_v0.3.1.zip
unzip bigsig_Linux_x86-64_v0.3.1.zip
chmod a=x ./bigsig
./bigsig -h
```


### Complile from source
```bash
git clone https://github.com/jianshu93/BigSig.git
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
  identify   Identify reads based on probability
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

An example to build and query BigSig database
```bash
bigsig construct -r ref_file_example.txt -b test -k 31 -mv 21 -s 10000000 -n 4 -t 24
bigsig identify -b test.mxi -q ./test_data/test.fastq.gz -n output -t 24 --high_mem_load

```

input format (gzippped or not for the file):

```bash
### ref_file_example.txt
Salmonella_enterica_salamae_VIII_74-1880	./refs/74_1880.fasta
Salmonella_enterica_salamae_b_Ar0560	./refs/Ar0560.fasta
Salmonella_enterica_houtenae_Type	./refs/NCTC10401.fasta
Salmonella_bongori_Type-NCTC12419	./refs/S_bongoriType_NCTC12419.fasta
```


## Results
With the default settings BigSiq will report reference sequences that share >35% of their k-mers with the query. Here is the output of a query with SRA reference SRR4098796 (L. monocytogenes lineage I) as query:
```
@m84137_250807_202000_s4/244911914/ccs	Ecoli_CE98,Ecoli_224,Ecoli_OB20	29	891	reject	3
@m84137_250807_202000_s4/244911147/ccs	no_significant_hits	0	4226	reject	0
@m84137_250807_202000_s4/256971887/ccs	Ecoli_CE98	152	762	accept	1
@m84137_250807_202000_s4/245240131/ccs	Ecoli_CE98,Ecoli_224,Ecoli_tEPEC,Ecoli_OB20	27	1454	reject	4
@m84137_250807_202000_s4/245503310/ccs	Ecoli_CE98	503	510	accept	1
@m84137_250807_202000_s4/245765017/ccs	Ecoli_224	42	1578	accept	1
@m84137_250807_202000_s4/244846078/ccs	Ecoli_CE98,Ecoli_OB20	20	1786	reject	2
```
Interpretation:

The identify command writes one line per read (or read pair) with six tab-separated fields: the read ID, an assignment/status (either a single reference, a comma-separated list of references, or special values like no_hits, no_significant_hits, or too_short), the number of supporting kmers/minimizers for the top hit(s), the total number of k-mers/minimizers used from that read after masking/down-sampling, a decision flag (accept or reject), and finally the number of tied top hits. Conceptually, accept means “this line is usable as-is for downstream counting”: either the read has a unique, statistically significant best match to one reference, or it has no hits / is too short but is still a clean, unambiguous outcome. In contrast, reject marks reads where the evidence is ambiguous or unreliable — either multiple references tie for the top score, or no hit passes the false-positive correction — so these are grouped together as “reject” in the summary counts rather than being credited to any specific reference.

BIGSIG models and controls Bloom-filter false positives for each reference in the index. For every reference, BIGSIG computes its theoretical Bloom filter false positive rate from the filter length, number of hash functions, and the number of kmers/minimizers stored for that reference, using the standard approximation for Bloom filters. During classification, the algorithm asks: “If this reference were not truly present, how many k-mer hits would I expect to see just from Bloom filter false positives?” and evaluates this with a binomial model. The --fp_correct parameter sets a per-read significance threshold on this probability (as −log10 p; the default 3.0 corresponds to p = 10⁻³). A read is marked accept only if its top hit has more supporting kmers/minimizers than expected under the false-positive-only model at this threshold; otherwise the read is labeled reject (including ambiguous multi-hit cases). In practice, BIGSIG indices are typically configured so the per-k-mer Bloom false positive rate is on the order of 10⁻³ or lower, and because each read contributes many independent kmers/minimizers, the effective per-read false positive rate is usually orders of magnitude smaller than the underlying Bloom filter rate.




This crate is largely inspired the following papers:

## Reference
1. Bradley, Phelim, et al. "Ultrafast search of all deposited bacterial and viral genomic data." Nature biotechnology 37.2 (2019): 152-159.
2. Bingmann, Timo, et al. "COBS: a compact bit-sliced signature index." String Processing and Information Retrieval: 26th International Symposium, SPIRE 2019, Segovia, Spain, October 7–9, 2019, Proceedings 26. Springer International Publishing, 2019.