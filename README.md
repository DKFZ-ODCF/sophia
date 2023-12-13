# SOPHIA Tool for Structural Variation Calling

SOPHIA is a Structural Variant(SV) detection algorithm based on the supplementary alignment(SA) concept of the aligner BWA-mem, combined with filters based on expert-knowledge to increase specificity. 

Currently, SOPHIA only is optimized for the hg38 assembly of the human genome. 
It uses a large panel of normals for filtering artifacts (most often due to mapping difficulties) and common SVs in the germline.
The parameters for filtering results are hand-tuned against the clinical gold standard FISH of V(D)J rearrangements.
Results from the hand-tuned parameter set were tested against hallmark findings from disease datasets where hallmark SVs were known (CDKN2A in various TCGA datasets, EGFR in TCGA-GBM, GFI1B, MYCN and PRDM6 in ICGC-PEDBRAIN-MB etc.) 

For a detailed description of the algorithm, please refer to Umut Topraks's dissertation at https://doi.org/10.11588/heidok.00027429, in particular chapter 2. Section 2.2.1 describes the method in more details.

SOPHIA is a very fast and resource-light algorithm. It uses 2GB RAM, 2 CPU cores and runs in ~3.5 hours for 50x coverage WGS, and can detect variants with a single pass of the input BAM files. No local assembly is done.

> This is a fork of the original [SOPHIA](https://bitbucket.org/utoprak/sophia/src/master/) bitbucket repository.

Sophia is included in the [SophiaWorkflow](https://github.com/DKFZ-ODCF/SophiaWorkflow) that uses the [Roddy Workflow Management Systems](https://github.com/TheRoddyWMS/Roddy).


### Citing

You can cite Sophia as follows:

    Integrative Analysis of Omics Datasets.
    Doctoral dissertation, German Cancer Research Center (DKFZ), Heidelberg.
    Umut Toprak (2019).
    DOI 10.11588/heidok.000274296

## Runtime Dependencies

The only dependency is Boost 1.70.0 (currently). E.g. you can do

```bash
conda create -n sophia boost=1.70.0
```

## Building

### Build-time Dependencies

* g++ >= 7
* Boost 1.70.0

### Dynamic Build

With Conda you can do

```bash
conda create -n sophia gxx_linux-64=8 boost=1.70.0
```

to create an environment to build the `sophia` and `sophiaAnnotate` binaries.

To build you need to do

```bash
source activate sophia
cd Release_sophia
build-sophia.sh

cd ../Release_sophiaAnnotate
build-sophiaAnnotate.sh
```

Note that the build-scripts are for when you manage your dependencies with Conda.


### Static Build

If you want to compile statically you need to install glibc and boost static libraries (currently, not possible with Conda) and do

```bash
source activate sophia

cd Release_sophia
STATIC=true build-sophia.sh

cd ../Release_sophiaAnnotate
STATIC=true build-sophiaAnnotate.sh
```

## Changes

* 35.1.0 (upcoming)
  * Minor: Nominally added support for hg38 (hg37 support remains)
  * Minor: Added `--assemblyname` option, defaulting to "hg37" when omitted (old behaviour)
    > WARNING: hg38 support was not excessively tested. In particular, yet hardcoded parameters may have to be adjusted.
  * Minor: Build script for `sophiaMref`
  * Patch: Code readability improvements, `.editorconfig` file, and `clang-format` configuration
  * Patch: Improved compilation instructions
  * Patch: Use `namespace::std` to get rid of `std::` noise in the code

* 9e3b6ed
  * Last version in [bitbucket](https://bitbucket.org/compbio_charite/sophia/src/master/)