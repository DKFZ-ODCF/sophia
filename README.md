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

### Tools

* `sophia` - The main tool for SV calling. It takes a BAM file as input and outputs a list of SVs in mref format.
* `sophiaAnnotate` - Tool for annotating SVs with gene information. It reads in an mref file created by `sophiaMref` and annotates the SVs in the input file with gene information.
* `sophiaMref` - The `sophiaMref` tool processes a list of gzipped control bed files and generates a reference that can be used by `sophiaAnnotate` for annotating structural variants with gene information.

For instructions on commandline parameters, invoke the tool with `--help`.

## Runtime Dependencies

The only dependency is Boost 1.82.0 (currently). E.g. you can do

```bash
conda create -n sophia boost=1.82.0
```

## Building

> Note that `make` will download one file from [StrTk](https://github.com/ArashPartow/strtk). If you want to delete an already downloaded file and download it again, run `make clean-all` before the compilation. See the `Makefile` for details.

### Dynamic Build

With [Conda](https://docs.conda.io/) you can do

```bash
conda create -n sophia gxx_linux-64=13 boost=1.82.0
```

to create an environment to build the SOPHIA binaries.

To build you need to do

```bash
source activate sophia
make -j 4
```

The binaries will be located in the top-level directory.


### Static Build

If you want to compile statically you need to install glibc and boost static libraries. Conda does not provide a statically compiled version of boost, though. Please refer to the [boost documentation]() for detailed instructions or in case of problems.

Shortly, to install boost statically (here without installing it system-wide) you need to do

```bash
./bootstrap.sh --prefix=build/
./b2 --prefix=build/ link=static runtime-link=static cxxflags="-std=c++11 -fPIC"
boost_lib_dir=$PWD/stage/lib
```

> NOTE: This requires that you have a C++11-compatible compiler installed.

After that you can do:

```bash
source activate sophia
cd "$repoRoot"
make -j 4 static=true boost_lib_dir=$boost_lib_dir all
```

### Development build

The development build produces non-optimized binaries (`-O0`) with debug symbols:

```bash
source activate sophia
cd "$repoRoot"
make -j 4 static=true boost_lib_dir=$boost_lib_dir develop=true all
```

## Changes

* 35.1.0 (upcoming)
  * Minor: Multiple assembly support 
    * Added `--assemblyname` option, defaulting to "hg37" when omitted (old behaviour)
    * Added support for hg38
    > WARNING: hg38 support was not excessively tested. In particular, yet hardcoded parameters may have to be adjusted.
  * Minor: Build system
    * Use `make` as build system. `Release_*` directories are removed.
    * Allow static building with `make STATIC=true boost_lib_dir=/path/to/boost/lib`
    * Build `sophiaMref`
  * Patch: Code readability improvements, `.editorconfig` file, and `clang-format` configuration

* 9e3b6ed
  * Last version in [bitbucket](https://bitbucket.org/compbio_charite/sophia/src/master/)

## How to run

### Sophia
Sophia takes a BAM file as input and outputs a list of breakpoints (BPs). The control and tumor BAM files have to be run separately. The outputs are merged in the annotation step.


```bash
samtools view -F 0x600 -f 0x001 {TUMOR_BAM_FILE} | \
  sophia \
  --assemblyname {ASSEMBLY} \
  --medianisize {MEDIAN_SIZE} \
  --stdisizepercentage {STD_PERC} \
  --properpairpercentage {PP_PERC} \
  --defaultreadlength {READ_LENGTH} | gzip --best > {PID-TUMOR}_{SOPHIA_VERSION}.bps.tsv.gz
```

#### Sophia annotation
Sophia annotation takes the output of Sophia and annotates the breakpoints with gene information. The output is a BEDPE file.

```bash
sophiaAnnotate \
  --tumorresults {PID-TUMOR}_{SOPHIA_VERSION}.bps.tsv.gz \
  --controlresults {PID-CONTROL}_{SOPHIA_VERSION}.bps.tsv.gz \
  --mref mergedMref_{SOPHIA_VERSION}_{NO_PIDs}_min3.bed.gz \
  --pidsinmref {PID_COUNT_IN_MREF} \
  --defaultreadlengthtumor {TUMOR_READ_LENGTH} \
  --defaultreadlengthcontrol {CONTROL_READ_LENGTH}
```

#### Sophia Master reference (Mref)
Sophia Mref takes a list of BP files from control samples and generates a background reference that can be used by Sophia annotation for annotating structural variants with the population information and filtering out common variants.

```bash
sophiaMref \
  --assemblyname {ASSEMBLY} \
  --gzins sophia_control_bps_paths.tsv \
  --version {SOPHIA_VERSION} \
  --defaultreadlength {READ_LENGTH} \
  --outputrootname mergedMref
```
The `sophia_control_bps_paths.tsv` contains full paths to the BP files from the control samples.

The `mergedMref` data can be further filtered based on the number of samples with a breakpoint (column 4 in the output file). The last column in the output file, which contains the PID indices, could be ignored (see below for details on the output formats).

```bash
module load htslib
awk -F '\t' '$4 > 2' mergedMref_{SOPHIA_VERSION}_{NO_PIDs}.bed | \
  cut -f1-9 | \
  bgzip > mergedMref_{SOPHIA_VERSION}_{NO_PIDs}_min3.bed.gz
```

## Output Format

### Sophia
Sophia outputs a list of breakpoints (BPs) in a BED format file. The columns are as follows:

1. Chromosome name
2. Position of the breakpoint
3. Position of the breakpoint + 1
4. CovType (pairedBreaksSoft,pairedBreaksHard,mateReadSupport,unpairedBreaksSoft,unpairedBreaksHard,shortIndelReads,normalSpans,lowQualSpansSoft,lowQualSpansHard,lowQualBreaksSoft,lowQualBreaksHard,repetitiveOverhangs)
5. leftCoverage,rightCoverage
6. suppAlignmentsDoubleSupport(primary,secondary,mate)
7. suppAlignments(primary,0,mate)
8. significantOverhangs

The supplementary alignments (SA) data structure in `suppAlignmentsDoubleSupport` and in `suppAlignments` have a similar data structure - eg 2:91647217-91647218|(0,1,2?/36);15:72741680|(0,2,2?/36);hs37d5:17985120-17985121_INV|(0,1,2?/36)

The SA data structure starts with `chr` and `position`. If the SA is fuzzy, an `extended position` is added to the coordinates. Then, if SA is inverted, a "_INV" str is added. Within the brackets, counts are added for 'primary support', 'secondary support' and 'mate support'. If the SA is semi-suspicious, a "?" is added and if it is suspicious, a "!" is added. And the number of discordants are added after "/". Based on the "M" in the CIGAR string, "|" is added before or after the coordinate positions.

A similar data structure is used for `suppAlignments` column as well. But the counts are only for 'primary support' and 'mate support'.

### Sophia annotation

### Sophia Master reference (Mref)

Sophia Mref outputs a list of merged breakpoints (BPs) in a BED format file. The columns are as follows:

1. Chromosome name
2. Position of the breakpoint
3. Position of the breakpoint + 1
4. Number of file indices (BP bed files used for generating the Mref databases)
5. Number of file indices with artifact ratios
6. Ratio of the number of file indices to the total number of PIDs (NUMPIDS)
7. Ratio of the number of file indices with artifact ratios to the total number of PIDs (PIDs are counted here if the sum of `artifactTotalRelaxed` and `eventTotalStrict` is greater than > 0. So in some cases, the AF ratio could be zero even when this ratio > 0.)
8. Average of artifact ratios if they exist, otherwise "NA"
9. If there are no supplementary alignments, a "." is added. If there are supplementary alignments, each one is finalized and printed if the total number of supplementary alignments is less than 11 or if the support of the alignment is greater than or equal to 20% of the number of file indices. The printed supplementary alignments are joined with a ";" separator.
10. The file indices are converted to strings and joined with a "," separator.

Artifact ratio: From each BPs, `artifactTotalRelaxed` is counted from the sum of `lowQualSpansSoft,lowQualBreaksSoft,repetitiveOverhangs`. And `eventTotalStrict` is calculated from the sum of `pairedBreaksSoft,pairedBreaksHard,unpairedBreaksSoft`. The artifact ratio is calculated by dividing the `artifactTotalRelaxed` by the sum of `eventTotalStrict` and `artifactTotalRelaxed`.



