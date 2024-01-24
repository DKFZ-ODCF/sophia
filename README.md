# SOPHIA Tool for Structural Variation Calling

SOPHIA is a Structural Variant(SV) detection algorithm based on the supplementary alignment(SA) concept of the aligner BWA-mem, combined with filters based on expert-knowledge to increase specificity. 

Currently, SOPHIA only is optimized for the hg37 assembly of the human genome.
It uses a large panel of normals for filtering artifacts (most often due to mapping difficulties) and common SVs in the germline.
The parameters for filtering results are hand-tuned against the clinical gold standard FISH of V(D)J rearrangements.
Results from the hand-tuned parameter set were tested against hallmark findings from disease datasets where hallmark SVs were known (CDKN2A in various TCGA datasets, EGFR in TCGA-GBM, GFI1B, MYCN and PRDM6 in ICGC-PEDBRAIN-MB etc.) 

For a detailed description of the algorithm, please refer to Umut Topraks's dissertation at https://doi.org/10.11588/heidok.00027429, in particular chapter 2. Section 2.2.1 describes the method in more details.

SOPHIA is a very fast and resource-light algorithm. 
It uses 2GB RAM, 2 CPU cores and runs in ~3.5 hours for 50x coverage WGS, and can detect variants with a single pass of the input BAM files. No local assembly is done.

> This is a fork of the original [SOPHIA](https://bitbucket.org/utoprak/sophia/src/master/) bitbucket repository.

Sophia is included in the [SophiaWorkflow](https://github.com/DKFZ-ODCF/SophiaWorkflow) that uses the [Roddy Workflow Management Systems](https://github.com/TheRoddyWMS/Roddy).


### Citing

You can cite Sophia as follows:

    Integrative Analysis of Omics Datasets.
    Doctoral dissertation, German Cancer Research Center (DKFZ), Heidelberg.
    Umut Toprak (2019).
    DOI 10.11588/heidok.000274296

### Tools

For instructions and short explanations of on commandline parameters, invoke each tool with `--help`.

#### `sophia` Tool

The main tool for SV calling.
`sophia` reads in a position-sorted BAM file and outputs a list of SVs breakpoints in a tab-separated BED format.

> NOTE: We have only tested BAM created by BWA-MEM.

A call to `sophia` may look like this:

```bash
samtools view -F 0x600 -f 0x001 /yourPositionSorted.bam \
  | sophia --assemblyname hg38 \
           --medianisize 323.0 \
           --stdisizepercentage 21.0 \
           --properpairpercentage 94.32 \
           --defaultreadlength 101 \
           --clipsize 10 \
           --basequality 23 \
           --basequalitylow 12 \
           --lowqualclipsize 5 \
           --isizesigma 5 \
           --bpsupport 3 \
  | gzip --best > out.breakpoints.mref.gz
```

A run on a typical Illumina X10 sample with 30x coverage takes about 5 GB of memory, 2 cores.
In extreme cases (like with chromothripsis) the runtime can jump to 120 hours, but usually is much shorter.

Note that currently only very specific genome references/assemblies for hg37 supported.
Please have a look at the [Hg37ChrConverter.cpp](src/Hg37ChrConverter.cpp) file, in which the sets of allowed chromosome names for "classic_hg37" is hard-coded.

The output is a BED file, which means the start and end positions are 0-based, and left-inclusive, right-exclusive. The 8 columns are:

| Column | Description                                                                                                                                                                                                                                                                                                                                | Format                                                                                                              |
|--------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------|
| 1      | chromosome identifier                                                                                                                                                                                                                                                                                                                      | string w/o \{:space:, :tab:, :newline, :carriage_return:\}                                                          |
| 2      | 0-based start (inclusive)                                                                                                                                                                                                                                                                                                                  | `\d+`                                                                                                               |
| 3      | 0-based end (exclusive). This is just start position + 1.                                                                                                                                                                                                                                                                                  | `\d+`                                                                                                               |
| 4      | Various counts of different types of coverage, specifically (1) pairedBreaksSoft, (2) pairedBreaksHard, (3) mateSupport, (4) unpairedBreaksSoft, (5) unpairedBreaksHard, (6) breaksShortIndel, (7) normalSpans, (8) lowQualSpansSoft, (9) lowQualSpansHard, (10) lowQualBreaksSoft, (11) lowQualBreaksHard, (12) repetitiveOverhangBreaks. | `\d+(,\d+){11}`                                                                                                     |
| 5      | Left and right coverage                                                                                                                                                                                                                                                                                                                    | `\d+,\d+`                                                                                                           |
| 6      | suppAlignmentsDoubleSupport(primary,secondary,mate)                                                                                                                                                                                                                                                                                        | `#`: missing BP information; `.`: no supplementary alignments; `$support(;$support)+`                               |
| 7      | suppAlignments(primary,0,mate)                                                                                                                                                                                                                                                                                                             | ibd.                                                                                                                |
| 8      | Sequences of significant overhangs, if present                                                                                                                                                                                                                                                                                             | `#`: missing BP information; `.`: if the overhang is empty; `${overhangSpec}+` the overhang information (see below) |

Lines starting with `#` are comments and can be ignored.

For columns 6 and 7 the `$support` has a complex format. See `Breakpoint::collapseSuppRange` and `SuppAlignment::print` for details.

For column 8, the overhang information is of the following format

* Each overhang specification has the form `>$id:$left$sequence$right\(\d+\);` where
  * `$id` consisting of the BP index an underscore `_` and a counter
  * `$left` is either empty or a pipe symbol `|`
  * `$sequence` is the sequence of the overhang
  * `$right` is either a pipe symbol `|` or empty (either `$left` or `$right` are a pipe symbol)
  * `\(\d+\)` is the number of times the overhang was observed
  > Developer note: See `Breakpoints::finalizeOverhangs` and `Alignment::printOverhang` for details.

* Columns 6, 7, and 8 are set to `#` if there is missing breakpoint information.

> Developer Note: See `Breakpoint::printBreakpointReport` for details.

#### `sophiaAnnotate` Tool

`sophiaAnnotate` is used for annotating SVs with gene information.
It reads in an output of `sophiaMref` and annotates the SVs in the input file with gene information.

```bash
sophiaAnnotate \
  --tumorresults $tumorSampleFile \
  --controlresults $controlSampleFile \
  --mref $mRef \
  --pidsinmref $pidsInMref \
  --bpfreq $bpFreq \
  --artifactlofreq $artifactLoFreq \
  --artifacthifreq $artifactHiFreq \
  --clonalitystrictlofreq $clonalityStrictLoFreq \
  --clonalitylofreq $clonalityLoFreq \
  --clonalityhifreq $clonalityHiFreq \
  --germlineoffset $germlineFuzziness \
  --defaultreadlengthtumor $tumorDefaultReadLength \
  --defaultreadlengthcontrol $controlDefaultReadLength \
  --germlinedblimit $germlineDbLimit \
  > $outputFile
```



#### `sophiaMref` Tool

The `sophiaMref` tool is used to create a reference file that can be used by `sophiaAnnotate` for annotating SVs with gene information.
Usually, you will only need to run it, if you adapt SOPHIA to a new genome assembly.

`sophiaMref` processes a list of gzipped BED files that were generated with the `sophia` tool.
From these it generates a reference that can be used by `sophiaAnnotate` for annotating structural variants with gene information. 

The file produced by `sophiaMref` is a BED file with the following columns (see `MrefEntry::printBpInfo` for details):

| Column | Description                                              | Format                                                                                                            |
|--------|----------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
| 1      | chromosome identifier                                    | string w/o \{:space:, :tab:, :newline, :carriage_return:\}                                                        |
| 2      | 0-based start (inclusive)                                | `\d+`                                                                                                             |
| 3      | 0-based end (exclusive). This is just start position + 1. | `\d+`                                                                                                             |
| 4      | number of fileIndices given in column 10                 | `\d+`                                                                                                             |
| 5      | number of fileIndicesWithArtifactsRatio                  | `\d+`                                                                                                             |
| 6      | fileIndices.size() / NUMPIDS                             | `\d\.\d+`                                                                                                         |
| 7      | fileIndicesWithArtifactsRatio / NUMPIDS                  | `\d\.\d+`                                                                                                         |
| 8      | average artifacts ratio                                  | `NA` if there are no artifacts ratios; otherwise `\d\.\d+`                                                        |
| 9      | suppAlignments                                           | `.` if there are no supplementary alignments; otherwise the same format as columns 6 and 7 of the `sophia` output |
| 10     | fileIndices                                              | Comma-separated list of file indices (`\d+`)                                                                       |

The "fileIndices" are the 0-based index into the list of gzipped control-BED input files given to `sophiaMref`.

The artifacts ratios are tracked in the moment, but not printed.
Files get an artifacts ratio only if a number of conditions are met (undocumented; see `MrefEntry::addEntry` for details).

Note that `sophiaMref` uses a lot of memory (e.g. 400 GB is a safe choice for human assembly), but usually will be only used for generating the reference files for a new genome assembly (which, currently, are hardcoded anyway).

`sophiaMref` expects that the input BED files match the filename pattern `.*/$realPidName.{1}$version.+`.
The `$version` is the value provided by the `--version` parameter that is only used to delimit the right end of the PID name.
For instance `/path/to/somePid1_35.1_bps.tsv.gz` would be a valid filename for the version `35.1` and the extracted PID will be `somePid`.


## Runtime Dependencies

The only dependency is Boost 1.82.0 (currently) and google-test 1.14.0 for testing. E.g. you can do

```bash
conda create -n sophia boost=1.82.0 gtest=1.14.0 backtrace=20220708 
```

`backtrace` is used to report useful stack traces in case of a crashes
`gtest` is only needed, if you run the tests, which is, however, the default if compile with `make` or `make all`.
To just build the binaries you can do `make binaries`

## Building

> Note that `make` will download one file from [StrTk](https://github.com/ArashPartow/strtk). Furthermore, `make` will download [rapidcsv](https://github.com/d99kris/rapidcsv) for TSV file parsing. If you want to delete an already downloaded file and download it again, run `make clean-all` before the compilation. See the `Makefile` for details.

### Dynamic Build

With [Conda](https://docs.conda.io/) you can do

```bash
conda create -n sophia gxx_linux-64=13 boost=1.82.0 gtest=1.14.0
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

#### Testing Assemblies

There is a new `GenericChrConverter` to handle arbitrary assemblies.
For development, different assemblies can be defined by adding `$assemblyName.tsv` files into the `resources/` directory.

> **NOTE**: For assembly name "classic_hg37" the behaviour of SOPHIA 35.0.0 is used. This is also the default if the `--assemblyname` parameter is omitted.

> **WARNING**: No mechanism to allow this for production has yet been implemented. This is a development feature only! The `GenericChrConverter` implementation is also less performant (yet) than the "classic_hg37" `Hg37ChrConverter`.

The `$assemblyName.tsv` file must be a TSV-separated with 4 columns and a header line with the following fields (names must match exactly):
  * `chromosome`: The chromosome name, which must **not** contain the following characters, because these characters are used as separators:
    * `\t`, `\n`, `\r`, ` ` (whitespace)
    * `,` (comma)
    * `;` (semicolon)
    * `|` (pipe) 
    * `(`, `)` (parentheses)
  
    > **NOTE**: The `:` (colon) symbol is allowed.

  * `size`: The length of the chromosome FASTA sequence in base pairs.
  * `compressedMref`: Whether the chromosome should be considered part of the compressed master-ref set of chromosomes. For instance, for "classic_hg37" these are all chromosomes including autosomes, gonosomes, decoys, unassigned (unplaced/random), EBV, but excluding the mitochondrion (extrachromosomal) and phix (technical). The string is converted to lower-case and matched against the following strings:
    * `true`, `yes`, `t`, `y`, or `1`: The chromosome is part of the compressed master-ref set of chromosomes.
    * `false`, `no`, `f`, `n`, or `0`: The chromosome is not part of the compressed master-ref set of chromosomes.
  * `category`: The following categories are allowed. See `GenericChrConverter::Category` for details. Categories are converted to lower-case and matched against the following strings:
    * `autosomal`: e.g. chr1, chr2, ...
    * `gonosomal`: e.g. chrX, chrY
    * `virus`: e.g. chrEBV
    * `decoy`: e.g. all chromosomes with a _decoy suffix or hs37d5
    * `unassigned`: sequences that belong to normal nuclear genome, but could not be positioned exactly, such as "unplaced", "random", "unlocalized" chromosomes in human assemblies, e.g. chrUn_gl000220
    * `extrachromosomal`: e.g. chrM
    * `technical`: Used for technical reasons, e.g. for calibrating the sequence, such as chrPhiX, lambda
    * `alt`: ALT contigs
    * `hla`: HLA contigs

> **NOTE**: The fact that these classes exist does not mean that they are used in the code.
If you want to know more then, currently, the only documentation of we can offer you for SOPHIA is the source code itself.


#### Running the tests

Currently, there are only very few tests that were added to the legacy code. To run them do

```bash
make test develop=true boost_lib_dir=$boost_lib_dir
```

This will create a `./testRunner` binary in the root directory and execute it.
Note that the `testRunner` is never linked statically.
It uses the `libgtest_main.so` library that supports a number of command line options.
See `testRunner --help` for details.

## Changes

* 35.1.0 (upcoming)
  * Minor: Generic assembly support
    * Added `--assemblyname` option, defaulting to SOPHIA 35.0.0's chromosome converter implementation "classic_hg37" when omitted.
    > WARNING: hg38 support was not excessively tested. In particular, yet hardcoded parameters may have to be adjusted. Furthermore, the runtime will be longer than for hg37 and also hg37 runtime has increase slightly (due to class polymorphism).
  * Minor: Build system
    * Use `make` as build system
    * `Release_*` directories with old build-scripts removed
    * Allow static building with `make static=true boost_lib_dir=/path/to/boost/lib`
    * Allow development build with `make develop=true`
    * Build `sophiaMref` (no build documentation before)
    * Build `testRunner` for running unit tests
  * Patch: Code readability improvements, documentation, `.editorconfig` file, and `clang-format` configuration
  * Patch: Added unit tests.
  * Patch: A `README.md` file that is worth its name and contains first documentation about the usage of the SOPHIA binaries.

* 35 (9e3b6ed)
  * Last version in [bitbucket](https://bitbucket.org/compbio_charite/sophia/src/master/)