# SOPHIA Tool for Structural Variation Calling

SOPHIA is a Structural Variant(SV) detection algorithm based on the supplementary alignment(SA) concept of the aligner BWA-MEM, combined with filters based on expert-knowledge to increase specificity. 

Currently, SOPHIA only is optimized for the hg37 assembly of the human genome.

> **NOTE**: We are preparing hg38 support. See the [Changes](#changes) section for details.

It uses a large panel of normals for filtering artifacts (most often due to mapping difficulties) and common SVs in the germline.
The parameters for filtering results are hand-tuned against the clinical gold standard FISH of V(D)J rearrangements.
Results from the hand-tuned parameter set were tested against hallmark findings from disease datasets where hallmark SVs were known (CDKN2A in various TCGA datasets, EGFR in TCGA-GBM, GFI1B, MYCN and PRDM6 in ICGC-PEDBRAIN-MB etc.) 

For a detailed description of the algorithm, please refer to [Umut Topraks's dissertation](https://doi.org/10.11588/heidok.00027429), in particular chapter 2. Section 2.2.1 describes the method in more details.

SOPHIA is a very fast and resource-light algorithm. 
It uses 2GB RAM, 2 CPU cores and runs in ~3.5 hours for 50x coverage WGS, and can detect variants with a single pass of the input BAM files. No local assembly is done.

Sophia is included in the [SophiaWorkflow](https://github.com/DKFZ-ODCF/SophiaWorkflow) that uses the [Roddy Workflow Management Systems](https://github.com/TheRoddyWMS/Roddy).


### Citing

You can cite the original version (35) of SOPHIA as follows:

    Integrative Analysis of Omics Datasets.
    Doctoral dissertation, German Cancer Research Center (DKFZ), Heidelberg.
    Umut Toprak (2019).
    DOI 10.11588/heidok.000274296

The code for the original version 35 can be found in the old [SOPHIA repository](https://bitbucket.org/utoprak/sophia/src/master/) bitbucket repository. 

The code here is a fork of that repository. If you use the newer versions here, please also include a reference to this repository in your citation -- in particular if you use SOPHIA for any other reference genome that the [1000 genomes reference](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
). See [Contributors](CONTRIBUTORS.md).

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
           --medianisizes 323.0 \
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

| Parameter                | Description                                                                                                                                                                                      |
|--------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--assemblyname`         | The name of the assembly. See [Reference Genomes / Assemblies](#reference-genomes--assemblies) for details.                                                                                      |
| `--mergedisizes`         | A file with just the median insert size in line 1 and just the standard deviation in line 3. See [Insert Size Distribution](#insert-size-distribution).                                          | 
| `--medianisize`          | The median insert size of the library. See [Insert Size Distribution](#insert-size-distribution).                                                                                                |
| `--stdisizepercentage`   | The standard deviation of the insert size of the library **in percent**. See [Insert Size Distribution](#insert-size-distribution).                                                              |
| `--isizesigma`           | The sigma value for the insert size distribution. See [Insert Size Distribution](#insert-size-distribution).                                                                                     |
| `--properpairpercentage` | Proper pair ratio as a percentage.                                                                                                                                                   |
| `--defaultreadlength`    | Default read length for the technology used in sequencing 101, 151, etc.                                                                                                                                                                         |
| `--clipsize`             | Minimum length of soft/hard clips in the alignment.                                                                                                                                                              |
| `--basequality`          | The minimum base quality for a base to be considered high quality.                                                                                                                               |
| `--basequalitylow`       | If 5 consecutive bases in a split read overhang have lower quality than this strict threshold, it will be low-quality clipped. |
| `--lowqualclipsize`      | Maximum length of a low quality split read overhang for discarding.                                                                                                                                       |
| `--bpsupport`            | Minimum number of reads supporting a discordant contig.                                                                                                                        |

##### Insert Size Distribution

The median insert size and its standard deviation are used for calculating the maximum insert size considered by the algorithm. 

For `--mergedisizes`, the actual value is calculated by $`min(4000.0, mean + isizeSigmaLevel * sd`$).

For `--medianisize`/`--stdisizepercentage`, the actual values is calculated by $`min(4000.0, median + isizeSigmaLevel * median * isizeStdPercentage * 0.01)`$.

If `--mergedisizes` is provided, `--medianisize` and `--stdisizepercentage` are ignored.

If none of these options are provided, instead a maximum insert size of 2000 is assumed.

##### Outputs

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
  --PIDS_IN_MREF $pidsInMref \
  --bpfreq $bpFreq \
  --artifactlofreq $artifactLoFreq \
  --artifacthifreq $artifactHiFreq \
  --clonalitystrictlofreq $clonalityStrictLoFreq \
  --clonalitylofreq $clonalityLoFreq \
  --clonalityhifreq $clonalityHiFreq \
  --germlineoffset $germlineFuzziness \
  --defaultreadlengthtumor $tumordefaultreadlength \
  --defaultreadlengthcontrol $controldefaulreadlength \
  --germlinedblimit $germlinedblimit \
  > $outputFile
```



#### `sophiaMref` Tool

The `sophiaMref` tool is used to create a reference file that can be used by `sophiaAnnotate` for annotating SVs with gene information.
Usually, you will only need to run `sophiaMref`, if you adapt SOPHIA to a new genome assembly.

`sophiaMref` processes a list of gzipped BED files that were generated with the `sophia` tool.
From these it generates a reference that can be used by `sophiaAnnotate` for annotating structural variants with gene information. 

The file produced by `sophiaMref` is a BED file suffixed with the following columns (see `MrefEntry::printBpInfo` for details):

| Column | Description                                              | Format                                                                                                            |
|--------|----------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
| 1      | chromosome identifier                                    | string w/o \{:space:, :tab:, :newline, :carriage_return:\}                                                        |
| 2      | 0-based start (inclusive)                                | `\d+`                                                                                                             |
| 3      | 0-based end (exclusive). This is just start position + 1. | `\d+`                                                                                                             |
| 4      | number of fileIndices given in column 10                 | `\d+`                                                                                                             |
| 5      | number of fileIndicesWithArtifactsRatio                  | `\d+`                                                                                                             |
| 6      | fileIndices.size() / NUM_PIDS                             | `\d\.\d+`                                                                                                         |
| 7      | fileIndicesWithArtifactsRatio / NUM_PIDS                  | `\d\.\d+`                                                                                                         |
| 8      | average artifacts ratio                                  | `NA` if there are no artifacts ratios; otherwise `\d\.\d+`                                                        |
| 9      | suppAlignments                                           | `.` if there are no supplementary alignments; otherwise the same format as columns 6 and 7 of the `sophia` output |
| 10     | fileIndices                                              | Comma-separated list of file indices (`\d+`). Maybe empty string.                                                 |

The "fileIndices" are the 0-based index into the list of gzipped control-BED input files given to `sophiaMref`.

Currently, the artifacts ratios are tracked, but not printed.
Files get an artifacts ratio only if a number of conditions are met (undocumented; see `MrefEntry::addEntry` for details).

Note that `sophiaMref` uses a lot of memory (e.g. 400 GB is a safe choice for human assembly), but usually will be only used for generating the reference files for a new genome assembly (which, currently, are hardcoded anyway).

`sophiaMref` expects that the input BED files match the filename pattern `.*/$realPidName.{1}$version.+`.
The `$version` is the value provided by the `--version` parameter that is only used to delimit the right end of the PID name.
For instance `/path/to/somePid1_35.1_bps.tsv.gz` would be a valid filename for the version `35.1` and the extracted PID will be `somePid`.


## Dependencies

If you built SOPHIA with dynamic libraries, then some libraries are runtime requirements, namely:

  * Boost 1.82.0
  * gtest 1.14.0
  * strtk 0.6.0

These may be runtime dependencies, if you choose a dynamic build.
The static build creates self-contained binaries that do not have any runtime dependencies.

You can install all dependencies for the dynamic build with [Conda](https://docs.conda.io/):

```bash
conda create -n sophia gxx_linux-64=13 boost=1.82.0 gtest=1.14.0
```

## Building

> Note that `make` will download [StrTk](https://github.com/ArashPartow/strtk) for string processing. If you want to delete an already downloaded file and download it again, run `make clean-all` before the compilation. See the `Makefile` for details.

### Dynamic Build

For compilation you additionally need the g++ compiler. So extend the Conda environment a bit:

```bash
conda create -n sophia gxx_linux-64=13 boost=1.82.0 gtest=1.14.0
```

Then you can do

```bash
source activate sophia
make -j 4 all
```

This will also run all unit tests (not many for this branch, though).
The binaries will be located in the top-level directory.

### Static Build

Static building produced 100% self-contained binaries that you can copy to any compatible OS independent of which libraries are installed there.

For static building all dependencies need to be available as static libraries (`.a` files). 
Specifically, libz, libm, glibc, and libstdc++ need to be available as static libraries. 
Fortunately, static versions of these libraries are installed by Conda.
You can use the same build environment as for the dynamic build. 

The only complications are the gtest library and the Boost libraries, both of which are not available as static builds from Conda.

The gtest library is only used for testing, and the `testRunner` binary is never build statically.

To install boost statically (here without installing it system-wide) you need to do the following:

```bash
# Download boost
wget https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.bz2
tar -xjf boost_1_82_0.tar.bz2

# Build b2
cd boost_1_82_0
./bootstrap.sh --prefix=build/

# Build boost
./b2 --prefix=build/ link=static runtime-link=static cxxflags="-std=c++11 -fPIC"
boost_lib_dir=$PWD/stage/lib
```

After that you can do:

```bash
source activate sophia
cd "$repoRoot"
make -j 4 static=true boost_lib_dir=$boost_lib_dir all
```

### Development build

> NOTE: Please adhere to the [Coding Conventions](CodingConventions.md) when developing in SOPHIA.

The development build produces non-optimized binaries (`-O0`) with debug symbols:

```bash
source activate sophia
cd "$repoRoot"
make -j 4 static=true boost_lib_dir=$boost_lib_dir develop=true all
```

### Running the tests

Currently, there are only very few tests that were added to the legacy code. To run them do

```bash
make test develop=true boost_lib_dir=$boost_lib_dir
```

This will create a `./testRunner` binary in the root directory and execute it.
Note that the `testRunner` is never linked statically.
It uses the `libgtest_main.so` library that supports a number of command line options.
See `testRunner --help` for details.


## Reference Genomes / Assemblies

> WARNING: This version of SOPHIA only supports a very specific reference genome for GRCh37/hg37.

The 35 version of SOPHIA was extensively tested on the [1000 genomes reference](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz) with a [phix](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1) sequence added.

The parser for chromosome names is very fragile, and also imprecise (e.g. assuming that a chromosome name starting with 'h' refers to 'hs37d5', the decoy sequences). The best is you create a referenc genome with exactly the chromosome names listed in [resources/hs37d5+phix.tsv](resources/hs37d5+phix.tsv).

## Changes

* 35.0.1
  * Minor: Build system
    * Use `make` as build system
    * `Release_*` directories with old build-scripts removed
    * Allow static building with `make static=true boost_lib_dir=/path/to/boost/lib`
    * Allow development build with `make develop=true`
    * Build `sophiaMref` (no build documentation before)
    * Build `testRunner` for running unit tests
  * Patch: A `README.md` file that is worth its name and contains first documentation about the usage of the SOPHIA binaries and input and output files.
  * Patch: Added unit tests.
  * Patch: Code readability improvements, documentation, `.editorconfig` file, and `clang-format` configuration
  
* 35 (9e3b6ed)
  * Forked from [bitbucket](https://bitbucket.org/compbio_charite/sophia/src/master/)