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

Note that currently only very specific genome references/assemblies for hg37 and hg38 are supported.
Please have a look at the [Hg37ChrConverter.cpp](src/Hg37ChrConverter.cpp) and [Hg38ChrConverter.cpp](src/Hg38ChrConverter.cpp) files, in which the sets of allowed chromosome names are (yet) hard-coded.

The output is a BED file, which means the start and end positions are 0-based, and left-inclusive, right-exclusive. The columns are:

| Column | Description               | Format                                                     |
|--------|---------------------------|------------------------------------------------------------|
| 1      | chromosome identifier     | string w/o \{:space:, :tab:, :newline, :carriage_return:\} |
| 2      | 0-based start (inclusive) | `\d+`                                                      |
| 3      | 0-based end (exclusive)   | `\d+`                                                      |
| 4      | ???                       | `\d+(,\d+)*`                                               |
| 5      | ???                       | `.\|\d+(,\d+)*`                                            |
| 6      | ???                       | `($chrName:$position)(;($chrName:$position))+`             |
| 7      | ???                       | `#`: missing BP information; `.`: ???                      |
| 8      | ???                       | ibd.                                                       |
| 9      | ???                       | ibd.                                                       |

For column 6.

* The `$chrName` is the chromosome name as given in the input BAM file.
* 

[//]: # (TODO describe the other columns completely)
[//]: # (TODO Describe the format of column 6, which has a more complex format)

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

Note that `sophiaMref` uses a lot of memory (e.g. 600 GB), but usually will be only used for generating the reference files for a new genome assembly (which are currently hardcoded anyway).

`sophiaMref` expects that the input BED files match the filename pattern `.*/$realPidName.{1}$version.+`.
The `$version` is the value provided by the `--version` parameter that is only used to delimit the right end of the PID name.
For instance `/path/to/somePid1_35.1_bps.tsv.gz` would be a valid filename for the version `35.1` and the extracted PID will be `somePid`.

[//]: # (TODO: What will be done with the extracted PID?)


## Runtime Dependencies

The only dependency is Boost 1.82.0 (currently). E.g. you can do

```bash
conda create -n sophia boost=1.82.0 gtest=1.14.0
```

`gtest` is only needed, if you run the tests (which is included in the default `make all` target).

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

#### Running the tests

Currently, there are only very few tests that were added to the legacy code. To run them do

```bash
make test develop=true boost_lib_dir=$boost_lib_dir
```

This will create a `./testRunner` binary in the root directory and execute it.
The `testRunner` is newer linked statically.
It uses the gtest_main library and supports a number of command line options.
See `testRunner --help` for details.

#### Chromosome names

For now, reference genomes are hard-coded, which is why this information is given here in the development section.

Chromosome names must not contain the following characters, because these characters are used as separators:

* `\t`, `\n`, `\r`, ` ` (whitespace)
* `,` (comma)
* `;` (semicolon)
* `|` (pipe)
* `(`, `)` (parentheses)

The `:` (colon) symbol can be used.

## Changes

* 35.1.0 (upcoming)
  * Minor: hg38 support 
    * Added `--assemblyname` option, defaulting to "hg37" when omitted (old behaviour)
    > WARNING: hg38 support was not excessively tested. In particular, yet hardcoded parameters may have to be adjusted. Furthermore, the runtime will be longer than for hg37 and also hg37 runtime has increase.
  * Minor: Build system
    * Use `make` as build system. `Release_*` directories are removed
    * Allow static building with `make static=true boost_lib_dir=/path/to/boost/lib`
    * Allow development build with `make develop=true`
    * Build `sophiaMref`
    * Build `testRunner` for running unit tests
  * Patch: Code readability improvements, documentation, `.editorconfig` file, and `clang-format` configuration
  * Patch: Added few unit tests.
  * Patch: A `README.md` file that is worth its name and contains first documentation about the usage of the SOPHIA tools.

* 9e3b6ed
  * Last version in [bitbucket](https://bitbucket.org/compbio_charite/sophia/src/master/)