# SOPHIA Tool for Structural Variation Calling

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
  * Last version in [bitbucket]](https://bitbucket.org/compbio_charite/sophia/src/master/)