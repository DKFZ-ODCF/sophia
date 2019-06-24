
# Runtime Dependencies

The only dependency is Boost 1.70.0 (currently). E.g. you can do

```bash
conda create -n sophia boost=1.70.0
```

# Building

## Build-time Dependencies

* g++ >= 7
* Boost 1.70.0

## Dynamic Build

With Conda you can do

```bash
conda install gxx_linux-64 boost boost=1.70.0
```

to create an environment to build the `sophia` and `sophiaAnnotate` binaries.

To build you need to do

```bash
cd Release_sophia
build-sophia.sh

cd ../Release_sophiaAnnotate
build-sophiaAnnotate.sh
```


## Static Build

If you want to compile statically you need to install glibc and boost static libraries (not possible with Conda in the moment) and do

```bash
cd Release_sophia
STATIC=true build-sophia.sh

cd ../Release_sophiaAnnotate
STATIC=true build-sophiaAnnotate.sh
```
