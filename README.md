# NSF BigData project Parallel Bioinformatics Library

New name: **Bliss** ??

Please take a look at our [Wiki](https://bitbucket.org/AluruLab/pbil/wiki/Home).

## Building

Build requirements are:

`cmake`, `boost_log`, `boost_system`, `boost_thread`, `boost_program-options`


Building:

```sh
git clone git@bitbucket.org:AluruLab/pbil.git
git submodule init
git submodule update

mkdir build
cd build
cmake ../
make
```

Running the tests:

```sh
ctest -T Test
```


### Configuring Eclipse:

Cmake typically uses a out-of-source build.  to generate eclipse compatible `.project` and `.cproject` files, supply
```
-G"Eclipse CDT4 - Unix Makefiles"
```
to cmake.

install ptp, egit, cmake ed.

