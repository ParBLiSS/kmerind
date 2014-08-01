# BLISS

Name can stand for:

- **B**asic **Li**fe **S**ciences **S**ubroutines
- **B**ioinformatics **Li**brary for **S**hort **S**equences

Please take a look at our [Wiki](https://bitbucket.org/AluruLab/bliss/wiki/Home).

## Building

Build requirements are:

`cmake`, `boost_log`, `boost_system`, `boost_thread`, `boost_program-options`


### Cloning:

```sh
git clone git@bitbucket.org:AluruLab/bliss.git
cd bliss
git submodule init
git submodule update
```

### Building

```sh
mkdir bliss-build
cd bliss-build
cmake bliss-src
make
```

### Running the tests

```sh
ctest -T Test
```

### Building the documentation

```sh
make doc
```


## Configuring Eclipse:

Cmake typically uses a out-of-source build.  to generate eclipse compatible `.project` and `.cproject` files, supply
```
-G"Eclipse CDT4 - Unix Makefiles"
```
to cmake.

install ptp, egit, cmake ed.