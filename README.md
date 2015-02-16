# BLISS

Name can stand for:

- **B**asic **Li**fe **S**ciences **S**ubroutines
- **B**ioinformatics **Li**brary for **S**hort **S**equences

Please take a look at our [Wiki](https://bitbucket.org/AluruLab/bliss/wiki/Home).

## Build

### Dependencies

Required:

- `g++` (version 4.7+)
- `cmake` (version 2.6+)


Optional libraries are:

`boost_log`, `boost_system`, `boost_thread`, `boost_program-options`

These are only needed if you intend to turn on boost log engine.


Optional tools include:

- `ccmake` (for graphical cmake configuration)
- `perl`, and perl packages `Term::ANSIColor`, `Getopt::ArgvFile`, `Getopt::Long`, `Regexp::Common` (for g++ error message formatting)


### Getting the source

```sh
git clone git@bitbucket.org:AluruLab/bliss.git
cd bliss
git submodule init
git submodule update
```

### Configuring for build

```sh
mkdir bliss-build
cd bliss-build
cmake ../bliss
```

alternatively, instead of `cmake ../bliss`, you can use

```sh
ccmake ../bliss
```


### CMake Parameters

The following are important parameters:

- `CMAKE_BUILD_TYPE`:  defaults to `Release`.
- `ENABLE_TESTING`:  `On` allows `BUILD_TEST_APPLICATIONS` to show, which enables building the test applications
- `BUILD_EXAMPLE_APPLICATION`: `On` allows applications in the `examples` directory to be built
- `LOG_ENGINE`: chooses which log engine to use.
- `LOGGER_VERBOSITY`: chooses the type of messages to prin.

- `ENABLE_SANITIZER`: turns on g++'s address or thread sanitizer.  Use `SANITIZER_STYLE` to configure.  This is for debugging
- `ENABLE_STLFILT`:  turns on g++ error message post processing to make them human readable.  Control verbosity via `STLFIL_VERBOSITY`

It is highly recommended that `ccmake` be used until you've become familiar with the available CMake options and how to specify them on the commandlinie.


### Compiling

```sh
make
```

### Running the tests

```sh
ctest -T Test
```

or 

```sh
make test
```


### Building the documentation

```sh
make doc
```

### Using Bliss
Please see  [Wiki](https://bitbucket.org/AluruLab/bliss/wiki/Home).



## Configuring Eclipse:

Cmake typically uses a out-of-source build.  to generate eclipse compatible `.project` and `.cproject` files, supply
```
-G"Eclipse CDT4 - Unix Makefiles"
```
to cmake.

Recommend that `ptp`, `egit`, and `cmake ed` also be installed.