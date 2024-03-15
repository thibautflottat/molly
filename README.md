# _molly_&mdash;read `xtc` files, fast

> **WARNING:** This library is unfinished and has not been tested to a
> sufficient degree. 
> 
> Please don't use it for any critical analyses. This repository cannot be
> considered as public, yet. Please inquire whether its address or contents may
> be shared on a case-per-case basis untill otherwise noted.
>
> This library is currently purposefully without a license to prevent
> dependence on its contents.

A reader for the Gromacs [xtc file format][xtc] implemented in pure Rust.

_molly_ is both a Rust library and exposes a set of bindings that allow access
to its functions from Python.

## installation

### As a library

To use _molly_ in a Rust project, add this repository to the dependencies in your `Cargo.toml`.

```toml
[dependencies]
molly = { git = "https://git.sr.ht/~ma3ke/molly" }
```

### The examples

A number of useful example programs can be found in the `examples` directory.
Some of these can be used to benchmark against other xtc reader
implementations, or to create test files.

> **NOTE:** I'm leaving these here for the moment, but ultimately, I will
> remove or fundamentally change many of these examples.

In order to access these, clone the repository and build them.

```console
git clone https://git.sr.ht/~ma3ke/molly
cd molly
cargo build --release --examples
target/release/examples/<name> [options]
```

Or directly run them them using 

```console
cargo run --release --example <name>
```

### As a Python module

A Rust compiler is required for building the Python bindings. 


To install the module, clone the repository, go into the bindings directory,
and install it using `pip`.

```console
git clone https://git.sr.ht/~ma3ke/molly
cd molly/bindings/python

# Perhaps you want to use/create a virtual environment first.
python3 -m venv venv
source venv/bin/activate

pip3 install .
```

## benchmarks

There's some benchmark scripts lying around here. May place them into a neat table at a later point. For now, many things are still changing, so there is no sense in making any hard conclusions about the performance.

Nevertheless, it looks like _molly_ is around 2&times; faster than
[_xdrf_][xdrf] (the widely-used Gromacs implementation), and around
4&times; faster than the [_chemfiles_ implementation][chemfiles].

[xtc]: https://manual.gromacs.org/current/reference-manual/file-formats.html#xtc
[xdrf]: https://gitlab.com/gromacs/gromacs/-/blob/d8d6543db04563cb15f71c90ffb5ed2fda092bce/src/gromacs/fileio/xdrf.h
[chemfiles]: https://chemfiles.org/

---

Marieke Westendorp, 2024
