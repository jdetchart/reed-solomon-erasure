# reed-solomon-erasure
[![Build Status](https://travis-ci.org/darrenldl/reed-solomon-erasure.svg?branch=master)](https://travis-ci.org/darrenldl/reed-solomon-erasure)
[![Build status](https://ci.appveyor.com/api/projects/status/47c0emjoa9bhpjlb/branch/master?svg=true)](https://ci.appveyor.com/project/darrenldl/reed-solomon-erasure/branch/master)
[![codecov](https://codecov.io/gh/darrenldl/reed-solomon-erasure/branch/master/graph/badge.svg)](https://codecov.io/gh/darrenldl/reed-solomon-erasure)
[![Coverage Status](https://coveralls.io/repos/github/darrenldl/reed-solomon-erasure/badge.svg?branch=master)](https://coveralls.io/github/darrenldl/reed-solomon-erasure?branch=master)
[![Crates](https://img.shields.io/crates/v/reed-solomon-erasure.svg)](https://crates.io/crates/reed-solomon-erasure)
[![Documentation](https://docs.rs/reed-solomon-erasure/badge.svg)](https://docs.rs/reed-solomon-erasure)
[![dependency status](https://deps.rs/repo/github/darrenldl/reed-solomon-erasure/status.svg)](https://deps.rs/repo/github/darrenldl/reed-solomon-erasure)

Rust implementation of Reed-Solomon erasure coding

> **This is a fork of [rust-rse/reed-solomon-erasure](https://github.com/rust-rse/reed-solomon-erasure)**,
> extended with a **BLS12-381 scalar field (`Fr`) backend** (`galois_prime`) for use in data availability
> and KZG polynomial commitment schemes. All original backends and features are preserved.

WASM builds are also available, see section **WASM usage** below for details

This is a port of [BackBlaze's Java implementation](https://github.com/Backblaze/JavaReedSolomon), [Klaus Post's Go implementation](https://github.com/klauspost/reedsolomon), and [Nicolas Trangez's Haskell implementation](https://github.com/NicolasT/reedsolomon).

Version `1.X.X` copies BackBlaze's implementation, and is less performant as there were fewer places where parallelism could be added.

Version `>= 2.0.0` copies Klaus Post's implementation. The SIMD C code is copied from Nicolas Trangez's implementation with minor modifications.

Version `>= 7.0.0` adds the BLS12-381 scalar field (`galois_prime`) backend (this fork).

See [Notes](#notes) and [License](#license) section for details.

---
 
## Fork additions (v7+)
 
### BLS12-381 scalar field backend (`galois_prime`)
 
In addition to the original `galois_8` and `galois_16` backends, this fork introduces a
`galois_prime` backend that performs Reed-Solomon erasure coding over the **scalar field of the
BLS12-381 elliptic curve** (`ark_bls12_381::fr::Fr`), via
[`ark-bls12-381`](https://github.com/arkworks-rs/curves) and
[`ark-ff`](https://github.com/arkworks-rs/algebra).
 
This enables two use cases that the original byte-oriented backends cannot address:
 
**Data Availability (DA)** — erasure-encode field elements directly, as required by DA layers
that commit to data using elliptic curve operations. Shards are native `Fr` elements, making them
directly compatible with commitment schemes without any re-encoding step.
 
**KZG polynomial commitments** — encode a polynomial's evaluations as shards, recover missing
evaluations from a subset, and feed them directly into KZG provers/verifiers without converting
between field representations.
 
#### Data packing
 
The backend exposes two ways to move between raw bytes and field elements:
 
- **`serialize` / `deserialize`** — each `Fr` element is represented as **32 bytes** (little-endian). Use this when working with field elements directly.
- **`from_data` / `into_data`** — each `Fr` element carries **31 bytes** of payload, staying safely below the modulus. Use this when encoding arbitrary byte data into field elements.
 
> **Note:** `ORDER` is currently set to `usize::MAX` as a placeholder — the true field order of
> BLS12-381 Fr is a 255-bit prime and does not fit in a `usize`. This does not affect correctness
> for normal encode/reconstruct operations but is a known limitation to be aware of.
---

## WASM usage

See [here](wasm/README.md) for details

## Rust usage
Add the following to your `Cargo.toml` for the normal version (pure Rust version)
```toml
[dependencies]
reed-solomon-erasure = "4.0"
```
or the following for the version which tries to utilise SIMD
```toml
[dependencies]
reed-solomon-erasure = { version = "4.0", features = [ "simd-accel" ] }
```
and the following to your crate root
```rust
extern crate reed_solomon_erasure;
```

NOTE: `simd-accel` is tuned for Haswell+ processors on x86-64 and not in any way for other architectures, set
environment variable `RUST_REED_SOLOMON_ERASURE_ARCH` during build to force compilation of C code for specific architecture (`-march` flag in
GCC/Clang). Even on x86-64 you can achieve better performance by setting it to `native`, but it will stop running on
older CPUs, YMMV.

## Example
```rust
#[macro_use(shards)]
extern crate reed_solomon_erasure;

use reed_solomon_erasure::galois_8::ReedSolomon;
// or use the following for Galois 2^16 backend
// use reed_solomon_erasure::galois_16::ReedSolomon;

fn main () {
    let r = ReedSolomon::new(3, 2).unwrap(); // 3 data shards, 2 parity shards

    let mut master_copy = shards!(
        [0, 1,  2,  3],
        [4, 5,  6,  7],
        [8, 9, 10, 11],
        [0, 0,  0,  0], // last 2 rows are parity shards
        [0, 0,  0,  0]
    );

    // Construct the parity shards
    r.encode(&mut master_copy).unwrap();

    // Make a copy and transform it into option shards arrangement
    // for feeding into reconstruct_shards
    let mut shards: Vec<_> = master_copy.iter().cloned().map(Some).collect();

    // We can remove up to 2 shards, which may be data or parity shards
    shards[0] = None;
    shards[4] = None;

    // Try to reconstruct missing shards
    r.reconstruct(&mut shards).unwrap();

    // Convert back to normal shard arrangement
    let result: Vec<_> = shards.into_iter().filter_map(|x| x).collect();

    assert!(r.verify(&result).unwrap());
    assert_eq!(master_copy, result);
}
```

## Backends
 
| Backend | Module | Field | Use case |
|---|---|---|---|
| GF(2⁸) | `galois_8` | Galois field 2⁸ | General erasure coding, storage |
| GF(2¹⁶) | `galois_16` | Galois field 2¹⁶ | General erasure coding, higher redundancy |
| BLS12-381 Fr | `galois_prime` | BLS12-381 scalar field (`Fr`) | Data availability, KZG commitments |

## Benchmark it yourself
You can test performance under different configurations quickly (e.g. data parity shards ratio, parallel parameters)
with standard `cargo bench` command.

## Performance
Version `1.X.X`, `2.0.0` do not utilise SIMD.

Version `2.1.0` onward uses Nicolas's C files for SIMD operations.

Machine: laptop with `Intel(R) Core(TM) i5-3337U CPU @ 1.80GHz (max 2.70GHz) 2 Cores 4 Threads`

Below shows the result of one of the test configurations, other configurations show similar results in terms of ratio.

|Configuration| Klaus Post's | >= 2.1.0 && < 4.0.0 | 2.0.X | 1.X.X |
|---|---|---|---|---|
| 10x2x1M | ~7800MB/s |~4500MB/s | ~1000MB/s | ~240MB/s |

Versions `>= 4.0.0` have not been benchmarked thoroughly yet

## Benchmarking
You can run benchmarks via `cargo bench`. To enable simd acceleration during benchmarks use `cargo bench --features simd-accel`.

## Changelog
[Changelog](CHANGELOG.md)

## Contributions
Contributions are welcome. Note that by submitting contributions, you agree to license your work under the same license used by this project as stated in the LICENSE file.

## Credits

#### BLS12-381 scalar field backend (`galois_prime`)
 
Added by [@jdetchart](https://github.com/jdetchart) — Reed-Solomon erasure coding over the BLS12-381
scalar field (`ark_bls12_381::fr::Fr`), for data availability and KZG polynomial commitment use cases.

#### Library overhaul and Galois 2^16 backend
Many thanks to the following people for overhaul of the library and introduction of Galois 2^16 backend

  - [@drskalman](https://github.com/drskalman)

  - Jeff Burdges [@burdges](https://github.com/burdges)

  - Robert Habermeier [@rphmeier](https://github.com/rphmeier)

#### WASM builds
Many thanks to Nazar Mokrynskyi [@nazar-pc](https://github.com/nazar-pc) for submitting his package for WASM builds

He is the original author of the files stored in `wasm` folder. The files may have been modified by me later.

#### AVX512 support
Many thanks to [@sakridge](https://github.com/sakridge) for adding support for AVX512 (see [PR #69](https://github.com/darrenldl/reed-solomon-erasure/pull/69))

#### build.rs improvements
Many thanks to [@ryoqun](https://github.com/ryoqun) for improving the usability of the library in the context of cross-compilation (see [PR #75](https://github.com/darrenldl/reed-solomon-erasure/pull/75))

#### no_std support
Many thanks to Nazar Mokrynskyi [@nazar-pc](https://github.com/nazar-pc) for adding `no_std` support (see [PR #90](https://github.com/darrenldl/reed-solomon-erasure/pull/90))

#### Testers
Many thanks to the following people for testing and benchmarking on various platforms

  - Laurențiu Nicola [@lnicola](https://github.com/lnicola/) (platforms: Linux, Intel)

  - Roger Andersen [@hexjelly](https://github.com/hexjelly) (platforms: Windows, AMD)

## Notes
#### Code quality review
If you'd like to evaluate the quality of this library, you may find audit comments helpful.

Simply search for "AUDIT" to see the dev notes that are aimed at facilitating code reviews.

#### Implementation notes
The `1.X.X` implementation mostly copies [BackBlaze's Java implementation](https://github.com/Backblaze/JavaReedSolomon).

`2.0.0` onward mostly copies [Klaus Post's Go implementation](https://github.com/klauspost/reedsolomon), and copies C files from [Nicolas Trangez's Haskell implementation](https://github.com/NicolasT/reedsolomon).

`7.0.0` onward adds the `galois_prime` backend — Reed-Solomon over the BLS12-381 scalar field
(`ark_bls12_381::fr::Fr`) — using the [Arkworks](https://github.com/arkworks-rs) algebra libraries.

The test suite for all versions copies [Klaus Post's Go implementation](https://github.com/klauspost/reedsolomon) as basis.

## License
#### Nicolas Trangez's Haskell Reed-Solomon implementation
The C files for SIMD operations are copied (with no/minor modifications) from [Nicolas Trangez's Haskell implementation](https://github.com/NicolasT/reedsolomon), and are under the same MIT License as used by NicolasT's project

#### TL;DR
All files are released under the MIT License
