![MSRV](https://img.shields.io/badge/MSRV-1.58.1-green?style=for-the-badge)

# Table of Contents
[[_TOC_]]

## Pre-Installation Steps
### Installing Rust
Installing the Rust compiler and `cargo` build system can be done by following instructions [here](https://doc.rust-lang.org/book/ch01-01-installation.html).
No root privilege is needed.
For Unix-like systems the following command is sufficient. 
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
This will install `rustup`, which is similar to `nvm` for NodeJS, which manages Rust toolchain versions and updates for you.
Make sure that you are on our version of Rust by running
```bash
rustup default 1.58.1
```
If your shell fails to resolve `rustc` or `rustup`, you should add `~/.cargo/bin` to your `PATH` variable.

If Rust is already installed run `rustup default 1.58.1` to make sure you're using our supported version of the Rust toolchain.

**IMPORTANT:** If you're on a M1 Mac be sure to set the `x86_64` target instead of the `arm64` target.

### Installing Other Dependencies
If not present you will also need to install:
- [wget](https://www.gnu.org/software/wget/)
- [unzip](http://infozip.sourceforge.net/UnZip.html)
- something equivalent to [build-essential](https://packages.ubuntu.com/impish/build-essential)
- [libclang](https://clang.llvm.org/doxygen/group__CINDEX.html)
- [boost](https://www.boost.org/)
- [autoconf](https://www.gnu.org/software/autoconf/)

Most if not all of these will already be installed on your system.

### Running Tests
This is a monorepo with multiple small, interdependent libraries. To make sure everything is working as expected run
```bash
cargo test
```
from the root of this repository before proceeding.
If there are any test failures please open a Gitlab Issue with the ~Bug label including:

## Building Binaries

### Local Environment
If you known what you're doing and followed previous instructions you can build all binaries with
```bash
cargo build --release
```
All binaries will be in `target/release/` relative to the root of the repo.
If you're even a little bit unsure do not attempt to build this project.
Instead use a pre-built image. These images contain all binaries in this repo along configuration data.

## Documentation
Documentation can be build with `cargo doc`
