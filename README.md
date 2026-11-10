# calc_tidalQ

`calc_tidalQ` calculates the complex degree-2 tidal Love numbers `h2` and
`k2`, their phase lags, and two corresponding tidal quality factors for
the solid Earth. It solves a radial boundary-value problem for a
spherically symmetric mantle with a uniform-density fluid core. The
input shear modulus is complex and should describe the mantle at the
requested forcing frequency. The model uses fixed Earth radii and
Earth–Moon constants.

This program is based on code developed for Korenaga (2025).
[*Tidal dissipation within Earth's solidifying magma ocean: I. Effects of inertia and lunar orbital eccentricity*](https://doi.org/10.1016/j.icarus.2025.116756),
*Icarus*, 442, 116756.

## Citation

Korenaga, J. (2026). *calc_tidalQ* (Version v1.0.0) [Computer software].
Zenodo. [10.5281/zenodo.24095696](https://doi.org/10.5281/zenodo.24095696)

```bibtex
@misc{korenaga2026calc_tidalq,
  author       = {Korenaga, Jun},
  title        = {{calc\_tidalQ}},
  year         = {2026},
  howpublished = {Zenodo},
  note         = {Version v1.0.0},
  doi          = {10.5281/zenodo.24095696},
  url          = {https://doi.org/10.5281/zenodo.24095696}
}
```

## Build and install

You need a C++17 compiler, GNU Make, `pkg-config`, and LAPACKE, the C
interface to LAPACK. Install the development dependency before building.

- Debian or Ubuntu — `sudo apt-get install liblapacke-dev pkg-config`
- macOS with Homebrew — `brew install openblas pkgconf`. OpenBLAS provides
  LAPACKE.

The Makefile looks for `lapacke` through `pkg-config` on Linux and
`openblas` on macOS. For a custom installation, set `LAPACK_CPPFLAGS`
and `LAPACK_LIBS` to the appropriate include and link flags.

```sh
make
make install PREFIX="$HOME/.local"
```

The default installation directory is `/usr/local/bin`. Set `PREFIX` for
a different installation prefix, or set `BINDIR` directly. Package
builders can stage the installation with `DESTDIR`. To select a compiler,
run `make clean` and rebuild with, for example, `make CXX=clang++`.

## Run the M₂ example

```sh
./calc_tidalQ -W1.4056e-4 -M2 -S examples/model_tidal_Q_model1_1.dat
```

The first row of `examples/tidalQ.1.expected.out` gives the expected
result for this input.

| Option | Meaning |
| --- | --- |
| `-W` | Tide-raising angular frequency in rad/s. Zero requires `-I`. |
| `-M` | Azimuthal order, 0, 1, or 2. The spherical harmonic degree is fixed at 2. |
| `-S` | Path to the mantle model. |
| `-I` | Disable inertia in the mantle equations and at the core–mantle boundary. |
| `-R` | Target radial spacing in meters. The default is 1000. The solver rounds the number of intervals and fits the grid to the two endpoints. |
| `-D` | Write a radial solution to `<root>.01.dat`. An existing file at that path is overwritten. |
| `-V` | Print the number of grid points and the actual spacing to standard error. |
| `-h`, `--help` | Show command-line help. |

Option values can be attached to the flag or passed separately. Quote
paths that contain spaces.

## Run the M₂ and nodal batch example

After building, run the paired sample inputs and compare the two rows:

```sh
./examples/calc_tidalQ.sh examples sample-results
diff -u examples/tidalQ.1.expected.out sample-results/tidalQ.1.out
```

The script writes M₂ first and the nodal tide second. Small floating-point
differences can occur across builds, especially in column 11, the scaled
equation residual.

With no arguments, the script uses the sample pair beside it and writes
to `sample-results`. For other inputs, pass `MODEL_DIR OUTPUT_DIR`. Name
each pair `model_tidal_Q_model<ID>_1.dat` for M₂ (`1.4056e-4` rad/s) and
`model_tidal_Q_model<ID>_2.dat` for the nodal tide (`1.065e-8` rad/s).
The script writes `tidalQ.<ID>.out` for each pair and replaces existing
files with those names. Set `TIDALQ` to use another executable.

The sample inputs were generated at periods of 44,700 s for M₂ and
5.9 × 10⁸ s for the nodal tide.

## Mantle model format

Each row has five whitespace-separated columns in SI units.

| Column | Quantity | Unit |
| --- | --- | --- |
| 1 | Radius | m |
| 2 | Density | kg/m³ |
| 3 | Bulk modulus | Pa |
| 4 | Real part of shear modulus | Pa |
| 5 | Imaginary part of shear modulus | Pa |

Start at the surface (6371000 m) and end at the core–mantle boundary
(3480000 m), in descending radius. Use positive density and bulk
modulus and a nonzero complex shear modulus. Blank lines and comments
beginning with `#` are allowed. To represent an internal discontinuity,
give two rows at the same radius with the shallower value first and the
deeper value second. The code interpolates between rows onto its
radial grid. Prepare the shear modulus for each frequency you want to
evaluate.

## Output

A successful calculation writes one whitespace-separated row to
standard output, in this order.

| Columns | Values |
| --- | --- |
| 1–2 | Angular frequency and azimuthal order |
| 3–4 | `Q_h` and `Q_k` |
| 5–6 | Real and imaginary parts of `h2` |
| 7–8 | Real and imaginary parts of `k2` |
| 9–10 | Phase lags of `h2` and `k2` in radians |
| 11–12 | Scaled residual of the discrete system and solve count (1) |

For a complex Love number `L`, the code reports
`atan2(-Im(L), Re(L))` as the lag and `abs(Re(L)) / abs(Im(L))` as
`Q`. For a negative real response with zero imaginary part, either
sign of π represents the same phase. If `abs(Im(L)) <= 1e-10`, it reports `1e10` as a fallback.
The radial equations and boundary conditions form a complex banded
system `A x = b`. LAPACK's `ZGBSV` solves it in one call. Column 11
reports `||A x - b||∞ / (||A||∞ ||x||∞ + ||b||∞)`, using the maximum
absolute row sum for the matrix norm. This measures how closely the
solution satisfies the discrete equations. Use a smaller `-R` to
assess radial discretization error.

The optional radial dump starts at the core–mantle boundary and ends at
the surface. Each row contains radius (m), the real and imaginary parts
of the six solution variables in solver order, and the real and
imaginary parts of shear modulus (Pa). Variables 1–2 are displacement
components (m), 3–4 are traction components (Pa), 5 is the potential
variable (m²/s²), and 6 is the associated gravity variable (m/s²).

## Reproduce a calculation

Record the exact code revision, compiler and build flags, LAPACK
implementation and version, mantle file and its provenance, and complete
command line.

The implementation follows the formulation cited in the source by
[Sabadini, Vermeersen, and Cambiotti (2016)](https://link.springer.com/book/10.1007/978-94-017-7552-6).
The core boundary calculation uses a P-wave speed of 8060 m/s just below
the core–mantle boundary, identified in the source as a value from
[Dziewonski and Anderson's PREM (1981)](https://doi.org/10.1016/0031-9201(81)90046-7).

## License

calc_tidalQ is released under the [MIT License](LICENSE). The license file
retains the original copyright notice for code shared with
[tomo2d](https://github.com/jun-korenaga/tomo2d).

The banded solve uses [LAPACK](https://github.com/Reference-LAPACK/lapack),
which has a modified BSD license.
