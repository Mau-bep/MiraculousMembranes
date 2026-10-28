# MiraculousMembranes

Simulations of triangulated membranes (bending, tension and volume energies,
beads, adaptive remeshing) built on a modified copy of
[geometry-central](https://geometry-central.net). The simulation program is
`main_cluster`; it reads a JSON config file. Everything is in
`projects/geometric-flow`.

All third-party code the build needs is in `deps/` (see below). The only things
to install by hand are a compiler toolchain and a few system headers for the
viewer.

## Requirements

Tested on Ubuntu 22.04 with g++ 11.4 and CMake 3.22.

| What | Version | Used for |
|---|---|---|
| C++ compiler with OpenMP (g++ or clang) | C++11 | everything |
| CMake | >= 3.10 | build |
| make, git | | build, clone |
| X11 and OpenGL development headers | | polyscope (the viewer) |

The viewer is always configured and built with the rest, so the X11/OpenGL
headers are needed on every machine, also one without a display. Nothing
opens a window during the build.

Ubuntu / Debian:

```bash
sudo apt install build-essential cmake git \
    libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libgl-dev
```

On other distributions install the equivalent X11 (Xrandr, Xinerama, Xcursor,
Xi) and OpenGL development packages. If one is missing, `cmake` stops with
"The RandR / Xinerama / Xcursor headers were not found", or `make` stops with
`X11/extensions/XInput2.h: No such file or directory` or
`GL/gl.h: No such file or directory`.

You do not need to install Eigen, OpenBLAS, SuiteSparse, Ipopt or Python to
build.

## Build

```bash
git clone git@github.com:Mau-bep/MiraculousMembranes.git
cd MiraculousMembranes/projects/geometric-flow
mkdir build && cd build
cmake ..
make -j$(nproc) main_cluster   # the simulation only: about 1.5 minutes on 8 cores
make -j$(nproc)                # everything (viewer, post-processing, tests): a few more minutes
```

The executables are in `build/bin`:

| Executable | Purpose |
|---|---|
| `main_cluster` | runs a simulation: `main_cluster <config.json>`; `--check` only validates the config |
| `main_visualize` | polyscope viewer |
| `PostProcessing`, `Processing_tube` | post-processing tools |
| `Test_*` | finite-difference, exactness and other checks |

CMake options:

| Option | Default | Meaning |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `Release` | `-O3 -march=native -ffp-contract=off` |
| `MM_USE_SYSTEM_EIGEN` | `OFF` | `ON` uses the system Eigen3 (>= 3.4) instead of `deps/eigen-3.4.0` |
| `SUITESPARSE` | `OFF` | geometry-central's optional SuiteSparse solvers (the project does not call them) |

`-march=native` makes the binary specific to the CPU type it was built on.
On a cluster, build on a node of the type that will run the jobs.
`-ffp-contract=off` keeps results reproducible between refactors; results are
bit-identical per machine, not across CPU types.

## Check the install

```bash
cd projects/geometric-flow
regression/check_install.sh        # uses ./build
```

It validates every config in `regression/` with `main_cluster --check` and runs
`Test_exactness` on them. It ends with `Install check passed.` The configs in
`regression/` name the checkout of whoever wrote them (`init_file`), and the
script runs on copies that point at yours.

To compare two builds on the same machine (before and after a change), keep
the dump of each and compare them:

```bash
regression/check_install.sh build_old /tmp/old.txt
regression/check_install.sh build /tmp/new.txt
cmp /tmp/old.txt /tmp/new.txt && echo identical
```

## Python scripts (optional)

The launchers and analysis scripts in `projects/geometric-flow/Scripts` are
plain Python 3:

- config and job generation: `numpy`, `jinja2`
- analysis and plots: `matplotlib`, `scipy`
- rendering (a few scripts only): `ovito`, `opencv-python`

`environment.yml` lists the first two groups for conda
(`conda env create -f environment.yml`). Its `prefix:` line is the path of the
environment on the cluster; edit or remove it on another machine.

## What is in `deps/`

| Folder | What | Notes |
|---|---|---|
| `geometry-central` | mesh data structures and remeshing | modified copy: triangle area, edge length, dihedral angle and edge-regularity terms with their gradients and Hessians are declared in `include/geometrycentral/surface/vertex_position_geometry.h` and implemented in `core/src/geometry.cpp`, which geometry-central compiles |
| `polyscope` | viewer, with its own glfw, glad, imgui, glm | needs the X11/OpenGL headers above |
| `eigen-3.4.0` | Eigen headers | copied from Ubuntu's `libeigen3-dev` 3.4.0-2ubuntu2. The code needs Eigen >= 3.4 (it uses `Eigen::Vector<T, N>`), so Eigen 3.3.x does not compile |
| `EigenRand-0.5.1` | random numbers for Eigen | header only |
| `nlohmann_json` | JSON parser for the config files | header only |
| `googletest` | unit tests | configured but not used by any target |

Not bundled any more: OpenBLAS, alglib, jsoncpp, taucs, Ipopt and HSL. No
active target uses them. They were only needed by the disabled
`NormalNLP`/`ShapeNLP` sources (`include/NormalNLP.hpp`, `include/ShapeNLP.hpp`),
which need Ipopt with HSL (HSL is licensed by STFC and cannot be redistributed).

## Layout

| Path | Contents |
|---|---|
| `projects/geometric-flow/src`, `include` | simulation code |
| `projects/geometric-flow/regression` | reference configs and `check_install.sh` |
| `projects/geometric-flow/Scripts` | launchers and analysis |
| `projects/geometric-flow/Templates` | job-script templates (SLURM; adapt modules and paths to your cluster) |
| `projects/geometric-flow/Meshes`, `input` | starting meshes |
| `projects/geometric-flow/Studies` | study-specific scripts and notes |
| `deps` | third-party code |
