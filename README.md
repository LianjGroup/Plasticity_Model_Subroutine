# enHill48 VUMAT

A Fortran material subroutine for Abaqus/Explicit implementing the evolving non-associated Hill48 plasticity model described in [1]. The yield function is fitted to directional flow stresses; a separate plastic potential is fitted to r-values. Both evolve with equivalent plastic strain.

**Any researcher using this subroutine must cite the original paper [1]: [doi:10.1016/j.ijsolstr.2017.04.007](https://doi.org/10.1016/j.ijsolstr.2017.04.007).**

This version includes isotropic linear elasticity and anisotropic plasticity. It does not include damage or fracture.

## Files

- `enHill48-VUMAT-v05.f`: stress update (`vumat`) and table interpolation (`ahard`).
- `input_template.inp`: QP1200 material definition. Include it in a complete Abaqus/Explicit model.

## Material input

The first eight user-material constants are:

```text
E, nu, 0, 0, 0, 0, 0, 0
```

Each following row contains:

```text
plastic strain, s0, s45, s90, sEB, r0, r45, r90
```

`s0`, `s45`, and `s90` are uniaxial flow stresses at 0°, 45°, and 90° to the rolling direction. `sEB` is the equibiaxial flow stress. The three r-values refer to the same uniaxial directions.

For `n` data rows, set `*User Material, constants=8+8n`. The supplied template contains 400 rows (`constants=3208`), with E=190000, nu=0.3, and density=7.85e-09. These values are compatible with mm–N–s–tonne units (stress in MPa).

Use at least two rows with strictly increasing plastic strain, starting at zero. The routine uses linear interpolation within the table and linear extrapolation outside it. Keep the data range wide enough to cover the simulation.

## Setup and output

Use three-dimensional solid elements (`ndir=3`, `nshr=3`). Align material axes 1, 2, and 3 with the rolling, transverse, and thickness directions. Shell and plane-stress interfaces are not supported by this implementation.

Allocate two state variables:

```abaqus
*Depvar
2,
```

| Output | Meaning |
| --- | --- |
| SDV1 | Accumulated equivalent plastic strain |
| SDV2 | Hill48 equivalent stress |

Request `S, SDV` in element field output. SDV1 is stored by the subroutine rather than in the built-in PEEQ variable; SDV2 differs from Mises stress.

## Run

With a supported Fortran compiler configured, submit a complete model from the Abaqus command environment:

```text
abaqus job=example input=example.inp user=enHill48-VUMAT-v05.f interactive
```

The source was written for Abaqus 6.11-1. Check the VUMAT interface against your installed version.

## Notes

- Stress integration uses an elastic predictor and Tangent Cutting Plane plastic correction. The iteration limit is normally 11 corrections, with no warning if convergence is not reached.
- Internal and inelastic energy outputs (`enerInternNew`, `enerInelasNew`) are not updated.
- Temperature and strain-rate effects are not included.
- Check single-element responses in the 0°, 45°, 90°, and equibiaxial loading directions before using the model. This README was prepared from the source and input template; no Abaqus job was run.

The source header credits Junhe Lian (February 2016) and requires the author's permission for distribution. Retain that attribution and restriction.

## Reference

[1] Lian, J., Shen, F., Jia, X., Ahn, D. C., Chae, D. C., Münstermann, S., and Bleck, W. (2018). An evolving non-associated Hill48 plasticity model accounting for anisotropic hardening and r-value evolution and its application to forming limit prediction. *International Journal of Solids and Structures*, 151, 20–44. [doi:10.1016/j.ijsolstr.2017.04.007](https://doi.org/10.1016/j.ijsolstr.2017.04.007)

```bibtex
@article{Lian2018enHill48,
  author  = {Lian, Junhe and Shen, Fuhui and Jia, Xiaoxu and
             Ahn, Deok Chan and Chae, Dong Chul and
             M{\"u}nstermann, Sebastian and Bleck, Wolfgang},
  title   = {An evolving non-associated Hill48 plasticity model accounting for anisotropic hardening and r-value evolution and its application to forming limit prediction},
  journal = {International Journal of Solids and Structures},
  year    = {2018},
  volume  = {151},
  pages   = {20--44},
  doi     = {10.1016/j.ijsolstr.2017.04.007}
}
```

