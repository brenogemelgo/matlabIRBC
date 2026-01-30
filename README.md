# cudaLBM boundary-condition generators

Target solver:
https://github.com/nd-92/cudaLBM

## Structure

- `densityBased/`   — generators for density-based formulations
- `pressureBased/`  — generators for pressure-based formulations
- `commonFunctions/` — shared symbolic and helper routines

## What you edit

For each formulation (`densityBased` or `pressureBased`), only these files
are intended to be modified:

- `neumannBoundaries.m`
- `staticBoundaries.m`

All other files should be treated as common infrastructure.

## Configuration

Inside both generator files:

### Phase-field coupling
```matlab
phase_field = true;   % or false
