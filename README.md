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

All other files are shared infrastructure.

## Configuration

Both `neumannBoundaries.m` and `staticBoundaries.m` expose the same controls.

### Phase-field coupling
```matlab
phase_field = true; % or false
```

### Boundary selection (wantNames)
`wantNames` defines which boundaries are generated.
Only the names listed in this array are emitted.

Example (`ssmd`):
```matlab
wantNames = [ ...
    "WEST_NORTH_FRONT", "EAST_NORTH_FRONT", ...
    ...
    "WEST_NORTH", "EAST_NORTH", ...
    "WEST_FRONT", "EAST_FRONT", ...
    "NORTH_FRONT", ...
    ...
    "NORTH", "FRONT"
];
```
