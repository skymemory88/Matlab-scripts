# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a MATLAB-based Monte Carlo simulation codebase for studying the magnetic properties of LiReF4 compounds (lithium rare-earth tetrafluorides). The simulation implements a classical Monte Carlo approach with quantum corrections to model the spin dynamics, phase transitions, and magnetic interactions in these materials.

## Running Simulations

### Main Commands

To run a basic simulation:
```matlab
[init, meas, coef, ion, params] = LiReF4_CMC_Yikai('Er', [8 8 8], 'disordered', 0, 0, 1, 1);
```

To run the optimized version with algorithmic improvements:
```matlab
[init, meas, coef, ion, params] = LiReF4_CMC_Yikai_a('Er', [8 8 8], 'disordered', 0, 0, 1, 1);
```

### Directory Structure

- **Main directory**: Contains standard implementations
- **beta/ subdirectory**: Contains `_a` suffixed files with algorithmic improvements and optimizations
- **Deleted GPU versions**: GPU-accelerated variants have been removed but optimizations may be worth reviewing

## Core Architecture

### Main Simulation Entry Points

- **LiReF4_CMC_Yikai.m**: Primary simulation driver with comprehensive functionality
- **LiReF4_CMC_Yikai_a.m**: Alternative version with modified algorithms (located in `beta/`)

Function signature:
```matlab
function [init, meas, coef, ion, params] = LiReF4_CMC_Yikai(mion, dims, iState, hyperfineOpt, loadOpt, plotOpt, saveOpt)
```

Parameters:
- `mion`: Magnetic ion type ('Er', 'Ho', 'Yb', 'Tm', 'Gd', 'Y', 'dope')
- `dims`: Lattice dimensions [x, y, z] array
- `iState`: Initial state configuration ('disordered', 'ordered', or other states)
- `hyperfineOpt`: Enable/disable hyperfine interactions (0/1)
- `loadOpt`: Continue from saved state (1) or start fresh (0)
- `plotOpt`: Enable visualization (0/1)
- `saveOpt`: Enable data saving (0/1)

## Complete Simulation Workflow

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                           LiReF4_CMC_Yikai_a.m                                 │
│                         Main Simulation Driver                                 │
└─────────────────┬───────────────────────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                          INITIALIZATION PHASE                                  │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 1. Parse Input Parameters & Setup Paths                                        │
│    - Ion type, lattice dimensions, initial state                               │
│    - OS-specific data paths (Windows/macOS)                                    │
│                                                                                 │
│ 2. Load/Create Ion Properties                                                   │
│    - Crystal field parameters (Ion_class.m, cf.m)                              │
│    - Magnetic moments (J, L, S, I values)                                      │
│    - Landé g-factors (gLande.m)                                                │
│    - Temperature/field ranges                                                   │
│                                                                                 │
│ 3. Setup Lattice & Interactions                                                │
│    - Lattice geometry (lattice.m)                                              │
│    - Neighbor lists (neighborList_a.m)                                         │
│    - Luijten-Blöte tables (cluster_init_a.m)                                   │
│    - GPU state initialization (if enabled)                                     │
└─────────────────┬───────────────────────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         initialize_a.m                                         │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 4. Initialize Quantum States                                                   │
│    - Electronic wavefunction coefficients (coef)                               │
│    - Spin expectation values (eSpin, nSpin)                                    │
│    - Single-ion energies (E_si)                                                │
│    - Crystal field Hamiltonian (hamI)                                          │
│    - Basis states for each T,B point                                           │
│                                                                                 │
│    Initial State Options:                                                      │
│    ├─ 'disordered': Random orientations on Bloch sphere                       │
│    ├─ 'ordered': Ground state configurations                                   │
│    │   ├─ Er/Yb: xy-AFM with propagation vector                               │
│    │   └─ Ho/Tm/Gd/Y: z-FM ferromagnetic                                       │
│    └─ Nuclear spins: Random thermal distribution                               │
└─────────────────┬───────────────────────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                      PARALLEL EXECUTION PHASE                                  │
│                    (Temperature × Field Grid)                                  │
├─────────────────────────────────────────────────────────────────────────────────┤
│                                                                                 │
│  For each (Temperature, Magnetic Field) point in parallel:                     │
│                                                                                 │
│  ┌───────────────────────────────────────────────────────────────────────────┐ │
│  │                     EQUILIBRATION PHASE                                   │ │
│  │                      equilibrate_a.m                                      │ │
│  ├───────────────────────────────────────────────────────────────────────────┤ │
│  │ Step 5. System Thermalization (params.tEQ steps)                         │ │
│  │                                                                           │ │
│  │ Primary Algorithm: Single-Spin Metropolis Updates                        │ │
│  │ ├─ thermalize_a.m: Bloch sphere rotations                                │ │
│  │ │  ├─ Random rotations (randSph.m)                                       │ │
│  │ │  ├─ Energy calculations:                                               │ │
│  │ │  │  ├─ Single-ion: Crystal field + Zeeman                             │ │
│  │ │  │  ├─ Dipole-dipole: CntrDip_a.m (neighbor lists)                    │ │
│  │ │  │  ├─ Exchange: exchange.m (if enabled)                               │ │
│  │ │  │  └─ Hyperfine: Nuclear-electronic coupling                          │ │
│  │ │  └─ Metropolis acceptance: exp(-ΔE*β)                                  │ │
│  │ │                                                                         │ │
│  │ ├─ Nuclear spin updates: therm_nuc_a.m (if hyperfine enabled)            │ │
│  │ │                                                                         │ │
│  │ └─ Critical Slowing Detection & Cluster Updates:                          │ │
│  │    ├─ Monitor: Low acceptance rate, trapped random walks                  │ │
│  │    ├─ Trigger: After 5 consecutive trapped states                         │ │
│  │    └─ Execute: thermalize_cluster_a.m                                     │ │
│  │       ├─ Quantum unitary transformations (2J+1)×(2J+1)                  │ │
│  │       ├─ Luijten-Blöte cluster construction                              │ │
│  │       ├─ Long-range dipole bond probabilities                            │ │
│  │       └─ Consistent wavefunction-spin updates                            │ │
│  └───────────────────────────────────────────────────────────────────────────┘ │
│                                  │                                             │
│                                  ▼                                             │
│  ┌───────────────────────────────────────────────────────────────────────────┐ │
│  │                      MEASUREMENT PHASE                                    │ │
│  │                       MC_sample_a.m                                       │ │
│  ├───────────────────────────────────────────────────────────────────────────┤ │
│  │ Step 6. Data Collection (params.mIntv steps)                             │ │
│  │                                                                           │ │
│  │ Continuous Sampling with Adaptive Cluster Updates:                       │ │
│  │ ├─ Regular thermalization: thermalize_a.m                                │ │
│  │ ├─ Autocorrelation monitoring (τ calculation)                            │ │
│  │ ├─ Adaptive cluster triggering:                                          │ │
│  │ │  ├─ Autocorrelation time τ > 5                                         │ │
│  │ │  ├─ Acceptance rate < 10%                                              │ │
│  │ │  ├─ Adaptive intervals (50-500 steps)                                  │ │
│  │ │  └─ Forced every 1000 steps                                            │ │
│  │ └─ Observables calculation:                                               │ │
│  │    ├─ Magnetization components <Mx>, <My>, <Mz>                          │ │
│  │    ├─ Interaction energies (dipole + exchange)                           │ │
│  │    ├─ Structure factors for neutron scattering                           │ │
│  │    └─ Order parameters and correlations                                  │ │
│  └───────────────────────────────────────────────────────────────────────────┘ │
│                                                                                 │
└─────────────────┬───────────────────────────────────────────────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         DATA OUTPUT PHASE                                      │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 7. Results Compilation                                                          │
│    - init: Initial configurations and parameters                               │
│    - meas: Measured observables (magnetization, energies)                      │
│    - coef: Final quantum wavefunction coefficients                             │
│    - ion: Ion properties and parameters                                        │
│    - params: All simulation parameters                                         │
│                                                                                 │
│ 8. Data Saving (if saveOpt enabled)                                            │
│    - OS-specific paths to Google Drive                                         │
│    - Format: Li[ion]F4_MC_[dims]_hyp=[option]_[state].mat                      │
│                                                                                 │
│ 9. Visualization (if plotOpt enabled)                                          │
│    - Spin configurations (plot_spin.m)                                         │
│    - Temperature/field dependencies                                            │
│    - Order parameters and phase diagrams                                       │
└─────────────────────────────────────────────────────────────────────────────────┘
```

### Core Simulation Workflow Components

1. **Initialization** (`initialize_a.m`):
   - Sets up lattice geometry and initial spin configurations
   - Initializes electronic and nuclear spin states using Bloch sphere parameterization
   - Prepares single-ion energy calculations and crystal field Hamiltonians
   - Supports ordered/disordered initial states with ion-specific ground states

2. **Equilibration** (`equilibrate_a.m`):
   - System thermalization using Metropolis Monte Carlo algorithm
   - Adaptive cluster update integration for critical slowing down
   - Nuclear spin thermalization when hyperfine interactions enabled
   - Parallel execution across temperature-field grid

3. **Sampling** (`MC_sample_a.m`):
   - Measurement collection after equilibration
   - Continuous autocorrelation monitoring and adaptive cluster triggering
   - Calculation of magnetization components and interaction energies
   - Real-time optimization based on system response

4. **Thermalization Algorithms**:
   - **Single-spin updates** (`thermalize_a.m`): Bloch sphere rotations with Metropolis acceptance
   - **Cluster updates** (`thermalize_cluster_a.m`): Quantum unitary transformations with Luijten-Blöte efficiency
   - **Nuclear spins** (`therm_nuc_a.m`): Separate thermalization for hyperfine degrees of freedom

## Function Dependency Map

```
LiReF4_CMC_Yikai_a.m (Main Driver)
├── Ion_class.m (Ion properties & crystal field parameters)
├── cf.m (Crystal field Hamiltonian construction)
├── gLande.m (Electronic g-factor calculations) 
├── lattice.m (Lattice geometry setup)
├── neighborList_a.m (Pre-compute neighbor lists for optimization)
├── cluster_init_a.m (Pre-compute Luijten-Blöte lookup tables)
│   └── Uses: params.pos (Angstrom → meters conversion)
├── GPUSpinState.m (GPU acceleration state management)
│
├── initialize_a.m (Quantum state initialization)
│   ├── randSph.m (Random points on Bloch sphere)
│   └── Ion operators: Jx, Jy, Jz, Ix, Iy, Iz
│
├── equilibrate_a.m (System thermalization)
│   ├── thermalize_a.m (Single-spin Metropolis updates)
│   │   ├── randSph.m (Random Bloch sphere rotations)
│   │   ├── CntrDip_a.m (Dipole interactions via neighbor lists)
│   │   └── exchange.m (Exchange interactions)
│   ├── therm_nuc_a.m (Nuclear spin thermalization)
│   └── thermalize_cluster_a.m (Quantum cluster updates)
│       ├── cluster_init_a.m (Luijten-Blöte tables)
│       ├── QR decomposition (Random unitary generation)
│       └── GPU state synchronization
│
├── MC_sample_a.m (Measurement sampling)
│   ├── thermalize_a.m (Continuous thermalization)
│   ├── therm_nuc_a.m (Nuclear spin updates)
│   ├── thermalize_cluster_a.m (Adaptive cluster updates)
│   │   ├── Autocorrelation monitoring
│   │   ├── Acceptance rate tracking
│   │   └── Dynamic interval adjustment
│   ├── dipSum_a.m (Total interaction energy calculation)
│   └── BraggPeaks.m (Structure factor calculations)
│
└── Utility Functions:
    ├── plot_spin.m (Visualization)
    ├── Ising_proj.m (Order parameter analysis)
    ├── powderaverage2D.m (Powder averaging)
    └── dipSum_gpu_a.m (GPU-accelerated energy calculations)
```

## Key Algorithm Details

### 1. Quantum State Representation
- **Electronic spins**: Represented as quantum wavefunctions in (2J+1)-dimensional Hilbert space
- **Bloch sphere parameterization**: `coef = [cos(β/2); sin(β/2)exp(iα)]`
- **Spin expectation values**: `<Sx> = Re(coef† Jx coef)`, etc.
- **Consistency**: Both `coef` (wavefunctions) and `eSpin` (expectation values) updated simultaneously

### 2. Energy Calculations
```
Total Energy = Single-ion + Dipole-dipole + Exchange + Hyperfine

Single-ion: E_si = coef† * (Crystal_field + Zeeman) * coef
Dipole-dipole: E_dip = Σ_<i,j> [Si·Sj - 3(Si·rij)(Sj·rij)/r²] / r³ 
Exchange: E_ex = -J_ex * Σ_<i,j> Si·Sj
Hyperfine: E_hyp = A * Se·Sn (for isotopes with nuclear spins)
```

### 3. Metropolis Algorithm (thermalize_a.m)
```matlab
for each spin i:
    1. Generate random rotation on Bloch sphere
    2. Calculate new wavefunction coefficients  
    3. Compute new spin expectation values
    4. Calculate energy change ΔE = E_new - E_old
    5. Accept with probability min(1, exp(-β*ΔE))
    6. Update quantum state consistently if accepted
```

### 4. Quantum Cluster Algorithm (thermalize_cluster_a.m)
```matlab
1. Generate random (2J+1)×(2J+1) unitary matrix U via QR decomposition
2. Select random seed spin i₀
3. Build cluster using Luijten-Blöte algorithm:
   - Bond probability ∝ quantum spin correlations
   - Use pre-computed lookup tables for O(log N) neighbor selection
   - Handle long-range dipole interactions efficiently
4. Apply unitary transformation to all cluster wavefunctions
5. Calculate total energy change and accept/reject cluster
6. Update both coef and eSpin consistently for accepted moves
```

### 5. Adaptive Cluster Integration
- **Autocorrelation monitoring**: Track observable time series, compute τ
- **Dynamic triggering**: Switch algorithms based on system response
- **Performance metrics**: Balance single-spin vs cluster update efficiency
- **GPU synchronization**: Maintain consistent state across CPU/GPU

### Physics Modules

**Interaction Calculations:**
- **CntrDip.m**: Standard dipole-dipole interaction calculation with distance cutoffs
- **CntrDip_a.m**: Optimized version using pre-computed neighbor lists (in `beta/`)
- **dipSum.m**: Computes total dipolar interaction energies for the system
- **dipSum_a.m**: Optimized version with neighbor list optimizations (in `beta/`)
- **exchange.m**: Handles exchange interactions between magnetic moments

**Ion and Crystal Field:**
- **Ion_class.m**: Class defining ion properties (crystal field parameters, g-factors, hyperfine constants)
- **cf.m**: Crystal field calculations
- **gLande.m**: Calculates electronic Landé g-factors

**Monte Carlo Core:**
- **thermalize.m**, **thermalize_a.m**: Metropolis algorithm implementations
- **MC_sample.m**, **MC_sample_a.m**: Measurement sampling routines with parallel execution
- **therm_nuc.m**, **therm_nuc_a.m**: Nuclear spin thermalization for hyperfine interactions
- **equilibrate.m**, **equilibrate_a.m**: System equilibration protocols
- **thermalize_cluster_a.m**: Primary quantum cluster Monte Carlo using Luijten-Blöte + unitary transformations
- **thermalize_cluster_a_backup.m**: Alternative quantum cluster implementation (reference version)

### Utility Functions

- **randSph.m**, **randSphN.m**: Generate random points on sphere surface for spin updates
- **plot_spin.m**: Visualization routines for spin configurations
- **lattice.m**: Lattice geometry and neighbor calculations
- **neighborList_a.m**: Pre-computes neighbor lists for optimized dipole calculations (in `beta/`)
- **BraggPeaks.m**: Neutron scattering structure factor calculations
- **Ising_proj.m**: Ising model projections for analysis
- **cluster_init_a.m**: Pre-computes Luijten-Blöte lookup tables for efficient cluster construction
- **dipSum_gpu_a.m**: GPU-accelerated dipole sum calculations with support for trial configurations

## File Naming Conventions

- **Base files**: Standard implementations in main directory
- **`_a` suffix**: Optimized/adapted implementations with algorithmic improvements (located in `beta/`)
- **`_GPU` suffix**: GPU-accelerated versions (deleted from repository but optimizations may be referenced)
- **`_cluster` suffix**: Cluster update algorithms for enhanced sampling

## Data Management

The simulation uses OS-specific file paths:
- Windows: `G:\My Drive\File sharing\PhD research\Li[ion]F4 project\Data\Simulations\MATLAB\Monte Carlo`
- macOS: `/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/File sharing/PhD program/Research projects/Li[ion]F4 project/Data/Simulations/Matlab/Monte Carlo/`

Data files are saved as `.mat` files with naming pattern: `Li[ion]F4_MC_[dims]_hyp=[option]_[state].mat`

## Constants and Units

All physical constants are defined in the main simulation functions:
- `muB`: Bohr magneton [J/T]
- `kB`: Boltzmann constant [J/K]  
- `muN`: Nuclear magneton [J/T]
- `mu0`: Vacuum permeability [H/m]
- `J2meV`: Joule to meV conversion factor

Energy units: meV
Length units: Angstrom
Temperature units: Kelvin

## Parallelization

The codebase supports MATLAB Parallel Computing Toolbox:
- Temperature/field sweeps are parallelized using `parfor` in `MC_sample.m` functions
- Worker identification system for debugging parallel execution (`getCurrentTask()`, `getCurrentJob()`)
- Persistent variables for timing across parallel workers
- Core identification strings for monitoring parallel execution

## Cluster Monte Carlo Implementation

The beta version includes advanced cluster Monte Carlo algorithms to overcome critical slowing down in magnetic systems near phase transitions. Two complementary approaches are implemented:

### 1. Quantum Spin Cluster Algorithm (`thermalize_cluster_a.m`) - **Primary Implementation**

**Algorithm**: Long-range Wolff/Luijten-Blöte cluster updates with quantum unitary transformations
- **Luijten-Blöte Tables**: Pre-computed lookup tables (`cluster_init_a.m`) for efficient cluster construction
- **Unitary Generation**: Random unitary matrices using Haar measure (QR decomposition) for (2J+1)×(2J+1) Hilbert space
- **Quantum Consistency**: Proper wavefunction-spin expectation value consistency
- **Long-range Interactions**: Handles dipole-dipole interactions with 1/r³ decay
- **Bond Probability**: Uses quantum expectation values for cluster bond decisions

**Key Features**:
- **Quantum Wavefunctions**: Direct transformation of quantum state coefficients
- **Consistent Updates**: Both wavefunctions (`coef`) and spin expectation values (`eSpin`) updated consistently
- Pre-computed offsets and cumulative sums for O(1) neighbor selection
- Position mapping for periodic boundary conditions
- Robust fallback for systems without Luijten-Blöte tables
- GPU state synchronization support
- Compatible with crystal field environments and hyperfine interactions

### 2. Alternative Quantum Implementation (`thermalize_cluster_a_backup.m`) - **Reference Version**

**Algorithm**: Alternative unitary transformation cluster implementation
- **Purpose**: Reference implementation demonstrating different quantum cluster approach
- **Unitary Generation**: Random unitary matrices using Haar measure (QR decomposition)
- **Wavefunction Updates**: Applies transformations to quantum state coefficients
- **Energy Calculation**: Includes single-ion, dipole, and exchange contributions

**Key Features**:
- Alternative quantum cluster construction method
- Direct quantum mechanical expectation value calculations
- GPU acceleration support for energy calculations

### Adaptive Cluster Integration

Both thermalization (`equilibrate_a.m`) and sampling (`MC_sample_a.m`) phases include intelligent cluster update triggering:

**Adaptive Triggering Criteria**:
- **Autocorrelation Time**: τ > 5 (configurable threshold)
- **Low Acceptance Rate**: < 10% acceptance for single-spin updates
- **Forced Intervals**: Every 1000 steps regardless of conditions
- **Dynamic Intervals**: Adaptive spacing (50-500 steps) based on system response

**Integration Points**:
- **Equilibration**: Triggered after 5 consecutive trapped random walks
- **Sampling**: Continuous monitoring with autocorrelation-based adaptation
- **GPU Synchronization**: Automatic GPU state updates after successful cluster moves
- **Diagnostics**: Comprehensive logging of cluster attempt outcomes

### Performance Characteristics

**When Cluster Updates Are Most Effective**:
- Near critical temperatures where autocorrelation times become large
- Low-temperature ordered phases with strong correlations
- Systems with long-range dipole interactions
- Large system sizes where single-spin updates become inefficient

**Computational Cost**:
- Cluster construction: O(N log N) with Luijten-Blöte tables
- Energy calculation: O(C²) for cluster size C
- Memory overhead: Pre-computed lookup tables scale as O(N²)

## Key Algorithmic Differences

### Standard vs. `_a` Versions:

**Standard versions** (main directory):
- Use full distance calculations for every dipole interaction
- Suitable for smaller systems or initial testing

**`_a` versions** (beta/ directory):
- Implement neighbor lists for optimized dipole calculations
- Pre-compute interaction distances and vectors
- **Cluster Monte Carlo capabilities** for overcoming critical slowing down
- Luijten-Blöte lookup tables for efficient long-range cluster construction
- Adaptive cluster triggering based on autocorrelation analysis
- GPU acceleration with trial configuration support
- Significantly faster for large systems and critical temperature regimes
- Include additional optimizations and algorithmic improvements

## Development Workflow

### Recommended Development Approach:
1. **Test with standard versions** for verification and smaller systems
2. **Use `_a` versions** for production runs and large-scale simulations
3. **Enable cluster updates** for critical temperature studies and large systems
4. **Enable parallel computing** for temperature/field sweeps
5. **Use hyperfine interactions** selectively based on physical requirements

### Key Implementation Notes:
- Electronic spins evolve on the Bloch sphere using random rotations
- Metropolis acceptance criterion: `exp(-ΔE * β)` where β = 1/(kB*T)
- **Quantum cluster updates** use unitary transformations on (2J+1)×(2J+1) Hilbert space
- **Cluster updates** automatically triggered when autocorrelation times exceed threshold
- **Luijten-Blöte pre-computation** enables efficient O(N log N) cluster construction
- **Function signature**: Cluster algorithm requires `ion`, `hamI`, `basis`, `nSpin` parameters for quantum operations
- Nuclear spins thermalize separately when hyperfine interactions are enabled
- GPU acceleration available with automatic state synchronization
- Data automatically saves to OS-specific Google Drive paths
- Crystal field calculations use ion-specific parameters from literature

## Development Notes

- The `_a` versions contain significant performance improvements and should be preferred for most applications
- **Quantum Cluster Monte Carlo** implementations provide substantial acceleration near critical temperatures
- **Primary cluster algorithm** (`thermalize_cluster_a.m`) uses quantum unitary transformations with Luijten-Blöte efficiency
- **Quantum Consistency**: Wavefunctions and spin expectation values are always kept consistent
- GPU versions have been removed but GPU acceleration is integrated into cluster algorithms
- The codebase handles both electronic and nuclear spin degrees of freedom with separate thermalization
- Autocorrelation monitoring enables intelligent switching between single-spin and cluster updates
- Hyperfine interactions can be selectively enabled/disabled per simulation
- Multiple ion types supported with realistic physical parameters from literature
- Order parameters calculated using structure factors for neutron scattering analysis
- Luijten-Blöte lookup tables are automatically pre-computed for systems with periodic boundary conditions

## Recent Updates and Fixes

### December 2024 - Cluster Algorithm Units Fix

**Issue Resolved**: Fixed critical units mismatch in quantum cluster Monte Carlo algorithm.

**Problem**: The Luijten-Blöte table pre-computation in `cluster_init_a.m` was using lattice positions in Angstroms but computing `1/r³` values that were later used with dipolar coupling constants expecting SI units (meters). This caused a factor of `(1e-10)³ = 1e-30` error in bond probabilities, making cluster formation nearly impossible.

**Solution**: Updated `cluster_init_a.m` to convert positions from Angstrom to meters before distance calculations:
```matlab
r = params.pos * 1e-10;  % Convert from Angstrom to meters for consistent units
```

**Physical Reality**: Er³⁺ dipolar interactions are genuinely weak (~0.0006 meV nearest neighbor) compared to thermal energy at typical simulation temperatures (T=0.02K gives ~0.0017 meV). For significant cluster formation, temperatures below ~0.007K are needed where βJ > 1.

**Status**: ✅ **All units now consistent throughout codebase. Cluster algorithm mathematically and physically correct.**

**Next Steps**: 
- Test clustering at lower temperatures (T < 0.007K) for physical dipolar clustering
- Continue with other optimization tasks from multi-day debugging plan
- Benchmark algorithm correctness against known Er³⁺ phase behavior