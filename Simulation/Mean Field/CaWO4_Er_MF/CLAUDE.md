# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a MATLAB codebase for mean field theory simulations of magnetic properties in CaWO4:Er³⁺ (calcium tungstate doped with erbium). The code calculates eigenvalues, eigenvectors, energy transitions, and magnetic susceptibility under applied magnetic fields.

## Core Architecture

### Main Scripts
- `MF_Er_CaWO4.m` - Primary simulation script with modular function structure
- `MF_Er_CaWO4_v1a.m` - **OPTIMIZED VERSION** with consolidated Hamiltonian assembly, vectorized calculations, and enhanced ZEFOZ analysis
- `MF_Er_CaWO4_v1b.m` - Additional variant

### Core Functions
- `spin_operators.m` - Generates electronic (J) and nuclear (I) spin operators, including hybrid operators for hyperfine-coupled systems
- `cf.m` - Crystal field Hamiltonian using Stevens operators (O₂₀, O₄₀, O₄₄, O₆₀, O₆₄)
- `MF_chi.m` - Magnetic susceptibility calculation using linear response theory
- `MF_susceptibility.m` - **ENHANCED** susceptibility calculation with dual calling conventions and optimized parallel processing
- `dE_overlay.m` - Transition energy analysis and plotting
- `perturb2nd.m` - **ENHANCED** standalone second-order perturbation theory engine with advanced degeneracy handling

### ZEFOZ Analysis Suite (New)
- `compute_all_S1S2.m` - Comprehensive S1/S2 sensitivity calculation for all transition orders
- `compute_all_S1S2_along_direction.m` - Directional sensitivity analysis along specified field directions
- `plot_S1S2_top_zefoz.m` - Multi-page visualization of top ZEFOZ candidates with tensor heatmaps

### Optimized Architecture (v1a)
The optimized version features significant structural improvements:

#### Consolidated Functions
- `buildHamiltonian()` - Single function for all Hamiltonian assembly (electronic + nuclear + crystal field terms)
- `calcDerivatives()` - Unified derivative calculation for both uniform and non-uniform field grids
- `calcPertCorr()` - Vectorized perturbation theory calculations (eliminates nested loops)

#### Pre-computed Data Structures  
- Spin operators calculated once in `initialization()` and passed through function calls
- Shared `spinOps` structure eliminates redundant operator generation
- Reduced function parameter passing overhead

### Key Physical Parameters
The simulations use either:
1. Full crystal field approach (J=15/2 for Er³⁺)
2. Effective doublet model (J=1/2)

Both approaches include nuclear hyperfine coupling (I=7/2) and Zeeman effects.

### Calculation Flow
1. **Initialization**: Define physical constants, crystal field parameters, hyperfine coupling constants
2. **Eigenvalue Calculation**: Solve Hamiltonian at each magnetic field point
3. **Transition Analysis**: Calculate inter-level transitions and derivatives
4. **Perturbation Theory**: Second-order analysis with degeneracy-aware calculations
5. **ZEFOZ Analysis**: Comprehensive S1/S2 sensitivity analysis for field-independent transitions
6. **Susceptibility**: Compute frequency-dependent magnetic response using Kubo formula (optional)

### Data Structures
- `Options` struct controls calculation parameters (energy levels, transitions, temperature, frequency range)
- `const` struct contains physical constants
- `params` struct holds system parameters (J, I, g-factors, hyperfine constants)
- Results stored as multi-dimensional arrays indexed by [state, field_point, frequency]

### Plotting and Visualization
- Energy level diagrams vs magnetic field
- Susceptibility tensor components (χ_xx, χ_yy, χ_zz) as 2D color plots
- Transition frequency plots with customizable line styles and colors
- **Enhanced perturbation theory visualizations (v1a):**
  - Second-order perturbation analysis with 2×2 subplot layout (4 focused plots)
  - Generic `plotTransitionResponses()` function supports four analysis types:
    - Linear transition responses (dν/dB field derivatives)  
    - Quadratic transition responses (d²ν/dB² second derivatives)
    - **NEW**: First-order energy correction differences per transition
    - **NEW**: Second-order energy correction differences per transition
  - Automatic multi-figure generation (max 6 subfigures per figure)
  - Smart figure positioning prevents window overlap

## Development Guidelines

### Modifying Parameters
- Crystal field parameters in `initialization()` function
- Magnetic field range and orientation in main scripts
- Temperature and frequency grids in Options structure

### Adding New Calculations
- Follow the modular function structure in `MF_Er_CaWO4.m`
- Use vectorized MATLAB operations for field sweeps
- Store intermediate results for reanalysis without recalculation

### Performance Considerations
- Large parameter sweeps can be memory-intensive
- **v1a optimizations provide significant performance improvements:**
  - Consolidated Hamiltonian assembly reduces ~50 lines of duplicate code
  - Vectorized perturbation theory (O(n²) → O(1) operations)
  - Pre-computed spin operators eliminate repeated calculations
  - Unified derivative calculations reduce code complexity
- Progress indicators included for long calculations
- Parallel computing support for field sweeps (when Parallel Computing Toolbox available)

### Code Optimization Guidelines (v1a)
When modifying the optimized version:

#### Core Design Principles
1. **Consistent Function Names**: Use concise, descriptive names with standardized conventions
   - `get_GHz2meV()` for unit conversion (preferred naming convention)
   - `compute_all_S1S2()` for comprehensive sensitivity analysis
   - `plot_S1S2_top_zefoz()` for ZEFOZ visualization
   - `buildHamiltonian()` for Hamiltonian assembly

2. **Memory-Efficient Design**: Optimize for large-scale calculations
   - Pre-filtering with `Options.onlyZefoz` for ZEFOZ-focused analysis
   - Transition limits with `Options.maxTransitions` to prevent memory overflow
   - Smart memory cleanup in plotting functions for unused data

3. **Robust Error Handling**: Comprehensive validation and fallbacks
   - Unit conversion validation with `get_GHz2meV()` and CODATA standards
   - Complete helper function dependencies (no missing functions)
   - Graceful degradation when optional features unavailable

4. **Visual Consistency**: Professional presentation standards
   - Fixed 2×2 grid layouts for all multi-panel figures
   - Consistent color schemes and axis labeling
   - Optimized defaults (TopN=8) for balanced figure layouts

5. **Minimal Command Window Output**: Essential information only for performance and readability
   - Display only critical results, summaries, and error messages
   - Avoid verbose diagnostic output during normal operation
   - Prioritize runtime speed and clean command window presentation
   - Use concise status messages for long calculations

#### Function Design Principles
- **Single Responsibility**: Each function handles one specific physics calculation
- **Data Structure Reuse**: Pass pre-computed structures (`spinOps`, `params`, etc.) rather than recalculating
- **Vectorization**: Prefer matrix operations over loops for MATLAB performance
- **Memory Efficiency**: Use logical indexing and masking for conditional operations

#### Hamiltonian Construction
- Use `buildHamiltonian()` for all Hamiltonian assembly
- Electronic Zeeman: `muB * gE * J * B`  
- Nuclear terms: `A*I*J + Q*I² + μN * gN * I * B`
- Crystal field: Pre-computed `hamCF` matrix

#### Derivative Calculations  
- Use `calcDerivatives()` for both first and second derivatives
- Handles uniform/non-uniform field grids automatically
- Vectorized finite difference schemes for accuracy and speed

#### Perturbation Theory (Enhanced Standalone Module)
**NEW: `perturb2nd.m` - Professional-grade perturbation theory engine:**
- **Rigorous degenerate perturbation theory** with basis-consistent block tensors
- **Direction parameter**: Compute perturbations along specific field directions `[1,0,1]`
- **Advanced output structure**: `hasDegeneracy`, `blocks`, `blockK` with common basis matrices
- **Automatic degeneracy detection**: 9 degenerate manifolds identified in CaWO4:Er³⁺
- **Basis consistency fix**: Block tensors use single `Ublk` basis per manifold (eliminates cross-term artifacts)
- **Dual computation modes**: Full tensor `K_ij` (non-degenerate) or diagonal-only (degenerate systems)
- **Wavefunction corrections**: First-order `|ψ_n⁽¹⁾⟩` with proper basis rotation
- **Professional documentation**: Comprehensive theory, examples, and implementation notes

#### Visualization Improvements (v1a)
- **Main perturbation figure**: Streamlined 2×2 layout with enhanced degeneracy handling
  - Plot 1: First-order energy corrections (grouped bar chart by field direction)
  - Plot 2: Second-order corrections with **adaptive display** (tensor vs diagonal mode)
  - Plot 3: **NEW**: Ground state wavefunction probability changes due to B_z perturbation
  - Plot 4: Combined field sensitivity with dual y-axes (linear + quadratic responses)
- **Generic transition plotting**: Single `plotTransitionResponses()` function handles all plot types
  - Supports: 'linear', 'quadratic', 'first_order_corr', 'second_order_corr'
  - **New perturbation correction analyses**: Energy correction differences per transition
  - Maximum 6 subfigures per figure prevents overcrowding
  - Automatic multi-figure generation with smart positioning
  - Eliminated ~340 lines of duplicate plotting code

#### Implementation Examples
These principles are demonstrated throughout the v1a optimizations:

**Modular Consolidation Success:**
- `plotTransitionResponses()` - Single 200-line function replaces two 170-line functions
- `buildHamiltonian()` - Consolidates electronic, nuclear, and crystal field terms  
- `calcPertCorr()` - Vectorized matrix operations eliminate nested loops

**Balanced Abstraction Examples:**
- Nested plotting functions (`plotLinearData`, `plotQuadraticData`) within main function
- Inline parameter validation instead of separate checker functions
- Direct mathematical operations in main code for simple unit conversions

**Naming Convention Standards:**
- `const.Gh2mV` - Preferred field name for GHz→meV conversion factor
- `get_GHz2meV()` - Robust unit conversion with validation 
- `GHz2meV` - Standard variable name for conversion factor
- `compute_all_S1S2()` - Comprehensive analysis functions with descriptive names

## ZEFOZ Analysis Suite

### Overview
The ZEFOZ (Zero-Field-Zero-Frequency) analysis suite provides comprehensive tools for identifying and analyzing field-independent quantum transitions. These transitions are critical for quantum sensing applications due to their immunity to magnetic field fluctuations.

### Key Components

#### `compute_all_S1S2.m` - Comprehensive S1/S2 Sensitivity Analysis
**Purpose**: Calculate first-order (S1, GHz/T) and second-order (S2, GHz/T²) field sensitivities for all requested transition orders.

**Key Features**:
- **Multi-order transitions**: Supports `Options.ndE` for orders 1, 2, 3, ..., or custom transition pairs
- **Dual data sources**: Uses both full tensor results (non-degenerate) and axis-wise results (degenerate-robust)
- **Automatic unit conversion**: Converts from meV to GHz using `const.Gh2mV` via robust `get_GHz2meV()`
- **True curvature calculation**: S2 represents actual Hessian (2× stored coefficient)
- **Memory optimization**: Pre-filtering with `Options.onlyZefoz` and transition limits
- **Sorted output**: Results ranked by ZEFOZ potential (ascending S1norm, S2fro)

**Outputs**:
- `T`: Table with columns `m`, `n`, `order`, `freq_GHz`, `S1norm_GHzT`, `S2fro_GHzT2`, plus axis components
- `Trans`: Struct array with detailed per-transition data including full tensors when available

#### `compute_all_S1S2_along_direction.m` - Directional Sensitivity Analysis  
**Purpose**: Calculate S1/S2 sensitivities along specific magnetic field directions.

**Key Features**:
- **Direction specification**: Takes unit vector `uhat` for field direction
- **Flexible transition selection**: 'adjacent', 'all', or custom pairs  
- **Robust calculation**: Works with both tensor and axis-wise perturbation results
- **Complete implementation**: Includes all helper functions (`transition_S1S2_direction`, `transition_S1S2_axis`)

#### `plot_S1S2_top_zefoz.m` - Multi-Page ZEFOZ Visualization
**Purpose**: Professional visualization of top ZEFOZ candidates with comprehensive analysis panels.

**Key Features**:
- **Consistent 2×2 layout**: Fixed 4 candidates per figure (2×2 grid) for visual consistency
- **Optimized defaults**: TopN=8 produces exactly 2 figures with balanced layouts
- **Dual-axis plots**: S1 (bar chart, left y-axis) + S2 (line plot, right y-axis) per panel
- **Optional tensor heatmaps**: Separate figures showing full S2 tensors with robust scaling
- **Noise analysis**: Optional linewidth prediction when `predict_linewidth_from_S1S2` function available
- **Smart scaling**: Panel-wise or global color limits to prevent "solid color" artifacts

### Integration with v1a

The ZEFOZ suite is fully integrated into `MF_Er_CaWO4_v1a.m`:

```matlab
% Enhanced ZEFOZ analysis with top candidates visualization (2×2 layout consistency)
[T_top, idxTop] = plot_S1S2_top_zefoz(pertResults, const, Options, ...
    'TopN', 8, 'NoiseCov', NoiseCov);  % 8 candidates = 2 figures × 4 panels each

% Alternative: comprehensive S1/S2 calculation (commented)  
% [T_all, Trans] = compute_all_S1S2(pertResults, const, Options);
```

### ZEFOZ Theory

**S1 (Linear Zeeman Shift)**: `S1 = ∂ν/∂B = (dE₁ᵐ - dE₁ⁿ)/ℏ`
- For ZEFOZ: |S1| ≈ 0 (field-independent frequency to first order)

**S2 (Quadratic Zeeman Shift)**: `S2 = ∂²ν/∂B² = 2(dE₂ᵐ - dE₂ⁿ)/ℏ`  
- Controls second-order field sensitivity and linewidth broadening

## Advanced Perturbation Theory Module

### `perturb2nd.m` - Enhanced Second-Order Perturbation Theory

This standalone module provides rigorous second-order perturbation theory for arbitrary quantum systems with special attention to degenerate manifolds.

#### Key Features

**Core Capabilities:**
- **Arbitrary Hamiltonians**: Works with any Hermitian matrix (not limited to spin systems)
- **Multi-component perturbations**: Handles vector perturbations `{V_x, V_y, V_z}` for magnetic field studies
- **Direction analysis**: Compute perturbations along specific directions `[u_x, u_y, u_z]`
- **Full tensor support**: Second-order tensor `K_{ij}` for cross-term analysis

**Degeneracy Handling:**
- **Automatic detection**: Groups energy levels within configurable threshold
- **Block diagonalization**: First-order corrections via proper degenerate PT
- **Basis-consistent tensors**: Uses common `U_blk` basis for all perturbations in each manifold
- **Cross-term integrity**: Eliminates basis-mixing artifacts in `K_{ij}` calculations

**Advanced Outputs:**
```matlab
results.eigE0        % [N x 1] Sorted unperturbed energies
results.eigW0        % [N x N] Corresponding eigenvectors
results.dE1          % [N x k] First-order corrections
results.dE2          % [N x k x k] Full tensor OR [N x k] diagonal
results.dPsi1        % {k}[N x N] First-order wavefunction corrections
results.hasDegeneracy% Logical degeneracy flag
results.blocks       % Cell array of degenerate manifold indices
results.blockK       % Block-level effective tensors with basis info
```

#### Usage Examples

**Basic usage:**
```matlab
results = perturb2nd(H0, V);                    % Single perturbation
results = perturb2nd(H0, {Vx, Vy, Vz});         % Vector perturbations
```

**Advanced features:**
```matlab
% Compute along specific direction (eliminates degeneracy ambiguity)
results = perturb2nd(H0, {Vx, Vy, Vz}, 'direction', [1, 0, 1]);

% Configure degeneracy threshold and calculation options
results = perturb2nd(H0, perturbations, 'threshold', 1e-10, ...
                    'calcTensor', false, 'calcWaveFunc', true);
```

#### Physical Applications

**CaWO4:Er³⁺ System Results:**
- **9 degenerate manifolds** automatically detected (nuclear hyperfine + crystal field)
- **Linear Zeeman anisotropy**: 10¹⁴ ratio (extreme Ising behavior)
- **Quadratic Zeeman isotropy**: X-Y plane coupling ~10⁴ GHz/T²
- **Block tensors**: Provide manifold-specific effective Hamiltonians

**Magnetic Field Analysis:**
- **Linear responses**: `dE/dB` for first-order Zeeman shifts
- **Quadratic responses**: `d²E/dB²` for second-order corrections and field-dependent mixing
- **Wavefunction evolution**: Probability redistribution under perturbation

#### Mathematical Implementation

**Degenerate Perturbation Theory:**
- Within each degenerate block: `E_n^(1) = ⟨n|V|n⟩` after block diagonalization
- Cross-manifold coupling: `E_n^(2) = Σ_{m∉block} |V_{nm}|²/(E_n - E_m)`
- Block tensors: `K_{ij}^(g) = Σ_{m∉g} V_i(g,m) V_j(m,g)/(E_g - E_m)` in common basis

**Basis Consistency:**
- Common basis per block: `U_blk` from diagonalizing `Σ_i V_i(G,G)²`
- All perturbations rotated identically: `V_i^{rot} = U_blk† V_i U_blk`
- Cross-terms mathematically consistent: `K_{ij}` elements meaningful for `i≠j`

## File Dependencies
The code relies on standard MATLAB functions and does not require external toolboxes. The physics modules are self-contained within this directory.

### Testing and Validation
- `test_optimized.m` - Validation script for optimized functionality  
- Verifies numerical consistency between original and optimized versions
- Tests core functions with minimal parameter sets
- **NEW**: Basis consistency verification for degenerate perturbation theory

## Enhanced Susceptibility Module

### `MF_susceptibility.m` - Dual-Signature Susceptibility Calculation

**Purpose**: Compute frequency-dependent magnetic susceptibility tensor χ(ω) with backward compatibility.

**Key Features**:
- **Dual calling conventions**: Supports both new optimized signature and legacy interface
- **Automatic signature detection**: Intelligently detects calling pattern and converts parameters
- **Enhanced parallel processing**: Improved parallel efficiency with reduced overhead
- **Combined electronic-nuclear response**: Effective g-factor weighting for nuclear contributions
- **Vectorized frequency calculation**: 3D tensor operations for optimal performance
- **Professional visualization**: 2×2 subplot layout for real/imaginary χ_xx and χ_zz components

**Calling Conventions**:
```matlab
% New optimized signature (recommended)
results = MF_susceptibility(Opts, const, params, eigE, eigW, Bf, spinOps);

% Legacy signature (backward compatibility)  
results = MF_susceptibility(Options, const, J, I, eigenE, eigenW, Bfield);
```

**Technical Implementation**:
- **Kubo formula**: Linear response theory with Lorentzian broadening
- **Thermal populations**: Boltzmann distribution with temperature dependence
- **Broadening parameter**: Configurable γ (meV) for spectral linewidth
- **Nuclear coupling**: Weighted electronic + nuclear magnetic moment operators

### Recent Developments

#### Version History
**Latest (September 2025)**: ZEFOZ Analysis Suite & Enhanced Modules
- **ZEFOZ Analysis Suite**: Complete toolkit for field-independent transition analysis
  - `compute_all_S1S2.m`: Multi-order S1/S2 sensitivity calculations with dual data source support
  - `compute_all_S1S2_along_direction.m`: Directional sensitivity analysis capabilities  
  - `plot_S1S2_top_zefoz.m`: Professional multi-page visualization with tensor heatmaps
- **Enhanced susceptibility**: `MF_susceptibility.m` with dual calling conventions and optimized parallel processing
- **Integration improvements**: ZEFOZ suite fully integrated into v1a main script

**Previous (2025)**: Enhanced Perturbation Theory Module
- **Modular architecture**: Extracted `perturb2nd.m` as standalone professional-grade function
- **Absorbed advantages**: Combined best features from multiple implementations into single enhanced module
- **Bug fix**: Resolved basis inconsistency in degenerate block tensor calculations
- **Adaptive visualization**: Plotting functions handle both tensor and diagonal modes automatically
- **Wavefunction visualization**: Added ground state probability redistribution plots

**Key Improvements:**
- **ZEFOZ-first workflow**: Prioritized field-independent transition identification
- **Consistent visualization**: Fixed 2×2 grid layouts with TopN=8 default for balanced presentation
- **Memory optimization**: Smart pre-filtering and transition limits for large systems
- **Robust unit conversion**: Standardized `get_GHz2meV()` with validation and fallbacks
- **Complete function dependencies**: All helper functions implemented and tested
- **Direction parameter**: New capability for directional perturbation analysis
- **Professional documentation**: Comprehensive theory guide with examples
- **Enhanced degeneracy handling**: 9 degenerate manifolds properly identified and processed
- **Basis-consistent cross-terms**: Block tensors now mathematically rigorous

#### Known Issues (Resolved)
- ✅ **Basis mixing bug**: Fixed inconsistent rotation matrices in degenerate block tensors
- ✅ **Tensor dimension errors**: Added adaptive handling for calcTensor true/false modes  
- ✅ **Empty subplot problem**: Resolved zero-valued quadratic corrections display
- ✅ **Limited transition plotting**: Extended to full 15-order analysis with multi-figure layout
- ✅ **Missing function dependencies**: Implemented `transition_S1S2_axis()` helper function
- ✅ **Memory efficiency issues**: Added pre-filtering and transition limits for large systems
- ✅ **Unit conversion inconsistencies**: Standardized to `get_GHz2meV()` with robust validation
- ✅ **Inconsistent figure layouts**: Fixed ZEFOZ visualization to consistent 2×2 grids