# Repository Guidelines

## Project Structure & Module Organization
The optimized workflow lives in `MF_Er_CaWO4_v1a.m`, which orchestrates initialization, Hamiltonian assembly, eigen-solvers, and ZEFOZ post-processing. `MF_Er_CaWO4.m` and `MF_Er_CaWO4_v1b.m` capture legacy and experimental variants; keep new studies close to whichever options block they extend. Reusable physics helpers sit in `spin_operators.m`, `cf.m`, `MF_susceptibility.m`, and `pert_2nd.m`, while ZEFOZ utilities (`compute_all_S1S2*.m`, `plot_S1S2_top_zefoz.m`) provide analysis and visualization. Store quick diagnostics in `test_v1a.m` and keep large figures or scratch notebooks outside the repo to avoid bloat.

## Build, Test, and Development Commands
Run the main simulation headless to surface runtime errors quickly:
```
matlab -batch "try, run('MF_Er_CaWO4_v1a.m'); catch ME, disp(getReport(ME)); exit(1); end; exit"
```
Swap scripts by changing the filename and tune the Options block before execution. For regression checks, call `matlab -batch "run('test_v1a.m')"` and expect the stage-by-stage log to complete without warnings. When iterating on ZEFOZ tooling, prepare inputs in a helper script and run `matlab -batch "run('plot_S1S2_top_zefoz.m')"` to validate plots.

## Coding Style & Naming Conventions
Stick to four-space indentation, uppercase `%%` section headers, and concise comments that explain physics intent. Use camelCase for functions like `buildHamiltonian` and for struct fields, while filenames stay domain-oriented (`MF_`, `compute_`, `pert_`). Constants belong in the `const` struct; new configuration values should extend `Options` so downstream helpers receive them automatically. Keep log strings ASCII to prevent encoding issues across platforms.

## Testing Guidelines
Maintain 	est_v1a.m as the smoke-test harness. Add focused sub-tests inside the 	ry block and reuse the light grids to keep runtime under a few minutes. When touching perturbation or susceptibility code, confirm pert_2nd.m returns tensors without warnings and capture any threshold tuning in comments. Record expected Zeeman ranges or ZEFOZ candidate counts when you change tolerances so regressions stay easy to spot.

## Commit & Pull Request Guidelines
Favor short, imperative commit subjects that highlight the touched module, e.g., `MF: tighten ZEFOZ filter` or `pert: cache spin operators`. Consolidate debug commits before sharing. Pull requests should summarize the physics motivation, list modified scripts, and include before/after metrics (runtime, ZEFOZ hits, susceptibility traces). Link to supporting notebooks or figures stored externally when referencing derived results.

## Agent Workflow Tips
Heavy parameter sweeps can exceed available memory; stage large field grids in scratch scripts and double-check `Options.parallel` before committing. Document new struct fields or data dependencies both at the call site and at the top of the helper you touched so future contributors can reproduce the setup.









