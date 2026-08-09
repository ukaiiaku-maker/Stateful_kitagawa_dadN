# Stateful-PD Kitagawa v2.8.1 summary-discovery hotfix

The v2.8 K360 failure calculation completed physically, but the gate runner predicted a sigma-directory tag with 12 significant digits while the solver wrote its established 6-digit display tag. The hotfix locates a completed summary by the immutable `sigma_a_MPa` value inside `summary.json`, requires exactly one match, and preserves the original v2.8 request identity so the completed job is reused without rerunning the solver.

The solver, parameter values, frozen v8.3 campaign, and physical output are unchanged.
