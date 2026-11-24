# Time integration on synthetic curves

Results: [docs/9_validation_time_integration.md](../../docs/9_validation_time_integration.md).

`radot integrate` on synthetic ¹⁷⁷Lu time-activity curves whose integral is known exactly: three shapes
(kidney, liver, lesion), three acquisition schedules (3, 4 and 6 time points), proportional noise (5 and
10 %), and the integration methods (fitted models with the noise estimated or known, model selection,
trapezoid, single time point). Per configuration: bias and spread of the TIA, failures, and the coverage
of the standard error. No data needed.

| File | Role |
|---|---|
| `tac_benchmarks.py` | Curves, schedules, methods, Monte Carlo, acceptance ranges (single source of truth) |
| `test_integrate_validation.py` | CI tests: exact without noise, true model unbiased with an honest standard error (300 realisations), page up to date |
| `results/integrate_results.json` | Full Monte Carlo (1000 realisations per configuration) |
| `report.py` | `--run`: full Monte Carlo (a few minutes); writes `docs/9_validation_time_integration.md` and its figures; `--check` (page text only) |
| `figures.py` | The figures of the page (matplotlib, needed only to regenerate them), `docs/images/validation_integrate/` |

```bash
pytest validation/integrate
python validation/integrate/report.py --run   # after a change of the fit, then commit results/ and the page
```
