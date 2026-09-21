# DynMux / Political Analysis: finalize checklist (target: submit Thu Sept 24)

State as of Sun Sept 20, ~10 pm ET. Repo `jfedgerton/dynamic_multiplex` at `fdac30f` (+ two one-line
patches on Roar: post/11 row-count assertion, empirical/10 bootstrap guard; will be committed and pushed
from Roar before the rerun is submitted). Everything below refers to that repo, not the Dropbox copy
(`/dynamic_multiplex` on Dropbox is the July 1.1.0 pipeline and must not be used).

## 0. What is settled (do not reopen)

| Decision | Outcome |
|---|---|
| Package bugs | weighted Jaccard (non-member strengths) and multislice null model fixed; v1.3.0; regression tests pass |
| Co-membership CI section | dropped; replaced by stability score + calibrated accuracy floor (pre-registered rule PASSED on Sept 20 run) |
| Overlap coupling | one sentence in main text + appendix row (K MAE); no overlap arm anywhere else |
| Second-stage resolution | default (gamma2 = 1); not a parameter |
| DynMux vs multislice | DynMux wins births/deaths (+0.06), abrupt rewiring (+0.10), node turnover (+0.08), robust across omega; multislice wins gradual switching, recurrence with period links, core/periphery, long series; Hungarian never wins on accuracy |
| Decision rule | one main-text subsection "When to use what" with `tab_decision_tree.tex` (post/16); mechanism detail in appendix |
| Multislice rows in Table 2 | adjacent links AND "same links as DynMux" |

## 1. Monday morning: rerun results (Claude)

- [ ] Roar rerun finished (`bash replication/run_all.sh status`); post-processing job exit 0
- [ ] `post/12` DECISION line = PASS for partition NMI (and ARI); if FAIL, stop and discuss before any writing
- [ ] Table 2 (`tab_metrics_wide.tex`): DynMux rows reproduce Sept 20 within seed noise; multislice rows updated
- [ ] `tab_decision_tree.tex`: DynMux leads on births/deaths, abrupt rewiring, node turnover survive
- [ ] Empirical stability table (`tab_app_empirical_stability.tex`) present for all 4 networks x 2 methods
- [ ] Calibration table copied into `r_code/inst/extdata` and `python_code/.../data`; packages reinstalled; regression tests pass
- [ ] Commit + push from Roar; Mac `git pull`

Deliverables to Jared by early afternoon:
- [ ] `main_edited.tex`, `appendix_edited.tex` regenerated with final numbers in every `%% [CLAUDE ...]` bullet
- [ ] `decision_subsection_draft.tex`: the "When to use what" subsection + table include (draft for you to adopt or rewrite)
- [ ] `tables_figures_v130.tgz`: all `manuscript/tables/*.tex` and `manuscript/figures/*` for Overleaf upload

## 2. Monday: Overleaf setup (Jared)

- [ ] Upload `main_edited.tex` as `main.tex` and `appendix_edited.tex` as `appendix 1.tex` (or diff them in first)
- [ ] Upload the new `tables/` and `figures/` (replace old ones; delete `fig_coverage_curve.png`, `fig_calibrated_interval.*`, `tab_app_coverage*`)
- [ ] Confirm `\input{tables/tab_decision_tree.tex}` and `figs/fig_stability_floor.pdf` compile
- [ ] Remove `app_coverage_robust.tex` from the project if it is still included anywhere

## 3. Monday-Tuesday: prose (Jared), in this order

1. [ ] Section 4 (uncertainty): stability score s, calibrated floor, decided/undetermined pairs, node stability. Bullets in the .tex give every number. Old CI claims are gone; do not carry any over.
2. [ ] Results text for Table 2: multislice numbers changed; "same links" row added; overlap row gone (one sentence pointing to appendix)
3. [ ] "When to use what" subsection: adopt or rewrite the draft; keep the asymmetric claim (DynMux when the node set or community set changes; multislice when nodes persist and per-layer signal is weak; Hungarian only when per-layer partitions must be preserved; unknown regime -> compare stability scores)
4. [ ] Abstract, intro, conclusion: replace any "DynMux dominates" framing with the asymmetric claim
5. [ ] Fix the truncated "among oth...oaches" sentence (flagged in main_edited.tex)
6. [ ] Fix eq:wjaccard (printed form is the pre-fix bug; bullet gives the corrected definition)
7. [ ] Appendix: new sections already scaffolded with headers + table includes (Jaccard vs overlap; omega sweep; selection rule; mechanism tests; empirical stability; leave-one-level-out). One paragraph each.
8. [ ] Placeholders: "Appendix X"/"XX" refs, Figure 3 include name, Warsaw Pact ccode-100 decision

## 4. Tuesday-Wednesday: verification (Claude, after your prose)

- [ ] Copy-edit pass (grammar/spelling only; cross-refs; `\ref`/`\cite` resolve)
- [ ] Every number in the text matches `manuscript/tables/*.tex` (list mismatches)
- [ ] No sentence still claims DynMux beats multislice in general
- [ ] README + `run_all.sh` instructions verified on a clean clone; `R CMD check --as-cran`; Python build
- [ ] Optional: Codex independent review (items 1-5 from the list I gave you) on a fresh clone, Monday; manuscript number check Wednesday

## 5. Thursday: submit (Jared)

- [ ] Final compiled PDF read-through
- [ ] GitHub tag `v1.3.0`; replication tarball from the tagged commit
- [ ] CRAN / PyPI release: your call whether before (Tue) or after submission; paper can cite the GitHub tag either way
- [ ] Submit

## Standing rules for everything above
- Never cancel Roar jobs; never rewrite Jared's prose (headers, refs, grammar/spelling, `%%` bullets only)
- Seeds: 123 for grids; per-task seeds documented in each script header
- Any new bug in the package or metrics found before Monday noon: report immediately; that is the only thing that would force a second rerun
