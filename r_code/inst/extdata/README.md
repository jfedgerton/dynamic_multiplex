# extdata

`stability_calibration_table.csv` — calibration lookup used by
`partition_stability()`. Columns: `level` (partition_nmi, partition_ari,
node_jaccard), `stab_lo`/`stab_hi` (stability bin), `n_calib` (calibration
fits in the bin), `acc_median`, `acc_q05` (the 5th-percentile accuracy floor),
`source`. Written by `replication/post/12_stability.R`; copied here by
`replication/run_all.sh post`. The bundled copy is provisional until the
front-to-end replication rerun regenerates it (partition_ari rows are added
then).
