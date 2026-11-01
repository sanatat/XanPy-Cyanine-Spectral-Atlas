# Archived exploratory cyanine structural matching

These scripts correspond to earlier stages 16–17. They require the cleaned FluoDB snapshot created on the original user's computer, **not distributed here**. Their strict chain bridge and terminal-group classifiers can exclude chemically valid cases and must not be interpreted as exhaustive database annotations. The most recent all-tag user-run screening found **zero qualifying matched single-substituent pairs** *under these algorithmic rules*. This does **not** falsify the 2026 substituent rule.

The Stage17 all-tags script imports Stage16 by path. Pass it explicitly:

```bash
python workflows/cyanine_identity/xanpy_stage17_alltags_xi_test.py \
  --stage16-script workflows/cyanine_identity/stage16_direct_xi_test.py \
  --help
```

See `README_stage16_source.md` for the upstream stage16 experimental protocol. These scripts are archived for transparency, **not run in the default CI or claimed to identify all cyanines without manual structure review**.
