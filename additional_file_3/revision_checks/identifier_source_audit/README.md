# Identifier and exclusion-scope checks

From the Additional file 3 root, run:

```text
python revision_checks/identifier_source_audit/run_audits.py --data ../additional_file_2 --output ./identifier_audit_outputs
```

The checks read the packaged HGNC reference, retained UniProt name fields, frozen candidate matrices, labels, strata, scores, configurations and source code. They produce mappings, fixed-ranking and review-universe comparisons, primary sampling checks and source-exclusion scopes. No model fitting or resource download is performed. NumPy and pandas are required; use the supplied analysis environment.

The outputs distinguish direct feature descendants from shared-resource paths. Top-K checks use fixed scores. The label-origin, background and evaluation scope of each analysis remains explicit in the resulting tables.
