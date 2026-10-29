# Current submission manifests

These manifests describe the final revised manuscript and Additional files 1–4. `submission_artifacts.tsv` records the exact submitted artifact identities; the reference manuscript PDF is not a newly published article or a new analysis run.

Extract Additional files 2 and 4 into the same parent directory. They restore one `additional_file_2` folder with 180 files. Its `sha256.json` (copied verbatim as `analysis_inputs_sha256.json`) verifies the other 179 files, including the two matrices in Additional file 4. `data_archive_members.tsv` maps every data file to its actual ZIP and digest. The checksum manifest belongs to Additional file 2, not the code archive.

`source_metadata/` and `source_roles/` are exact copies of selected metadata and role files read from Additional file 2. The source-data overview retains their documented resource descriptions. F/L/V assignments depend on the cancer, axis and record; broad source-family membership does not imply independent measurements. PDAC reproducibility support uses GEO cross-cohort support; the LUAD retained field is an internal TCGA/DepMap composite and is not independent external replication. The exact field definitions are in `source_metadata/reproducibility_support_field_dictionary.tsv`.

The data archives and upstream raw downloads remain distributed through the manuscript supplementary files and documented original providers. The repository publishes the matching analysis code archive and these retrieval/identity manifests. No new full-model fit, full raw-source reconstruction or new resource acquisition is claimed by this synchronization.
