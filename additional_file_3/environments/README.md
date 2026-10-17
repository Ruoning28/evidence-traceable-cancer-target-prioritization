# Environment separation

Root `requirements.txt` and `environment.json` retain the analysis environment. The plotting locks here specify the separately verified Windows x86-64 rendering environment. Do not merge these pins into the analysis environment.

- `figures-win-64.explicit.txt`: complete transitive Conda dependency closure, exact URLs/builds and MD5 package checksums.
- `figures-conda-packages.json`: package versions/builds, retained SHA-256/MD5 and licence metadata.
- `figure2-pip-win64.lock.txt`: SHA-256-locked Matplotlib 3.10.0 CPython 3.13 Windows x86-64 wheel; install with `--no-deps` into a separate target on the locked plotting runtime.
- `plotting-profiles.json`: Python/library/FreeType/font fingerprints checked before drawing.
- `figure-reference.json`: manuscript native drawing-stream and RGB hashes, excluding creation timestamps.

Rebuild and per-figure commands: `../figure_code/README.md`. Network access is needed unless the exact packages are cached. Locks use public package URLs, without local author-directory paths. Arial is checked but not redistributed.

The profile is platform-specific; pixel equality is not claimed for another OS, architecture, font or FreeType build. Rendering validation does not mean all upstream models were refitted.
