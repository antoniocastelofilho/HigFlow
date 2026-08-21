# Third-party sources

Dependencies that are built from source rather than installed from a package
manager. Nothing here is part of HigFlow.

| File | What it is |
|---|---|
| `fetch-petsc.sh` | Downloads and verifies the PETSc 3.14.0 release archive |
| `libfyaml-master.zip` | libfyaml source snapshot |

## PETSc

The release archive is **not** stored in this repository. Fetch it with:

```bash
bash bibliotecas/fetch-petsc.sh
```

The script downloads `petsc-3.14.0.tar.gz` from the Argonne release snapshots
and checks it against a recorded SHA-256, so a corrupted or substituted
download fails loudly instead of producing a subtly different build. Running it
again when a valid archive is already present does nothing.

Both `install_higflow_ubuntu22` and the Singularity build script call it before
extracting.

## libfyaml

`libfyaml-master.zip` is still committed. It is a snapshot of the upstream
`master` branch rather than a tagged release, so there is no stable URL that
would reproduce this exact archive — replacing it with a download means also
choosing a release version to pin to, and that changes which library the
project builds against. That is a dependency decision for the maintainers, not
a cleanup, so the file stays until it is made.

## Ignored contents

`.gitignore` excludes downloaded archives (`*.tar.gz`, `*.part`) and the
directories they extract into (`petsc-*/`, `libfyaml-master/`). The directory
itself is tracked so that this file and `fetch-petsc.sh` are versioned.
