Fixes #

## Before merging into `dev`-branch, please make sure that the following points are checked:

- [ ] All pre-commit tests passed locally with: `pre-commit run -a`
- [ ] File `CHANGELOG.md` was updated
- [ ] The docs were updated

If packages were modified:
- [ ] File `poetry.lock` was updated with: `poetry lock`
- [ ] A new env was successfully set up

If data flow was adjusted:
- [ ] Data pipeline run finished successfully with: `snakemake -jX`
- [ ] Esys appdata was created successfully with: `snakemake -jX make_esys_appdata`

  (with `X` =  desired number of cores, e.g. 1)
