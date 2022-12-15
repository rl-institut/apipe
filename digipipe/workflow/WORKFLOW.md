# Workflow

**TBD**

## Snakefiles and config

- The global workflow is defined in the main
  [Snakefile](../workflow/Snakefile).
- It includes the module Snakefiles from the data store located at
  - [store/preprocessed/module.smk](../store/preprocessed/module.smk) and
  - [store/datasets/module.smk](../store/datasets/module.smk)
- In each of these modules, the rules as well as the config from the contained
  datasets are imported.
  