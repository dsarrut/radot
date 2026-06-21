# radot Development Checklist

When implementing new tools, commands, or functions in `radot`, please adhere to the following checklist to ensure consistent handling of image metadata and parameters.

## Metadata Handling

- [ ] **Check for Metadata**: If a function depends on parameters such as image modality (CT, PET, SPECT), voxel physical unit (HU, SUV, Bq/mL, Bq), patient weight, or injected activity, check if these can be read from the image metadata header first.
- [ ] **Require Explicit Values as Fallback**: If the required metadata does not exist in the image header, the function **must** raise an explicit, clear error (e.g. `ValueError` or similar) informing the user that the value is missing.
- [ ] **CLI Configuration**: 
  - Do not mark CLI options like `--type`, `--unit`, `--weight`, or `--activity` as `required=True`.
  - Avoid hardcoding static defaults (like `default="ct"`) for parameters that could otherwise be resolved from metadata, as this prevents metadata resolution.
  - Let options default to `None` and resolve them dynamically from the image metadata inside the code.
- [ ] **Unit Tests**: Always include tests that verify:
  - Successful execution when the parameters are omitted on the command line but present in the image metadata.
  - Correct failure and explicit error messages when parameters are omitted and no metadata is present in the image.
