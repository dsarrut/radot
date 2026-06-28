# RADOT Logging & Printing Guide

This document defines how diagnostic messages, operational progress logs, and warnings are managed and configured in the `radot` library.

## 1. Programmatic Python API

By default, the `radot` library initializes a package-level logger that sends messages to standard error (`sys.stderr`) at the `INFO` level. This ensures you see progress logs in your scripts without standard output (`sys.stdout`) pollution.

### Changing Log Level Programmatically

You can adjust the logging level directly via the package-level `set_log_level()` helper:

```python
import radot

# Silence all diagnostic logs (only show WARNING or ERROR)
radot.set_log_level("WARNING")

# Set to DEBUG to see detailed resampling and database downloads
radot.set_log_level("DEBUG")

# Restore standard level
radot.set_log_level("INFO")
```

### Accessing the Package Logger

If you want to attach custom handlers (e.g. logging to a file), you can retrieve the standard logger directly:

```python
import logging
import radot

logger = logging.getLogger("radot")
logger.addHandler(logging.FileHandler("dosimetry.log"))
```

---

## 2. Command-Line Interface (CLI)

When running the `radot` command-line utility, all progress logs and status updates are printed to `sys.stderr` so that standard output remains clear for numerical piping or redirection.

### CLI Verbosity Controls

You can control verbosity using the global flags:

```bash
# Run silently (only output error reports)
radot dose spect.nii.gz -l labels.nii.gz -m hanscheid2018 --isotope lu177 --quiet

# Run with verbose debugging details (forces logger to DEBUG level and shows full python tracebacks on failures)
radot dose spect.nii.gz -l labels.nii.gz -m hanscheid2018 --isotope lu177 --verbose
```
