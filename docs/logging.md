# RADOT Logging & Printing Guide

This document defines how diagnostic messages, operational progress logs, and warnings are managed and configured in the `radot` library.

## 1. Programmatic Python API

By default, the `radot` library initializes a package-level logger that sends messages to standard error (`sys.stderr`) at the `INFO` level. This ensures you see progress logs in your scripts without standard output (`sys.stdout`) pollution.

### Changing Log Level Programmatically

You can adjust the logging level directly via the package-level `set_log_level()` helper:

```python
import radot

# Silence all diagnostic logs and standard warnings (only show CRITICAL/ERROR logs)
radot.set_log_level("QUIET")

# Alternatively, set standard python logging levels
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

You can control verbosity using the following options:

*   **Default (`INFO` level):** Shows standard operational logs (like reading images, conversion stages, and database downloads) printed to `sys.stderr`, followed by the final dosimetry table on `sys.stdout`.
*   **`-q` / `--quiet` (`ERROR` level):** Silences all progress/status outputs, printing only critical errors to `sys.stderr` (and output files if `-o` is used).
*   **`-v` / `--verbose` (`DEBUG` level):** Shows detailed debugging information (e.g. detailed pixel grid resampling details, API downloads) and prints full Python tracebacks in case of exceptions.

> [!WARNING]
> The `--verbose` and `--quiet` options are mutually exclusive. Passing both flags simultaneously will raise a `UsageError`.

### Streams Separation
To facilitate scripting and piping:
*   **`stdout`:** Receives only programmatic outputs (like the final pretty-printed dosimetry summary table or JSON results if redirection is not used).
*   **`stderr`:** Receives all logging entries formatted as `[LEVEL] radot: message` (e.g. `[INFO]`, `[DEBUG]`, `[WARNING]`).

