import logging
import sys
from typing import Union

# Root package logger
logger = logging.getLogger("radot")
logger.setLevel(logging.INFO)

class DynamicStderrStream:
    def write(self, data):
        sys.stderr.write(data)
    def flush(self):
        sys.stderr.flush()

# Default console stream handler pointing to sys.stderr dynamically
_default_handler = logging.StreamHandler(DynamicStderrStream())
_default_handler.setFormatter(
    logging.Formatter("[%(levelname)s] %(name)s: %(message)s")
)
logger.addHandler(_default_handler)


def set_log_level(level: Union[int, str]) -> None:
    """
    Set the log level of the radot library logger.
    Supports standard logging levels, as well as the 'QUIET' string helper.
    """
    if isinstance(level, str):
        clean_level = level.upper().strip()
        if clean_level == "QUIET":
            target_level = logging.ERROR
        else:
            # Map standard logger strings
            target_level = getattr(logging, clean_level, logging.INFO)
    else:
        target_level = level

    logger.setLevel(target_level)
    for handler in logger.handlers:
        handler.setLevel(target_level)


def configure_logging(verbose: bool = False, quiet: bool = False) -> None:
    """
    Configure the logging handlers and level for CLI or scripts.
    - verbose: set level to DEBUG
    - quiet: set level to ERROR
    """
    if verbose:
        set_log_level(logging.DEBUG)
    elif quiet:
        set_log_level(logging.ERROR)
    else:
        set_log_level(logging.INFO)
