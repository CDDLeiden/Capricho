"""Configure loguru logger for the project."""

import sys

from loguru import logger

# Compact format for INFO and above: level + message only
COMPACT_FORMAT = "<level>{level: <8}</level> | <level>{message}</level>"

# Verbose format for DEBUG and below: timestamp + source location + message
VERBOSE_FORMAT = (
    "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | "
    "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
    "<level>{message}</level>"
)


def setup_logger(level="INFO", out_file=None, _sink=None, _verbose_sink=None):
    """Initialize the logger with a two-tier format.

    INFO and above use a compact format (no timestamp, no source location).
    DEBUG and below use a verbose format with timestamps and source locations.

    Args:
        level: Minimum log level to display.
        out_file: Optional file path for log output.
        _sink: Override sink for testing (replaces stderr for compact format).
        _verbose_sink: Override sink for testing (replaces stderr for verbose format).
    """
    logger.remove()

    compact_sink = _sink if _sink is not None else sys.stderr
    colorize = _sink is None  # Only colorize when writing to stderr

    level_upper = level.upper() if isinstance(level, str) else level

    # Always add compact format for INFO+
    logger.add(
        compact_sink,
        format=COMPACT_FORMAT,
        level="INFO",
        colorize=colorize,
        filter=lambda record: record["level"].no >= logger.level("INFO").no,
    )

    # Add verbose format for DEBUG/TRACE when level is below INFO
    if level_upper in ("DEBUG", "TRACE"):
        verbose_sink = _verbose_sink if _verbose_sink is not None else compact_sink
        verbose_colorize = _verbose_sink is None and _sink is None
        logger.add(
            verbose_sink,
            format=VERBOSE_FORMAT,
            level=level_upper,
            colorize=verbose_colorize,
            filter=lambda record: record["level"].no < logger.level("INFO").no,
        )

    if out_file:
        logger.add(out_file, format=VERBOSE_FORMAT, level=level_upper)


def set_log_level(level):
    """Update the logger's level."""
    setup_logger(level=level)


setup_logger()  # `logger` is configured whenever it gets imported from this module
