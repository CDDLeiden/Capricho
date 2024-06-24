"""Configure luguru logger for the project."""

import sys

from loguru import logger


def setup_logger(level="INFO"):
    """Initialize the logger with a default log level and format"""
    colorful_format = (
        "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | "
        "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - "
        "<level>{message}</level>"
    )
    logger.remove()
    logger.add(sys.stderr, format=colorful_format, level=level, colorize=True)


def set_log_level(level):
    """Update the logger's level."""
    setup_logger(level=level)


setup_logger()  # `logger` is configured whenever it gets imported from this module
