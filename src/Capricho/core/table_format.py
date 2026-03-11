"""Format pandas DataFrames for console display in markdown or CSV. (optional per-row ANSI color cycling)"""

import csv
from io import StringIO
from itertools import cycle

import pandas as pd

RESET = "\033[0m"
COLORS = [
    "\033[97m",  # bright white
    "\033[94m",  # bright blue
    "\033[91m",  # bright red
    "\033[92m",  # bright green
    "\033[93m",  # bright yellow
    "\033[96m",  # bright cyan
]


def _is_markdown_separator(line: str) -> bool:
    """Check if a line is a markdown table separator (e.g. |:---|:---|)."""
    stripped = line.strip()
    return bool(stripped.startswith("|") and "---" in stripped)


def _colorize_pipe_line(line: str) -> str:
    """Colorize cells in a pipe-delimited line."""
    parts = line.split("|")
    # parts[0] is before first |, parts[-1] is after last |
    # parts[1:-1] are the actual cells
    if len(parts) < 3:
        return line
    colored = [parts[0]]
    color_iter = cycle(COLORS)
    for cell in parts[1:-1]:
        color = next(color_iter)
        colored.append(f"{color}{cell}{RESET}")
    colored.append(parts[-1])
    return "|".join(colored)


def _colorize_columns_pipe(text: str) -> str:
    """Apply per-column ANSI colors for markdown pipe tables."""
    lines = text.split("\n")
    result = []
    for line in lines:
        if not line.strip() or _is_markdown_separator(line):
            result.append(line)
        else:
            result.append(_colorize_pipe_line(line))
    return "\n".join(result)


def _colorize_columns_csv(text: str) -> str:
    """Apply per-column ANSI colors for CSV format."""
    reader = csv.reader(StringIO(text))
    result_lines = []
    for row in reader:
        color_iter = cycle(COLORS)
        colored_cells = []
        for cell in row:
            color = next(color_iter)
            colored_cells.append(f"{color}{cell}{RESET}")
        result_lines.append(",".join(colored_cells))
    return "\n".join(result_lines)


def format_dataframe(df: pd.DataFrame, fmt: str, colorize: bool = False) -> str:
    """Format a DataFrame for console display.

    Args:
        df: DataFrame to format.
        fmt: Output format - "markdown" or "csv".
        colorize: If True, apply ANSI color cycling to columns.

    Returns:
        Formatted string representation of the DataFrame.
    """
    if fmt == "csv":
        text = df.to_csv(index=False)
    elif fmt == "markdown":
        text = df.to_markdown(index=False, tablefmt="pipe")
        if text is None:
            text = ""
    else:
        raise ValueError(f"Unsupported format: '{fmt}'. Use 'markdown' or 'csv'.")

    if colorize and text.strip():
        if fmt == "csv":
            text = _colorize_columns_csv(text)
        else:
            text = _colorize_columns_pipe(text)

    return text
