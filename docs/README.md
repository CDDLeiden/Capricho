# CAPRICHO Documentation

This directory contains the Sphinx documentation for CAPRICHO.

## Building the Documentation

Install CAPRICHO with documentation dependencies:

```bash
pip install -e ".[dev]"
```

Build the HTML documentation:

```bash
cd docs
make html
```

The built documentation will be in `build/html/`. Open `build/html/index.html` in a browser to view.

## Structure

- `source/` - Documentation source files
- `build/` - Built documentation (created when building)
- `Makefile` - Build commands for Unix/Linux/macOS
- `make.bat` - Build commands for Windows