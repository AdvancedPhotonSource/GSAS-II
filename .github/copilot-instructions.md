# GSAS-II - Crystallography Software Package

GSAS-II is a comprehensive Python package for analysis of x-ray and neutron diffraction data, including single-crystal, powder, and time-of-flight data. It includes both GUI (wxPython) and scriptable interfaces, with compiled Fortran extensions for performance-critical computations.

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Working Effectively

### Bootstrap and Build the Repository
```bash
# Install core build dependencies
python -m pip install --upgrade pip setuptools wheel
python -m pip install meson-python ninja numpy cython scipy pycifrw

# Build the package (takes ~16 seconds, NEVER CANCEL, set timeout to 30+ minutes)
python -m build -wnx --no-isolation

# OR install in editable mode for development (takes ~12 seconds, NEVER CANCEL, set timeout to 30+ minutes)  
python -m pip install -e . --no-build-isolation
```

### Install Optional Dependencies
```bash
# For GUI functionality
python -m pip install wxpython matplotlib pyopengl

# For additional features
python -m pip install pillow h5py imageio requests gitpython pybaselines

# For testing and development
python -m pip install pytest
```

### Run Tests
```bash
# Add compiled executables to PATH before running tests
export PATH="$PATH:$(pwd)/build/cp312/sources"

# Alternative: Copy executables to a location in PATH (one-time setup)
# mkdir -p ~/.local/bin
# cp build/cp312/sources/{LATTIC,convcell} ~/.local/bin/

# Run all tests (some require network connectivity, ~3-4 seconds for working tests)
python -m pytest tests/ -v

# Run specific test modules that work offline
python -m pytest tests/test_lattice.py tests/test_nistlat.py tests/test_elm.py -v
```

### Validation Scenarios
Always run these validation steps after making changes:

1. **Basic Import Test**:
   ```bash
   python -c "import GSASII; print('GSAS-II import successful')"
   ```

2. **Scriptable Interface Test**:
   ```bash
   python -c "
   import GSASII.GSASIIscriptable as G2sc
   gpx = G2sc.G2Project(newgpx='/tmp/test.gpx')
   print('GSASIIscriptable working correctly')
   "
   ```

3. **Binary Extensions Test**:
   ```bash
   python -c "
   import GSASII.GSASIIlattice as G2lat
   import GSASII.GSASIIspc as G2spc
   print('Compiled extensions working')
   "
   ```

4. **Run Core Functionality Tests**:
   ```bash
   export PATH="$PATH:$(pwd)/build/cp312/sources"
   python -m pytest tests/test_lattice.py::test_Brav -v
   ```

## Build System Details

### Timing Expectations
- **NEVER CANCEL**: Build takes ~16 seconds on typical hardware but may take up to 30 minutes on slower systems
- **NEVER CANCEL**: Editable install takes ~12 seconds but may take up to 30 minutes  
- **NEVER CANCEL**: Test suite takes 3-4 seconds for working tests but up to 10 minutes with all tests
- Set explicit timeouts of 30+ minutes for build commands and 15+ minutes for test commands

### Build Architecture
- Uses **meson** build system with **f2py** for Fortran compilation
- Compiles Fortran extensions: pyspg, pydiffax, pypowder, pytexture, pack_f, histogram2d, unpack_cbf
- Builds standalone executables: LATTIC, convcell
- Requires gfortran compiler

### Build Troubleshooting
- If binaries are missing: run `python -m pip install -e . --no-build-isolation` 
- For PATH issues with executables: `export PATH="$PATH:$(pwd)/build/cp312/sources"`
- Build output location: `./build/cp312/sources/` (for Python 3.12)

## Running Applications

### Scriptable Interface (Recommended for Development)
```bash
python -c "
import GSASII.GSASIIscriptable as G2sc
gpx = G2sc.G2Project(newgpx='project.gpx')
# Add your GSAS-II scripting code here
"
```

### GUI Application
```bash
# Note: GUI requires display/X11 - will not work in headless environments
python -m GSASII
```

## Development Workflow

### Common Development Tasks
Always validate changes by running:
1. Build the package: `python -m build -wnx --no-isolation` 
2. Install in editable mode: `python -m pip install -e . --no-build-isolation`
3. Add executables to PATH: `export PATH="$PATH:$(pwd)/build/cp312/sources"`
4. Run working tests: `python -m pytest tests/test_lattice.py tests/test_elm.py -v`
5. Test core functionality: Run the validation scenarios above

## Troubleshooting

### Common Issues
- **"binary load error: pyspg not found"**: Run editable install to compile binaries
- **"FileNotFoundError: 'LATTIC'"**: Add `./build/cp312/sources` to PATH  
- **Network connectivity test failures**: Normal in restricted environments, focus on offline tests
- **wxPython GUI issues**: GUI requires display, use scriptable interface in headless environments

### Dependencies Not Found
Install the specific missing dependency:
```bash
# For seekpath k-vector functionality
python -m pip install seekpath

# For complete GUI experience  
python -m pip install wxpython matplotlib pyopengl pillow h5py
```

## Key Projects in Codebase

### Main Modules
- **GSASII/**: Main package directory
  - **GSASIIscriptable.py**: Primary API for programmatic use
  - **GSASIIGUI.py**: Main GUI application entry point
  - **GSASIIdata.py**: Core data structures and file I/O
  - **GSASIImath.py**: Mathematical routines for crystallography
  - **GSASIIlattice.py**: Lattice parameter and space group operations
  - **GSASIIspc.py**: Space group and symmetry operations

### Build Files  
- **meson.build**: Main build configuration
- **sources/meson.build**: Fortran extension build configuration
- **pyproject.toml**: Package metadata and dependencies

### Testing
- **tests/**: Self-contained test suite
  - **test_lattice.py**: Lattice and crystallography tests (always work offline)
  - **test_nistlat.py**: NIST lattice utility tests (require executables in PATH)
  - **test_scriptref.py**: Scripting interface integration tests

### Documentation  
- **docs/**: Sphinx documentation source
- **README.md**: Project overview and links

Always check binary compilation status and add executables to PATH before running tests or using advanced functionality.

## Code Style and Quality

### Linting and Formatting
The project uses **ruff** for linting and code quality checks:
```bash
# The project is configured with ruff in pyproject.toml
# Key conventions:
# - Python 3.10+ syntax
# - NumPy-style docstrings preferred
# - Imports organized with isort conventions
# - Many flake8-style checks enabled (see pyproject.toml [tool.ruff.lint])
```

### Common Style Guidelines
- **Avoid excessive changes**: Don't reformat code unnecessarily
- **Follow existing patterns**: Match the style of surrounding code
- **Legacy code**: Much of the codebase predates modern Python conventions; gradual improvements are acceptable
- **Type hints**: Not universally used; add them to new code when beneficial
- **Documentation**: Update inline comments when making complex changes

### Configuration Files
- `pyproject.toml`: Contains ruff, pytest, and mypy configuration
- `meson.build`: Build system configuration for Fortran extensions

## Architecture Overview

### Code Organization
GSAS-II has a dual-interface architecture:

1. **GUI Interface** (`GSASIIGUI.py`, `*GUI.py` files):
   - wxPython-based graphical interface
   - Files ending in `GUI.py` contain GUI-specific code
   - Requires display/X11; won't work in headless environments

2. **Scriptable Interface** (`GSASIIscriptable.py`):
   - Python API for programmatic use
   - Preferred for automated workflows and testing
   - Works in headless environments
   - Entry point: `G2Project` class

3. **Core Modules** (no `GUI` suffix):
   - `GSASIIdata.py`: Data structures and file I/O
   - `GSASIImath.py`: Crystallography calculations
   - `GSASIIlattice.py`, `GSASIIspc.py`: Space group operations
   - Can be used independently of GUI

### Fortran Extensions
Performance-critical code is in Fortran (in `sources/`):
- Compiled to Python extensions via f2py
- Binary modules: pyspg, pydiffax, pypowder, pytexture, pack_f, histogram2d
- Standalone executables: LATTIC, convcell
- Require rebuilding after changes to Fortran source

## CI/CD and Testing

### Continuous Integration
Workflows in `.github/workflows/`:
- `smoke_test.yml`: Builds and tests on Ubuntu, macOS, Windows with Python 3.10-3.13
- `bilbao-test.yml`: Tests Bilbao symmetry server integration
- Uses pixi for environment management
- Tests run on push and pull requests

### Testing Strategy
- **Fast tests**: `test_lattice.py`, `test_elm.py` (work offline, ~3-4 seconds)
- **Network tests**: `test_kvec.py`, `test_bilbao.py` (require internet)
- **Integration tests**: `test_scriptref.py` (test scriptable interface)
- Run tests early and often during development

### Test Execution
```bash
# After building, always set PATH
export PATH="$PATH:$(pwd)/build/cp312/sources"

# Run fast offline tests
python -m pytest tests/test_lattice.py tests/test_elm.py -v

# Run all tests (some may fail without network)
python -m pytest tests/ -v
```

## Contributing Best Practices

### Before Making Changes
1. Understand the dual architecture (GUI vs scriptable)
2. Check if changes affect Fortran code (requires rebuild)
3. Identify which tests validate your changes
4. Check if binaries/executables are in PATH

### Making Changes
1. Keep changes minimal and focused
2. Test with scriptable interface when possible (faster iteration)
3. Rebuild only when necessary (Fortran changes, dependency updates)
4. Run relevant tests, not entire suite
5. Update documentation if changing public APIs

### Debugging
- For binary issues: Check build output, verify PATH
- For import errors: Ensure editable install completed successfully
- For GUI issues: Try scriptable interface first to isolate problem
- For test failures: Check if executables are in PATH

## Common Patterns and Conventions

### Working with Space Groups
```python
import GSASII.GSASIIspc as G2spc
# Space groups are central to crystallography
# Use GSASIIspc for symmetry operations
```

### Working with Lattice Parameters
```python
import GSASII.GSASIIlattice as G2lat
# Lattice calculations and transformations
```

### Scriptable Interface Pattern
```python
import GSASII.GSASIIscriptable as G2sc
# Always start with G2Project
gpx = G2sc.G2Project(newgpx='/path/to/project.gpx')
# Add histograms, phases, do refinements
```

### File I/O
```python
import GSASII.GSASIIdata as G2data
# Standardized data structures and file operations
```

## Gotchas and Important Notes

### Build System
- **Never cancel builds prematurely**: Meson builds can take 12-30 seconds, but may need up to 30 minutes on slow systems
- **Build artifacts location**: `./build/cp3XX/sources/` where XX is Python version
- **Missing executables**: Most common issue; solution is `export PATH`

### Dependencies
- **NumPy**: Required; build system needs it at build time
- **wxPython**: Optional; only for GUI
- **Fortran compiler**: Required; gfortran recommended
- **pixi**: Used in CI/CD; optional for local development

### Platform-Specific Issues
- **Windows**: May need Visual Studio Build Tools for compilation
- **macOS**: XCode command line tools required
- **Linux**: Usually works out-of-box with build-essential

### Testing Environment
- Some tests require network connectivity (Bilbao server, k-vector tests)
- GUI tests require display; use scriptable interface in headless environments
- Executables must be in PATH for certain tests (LATTIC, convcell)