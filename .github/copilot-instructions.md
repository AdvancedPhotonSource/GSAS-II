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
python -m pip install pytest nox
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

# Using nox for testing workflow
python -m nox -s tests
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

### Linting and Code Quality
```bash
# Using nox (recommended)
python -m nox -s pylint    # Static code analysis
```

### Documentation
```bash
# Build documentation  
python -m nox -s docs

# Build and serve docs locally
python -m nox -s docs -- --serve
```

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