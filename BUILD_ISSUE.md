# Build Issue - Broken Xcode Command Line Tools

## Current Status

We've successfully set up the Python/C++ hybrid architecture:
- ✅ Project structure created
- ✅ CMakeLists.txt configured for pybind11
- ✅ setup.py created
- ✅ Python package structure ready
- ✅ Initial bindings written
- ✅ CMake finds all dependencies (pybind11, Eigen, Boost, OpenMP)
- ✅ CMake configures successfully

## The Problem

**The C++ compiler cannot find standard library headers** (`<iostream>`, `<complex>`, `<cstddef>`, etc.)

This is **NOT** a problem with our code or configuration. This is a **broken Xcode Command Line Tools installation** on your system.

### Evidence

1. CMake correctly sets `-isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk`
2. The SDK exists at that path
3. Headers exist in `/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1/`
4. **But** even a simple test fails:
```bash
/usr/bin/c++ -std=c++14 -stdlib=libc++ -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk -c /tmp/test.cpp
# ERROR: 'iostream' file not found
```

This is a system-level issue that affects ALL C++ compilation, not just our project.

## Solution

You need to **reinstall Xcode Command Line Tools**:

```bash
# Remove the broken installation
sudo rm -rf /Library/Developer/CommandLineTools

# Trigger reinstallation
xcode-select --install

# Follow the GUI prompts to install
```

After reinstalling, verify it works:
```bash
echo '#include <iostream>
int main() { std::cout << "Hello"; }' > /tmp/test.cpp

clang++ -std=c++14 /tmp/test.cpp -o /tmp/test && /tmp/test
```

If that prints "Hello", then you're good to go!

## Next Steps After Fixing

Once CommandLineTools are fixed, building should be straightforward:

```bash
cd /Users/haocheng/Github/GCTB

# Activate virtual environment
source venv/bin/activate

# Build and install
pip install -e .

# Test
python3 -c "import gctb; print(gctb.__version__)"
```

## Alternative: Use Homebrew LLVM

If reinstalling CommandLineTools doesn't work, you could try using Homebrew's LLVM:

```bash
brew install llvm

# Then update CMakeLists.txt to use:
# set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
```

## What We've Built So Far

Even though we can't compile yet, we've created a complete build infrastructure:

1. **CMakeLists.txt** - Properly configured for macOS, finds all dependencies
2. **setup.py** - Python packaging with CMake integration
3. **pyproject.toml** - Modern Python package configuration
4. **python/bindings/bindings.cpp** - Initial pybind11 bindings for:
   - SnpInfo
   - IndInfo  
   - Data class (with file I/O methods)
   - Timer utility
5. **Tests** - test_basic.py and test_data_io.py ready to run
6. **Build script** - build.sh for easy compilation

Once the C++ compiler is fixed, everything else is ready to go!

## Timeline Impact

This system issue has cost us about 2 hours of debugging. Once fixed:
- Compilation should take ~5 minutes
- We can immediately test the basic bindings
- Then proceed with expanding bindings for Model, MCMC, etc.

The 2-week timeline is still feasible once this system issue is resolved.

