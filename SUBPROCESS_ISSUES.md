# Subprocess Issues with Python C++ Extensions

**Critical Information for Testing and Development**

---

## ⚠️ **The Problem**

When testing GCTB Python package, using `subprocess.run()` or `subprocess.Popen()` to launch Python scripts that import the C++ extension **causes crashes**.

### Symptoms:

```python
import subprocess
result = subprocess.run(['python3', '-c', 'import gctb'])
# Result: Exit code -5 (SIGTRAP)
# "Python quit unexpectedly"
```

**Error:**
- Exit code: -5, -6, or 133
- Signal: SIGTRAP or SIGABRT
- Message: "Python quit unexpectedly" or no error message
- All subprocess attempts fail

---

## 🔍 **Root Cause**

### Why This Happens:

The GCTB Python package uses:
1. **pybind11** - C++ to Python bindings
2. **C++ shared library** (`_core.so`) with:
   - OpenMP (libomp.dylib)
   - Eigen (C++ template library)
   - Complex memory management
3. **Native libraries** that don't play well with subprocess

### Technical Details:

#### 1. **Library Initialization Order**

When using subprocess:
```
subprocess.run() 
  → fork() 
    → Python interpreter 
      → import gctb 
        → load _core.so 
          → libomp.dylib initialization
            → ❌ CRASH
```

Something in this chain fails, likely:
- OpenMP thread initialization after fork()
- Memory mapping issues
- Signal handler conflicts

#### 2. **Signal Handling Conflicts**

- C++ libraries (especially OpenMP) set signal handlers
- subprocess interferes with signal handling
- Results in SIGTRAP (signal 5) or SIGABRT (signal 6)

#### 3. **macOS ARM64 Specific Issues**

- M1/M2 processors have additional constraints
- Rosetta translation layer can cause issues
- macOS sandbox restrictions on forked processes

#### 4. **OpenMP + fork() = Problems**

OpenMP is not fork-safe:
```c++
// Inside C++ code
#pragma omp parallel for
// If parent process forks after OpenMP init, child may crash
```

---

## ❌ **What DOESN'T Work**

### 1. **Python subprocess.run()**

```python
# ❌ FAILS - Will crash
import subprocess
subprocess.run(['python3', '-m', 'gctb.cli', 'bayes', ...])
# Exit code: -5
```

### 2. **Python subprocess.Popen()**

```python
# ❌ FAILS - Will crash
import subprocess
p = subprocess.Popen(['python3', '-c', 'import gctb'])
p.wait()  # Returns -5
```

### 3. **Changing OMP_NUM_THREADS Mid-Process**

```python
# ❌ DOESN'T WORK - Too late!
import os
import gctb

os.environ['OMP_NUM_THREADS'] = '4'  # Ignored!
# OpenMP already initialized when gctb imported
```

### 4. **multiprocessing.Pool with C++ Extension**

```python
# ❌ FAILS - Same fork() issues
from multiprocessing import Pool
import gctb

def work(x):
    return gctb.run_analysis(x)

pool = Pool(4)  # Will crash
pool.map(work, data)
```

---

## ✅ **What WORKS**

### 1. **Shell Scripts (Recommended)**

```bash
#!/bin/bash
# ✅ WORKS PERFECTLY

# Each run in fresh process
OMP_NUM_THREADS=1 python3 -m gctb.cli bayes --bfile data --pheno pheno.txt
OMP_NUM_THREADS=2 python3 -m gctb.cli bayes --bfile data --pheno pheno.txt
OMP_NUM_THREADS=4 python3 -m gctb.cli bayes --bfile data --pheno pheno.txt
```

**Why it works:**
- Fresh Python interpreter each time
- No fork() involved
- Clean library initialization
- Environment set BEFORE Python starts

### 2. **bash -c for Programmatic Testing**

```bash
# ✅ WORKS
OMP_NUM_THREADS=1 bash -c 'python3 -m gctb.cli bayes ...'
OMP_NUM_THREADS=2 bash -c 'python3 -m gctb.cli bayes ...'
```

**Why it works:**
- `bash -c` creates fresh shell
- No Python subprocess
- Proper process isolation

### 3. **Direct Python API (Single Process)**

```python
# ✅ WORKS - If you don't need to change threads
import os
os.environ['OMP_NUM_THREADS'] = '4'  # MUST be before import!

import gctb

# Now use the API directly
data = gctb.Data()
# ... run analysis ...
```

**Important:** Environment MUST be set BEFORE importing gctb!

### 4. **System Calls via os.system()**

```python
# ✅ WORKS (but not recommended, use shell script instead)
import os
os.system('OMP_NUM_THREADS=4 python3 -m gctb.cli bayes ...')
```

---

## 📋 **Best Practices**

### For Testing C++ Extensions:

#### ✅ **DO:**

1. **Use Shell Scripts:**
```bash
#!/bin/bash
for threads in 1 2 4 8; do
    OMP_NUM_THREADS=$threads python3 script.py
done
```

2. **Set Environment Before Import:**
```python
import os
os.environ['OMP_NUM_THREADS'] = '4'
# NOW safe to import
import gctb
```

3. **Use Fresh Processes:**
```bash
# Each test in new process
bash -c 'OMP_NUM_THREADS=1 python3 test.py'
bash -c 'OMP_NUM_THREADS=2 python3 test.py'
```

#### ❌ **DON'T:**

1. **Don't Use subprocess with C++ Extensions:**
```python
# ❌ NO!
subprocess.run(['python3', '-c', 'import gctb'])
```

2. **Don't Change OpenMP Settings After Import:**
```python
# ❌ NO!
import gctb
os.environ['OMP_NUM_THREADS'] = '4'  # Too late!
```

3. **Don't Use multiprocessing.Pool:**
```python
# ❌ NO!
from multiprocessing import Pool
pool.map(gctb_function, data)  # Will crash
```

---

## 🔧 **Debugging Tips**

### If You Get Subprocess Crashes:

1. **Check Exit Code:**
```python
result = subprocess.run([...], capture_output=True)
print(f"Exit code: {result.returncode}")
# -5 = SIGTRAP (library initialization issue)
# -6 = SIGABRT (assertion failure)
# 133 = Threading warning (less severe)
```

2. **Test Direct Import:**
```python
# If this crashes, problem is in the library itself
python3 -c "import gctb; print('OK')"
```

3. **Test with Minimal Case:**
```python
# Simplest possible case
import sys
sys.path.insert(0, 'python')
import gctb
print("Imported successfully")
```

4. **Check Library Loading:**
```bash
# macOS - see what libraries are loaded
otool -L python/gctb/_core*.so

# Should show:
#   /opt/homebrew/opt/libomp/lib/libomp.dylib
#   /usr/lib/libc++.1.dylib
```

---

## 📚 **Technical Background**

### Why pybind11 + subprocess Is Problematic:

1. **Shared Library State:**
   - C++ libraries maintain global state
   - fork() duplicates this state
   - Child process inherits parent's library state
   - Can cause conflicts, especially with threads

2. **Thread-Local Storage:**
   - OpenMP uses thread-local storage
   - fork() only copies calling thread
   - Other threads' TLS is lost
   - Results in crashes when accessed

3. **Signal Handlers:**
   - C++ libraries register signal handlers
   - subprocess may interfere with these
   - Conflicting handlers → SIGTRAP

4. **Memory Management:**
   - C++ uses complex memory allocators
   - fork() duplicates memory state
   - Allocators may be in inconsistent state
   - Can cause crashes on first allocation

### macOS Specific Issues:

1. **Code Signing:**
   - macOS requires signed binaries for some operations
   - subprocess may trigger signing checks
   - Can fail with unsigned/development code

2. **System Integrity Protection (SIP):**
   - Restricts library injection
   - Can interfere with dynamic loading
   - May cause crashes in subprocess

3. **ARM64 Architecture:**
   - M1/M2 have stricter memory ordering
   - Some x86 assumptions don't hold
   - Can expose latent bugs in fork() paths

---

## 🎯 **Real-World Examples**

### Example 1: OpenMP Scaling Test (What We Did)

**❌ First Attempt (Failed):**
```python
# test_openmp_scaling.py
for num_threads in [1, 2, 4, 8]:
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = str(num_threads)
    
    result = subprocess.run(
        ['python3', '-m', 'gctb.cli', 'bayes', ...],
        env=env
    )
    # Result: All crashed with exit code -5
```

**✅ Solution (Worked):**
```bash
# test_openmp_comprehensive.sh
OMP_NUM_THREADS=1 bash -c 'python3 -m gctb.cli bayes ...'
OMP_NUM_THREADS=2 bash -c 'python3 -m gctb.cli bayes ...'
OMP_NUM_THREADS=4 bash -c 'python3 -m gctb.cli bayes ...'
OMP_NUM_THREADS=8 bash -c 'python3 -m gctb.cli bayes ...'
# Result: All worked perfectly
```

### Example 2: Benchmarking (Similar Issue)

**❌ Don't Do This:**
```python
import subprocess
import time

times = []
for threads in [1, 2, 4]:
    start = time.time()
    subprocess.run(['python3', 'run_analysis.py', str(threads)])
    times.append(time.time() - start)
# Crashes
```

**✅ Do This:**
```bash
#!/bin/bash
for threads in 1 2 4; do
    echo "Testing $threads threads"
    time OMP_NUM_THREADS=$threads python3 run_analysis.py
done
```

---

## 🐛 **Common Error Messages**

### 1. "Python quit unexpectedly"
```
Exit code: -5
Signal: SIGTRAP
```
**Cause:** subprocess + C++ extension  
**Solution:** Use shell script

### 2. "resource_tracker: leaked semaphore"
```
Exit code: 133
Warning: leaked semaphore objects
```
**Cause:** Threading cleanup issue  
**Impact:** Cosmetic only, analysis works  
**Solution:** Ignore or suppress warning

### 3. "Abort trap: 6"
```
Exit code: -6  
Signal: SIGABRT
```
**Cause:** Assertion failure in C++ code  
**Solution:** Check if OpenMP properly initialized

---

## 📖 **Further Reading**

### Python Documentation:
- subprocess module: https://docs.python.org/3/library/subprocess.html
- multiprocessing notes on fork: https://docs.python.org/3/library/multiprocessing.html#contexts-and-start-methods

### OpenMP Documentation:
- Thread affinity: https://www.openmp.org/spec-html/5.0/openmpse52.html
- Fork safety: OpenMP is NOT fork-safe by design

### pybind11 Documentation:
- FAQ on threading: https://pybind11.readthedocs.io/en/stable/advanced/misc.html#module-destructors

---

## ✅ **Summary**

### The Rule of Thumb:

**If your Python package has C++ extensions (especially with OpenMP or threads):**

✅ **DO:**
- Use shell scripts for testing
- Set environment before import
- Use fresh processes for each test

❌ **DON'T:**
- Use subprocess.run() or subprocess.Popen()
- Use multiprocessing.Pool
- Change thread settings after import

### Quick Reference:

```bash
# ✅ Testing pattern that always works:
#!/bin/bash
for setting in option1 option2 option3; do
    SETTING=$setting python3 -m package.cli command
done
```

```python
# ✅ Single-process pattern:
import os
os.environ['SETTING'] = 'value'  # BEFORE import!
import package
# Use package API
```

---

## 🎓 **Lessons Learned**

1. **subprocess + C++ extensions = Bad combination**
2. **OpenMP must be configured before library loads**
3. **Shell scripts are more reliable than Python subprocess**
4. **Fresh processes avoid state conflicts**
5. **Always test import in subprocess before relying on it**

---

**This document applies to any Python package using:**
- pybind11
- Cython
- ctypes with complex C++ libraries
- OpenMP or threading libraries
- Native extensions with global state

**When in doubt: Use a shell script!** 🎯

---

**Last Updated:** 2025-11-07  
**Platform:** macOS ARM64 (M1/M2)  
**Python Version:** 3.13  
**Package:** GCTB v1.0.0

