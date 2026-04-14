---
name: Build with ninja -j4
description: User prefers building with ninja -j4 directly in build dir instead of cmake --build
type: feedback
---

Use `cd build && ninja -j4` to compile, not `cmake --build build`.

**Why:** User explicitly corrected to use ninja -j4. Likely to avoid OOM from too many parallel jobs and for direct control.

**How to apply:** Always use `ninja -j4` in the build directory for compilation.
