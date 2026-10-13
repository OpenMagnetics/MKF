#pragma once
// The ONLY way MKF includes Eigen (ABT #1850). Every Eigen symbol MKF instantiates is private to the binary MKF is
// linked into; never include <Eigen/...> directly in MKF.
//
// Why. Eigen derives how it allocates a dynamic matrix from the instruction set of the translation unit: built for
// AVX / AVX-512 (MKF's -march=native) it over-allocates and keeps the malloc pointer 8 bytes below the block, released
// by free(block[-1]); built without, plain malloc/free. Eigen 3.4 offers no setting that decouples the two
// (EIGEN_MAX_ALIGN_BYTES=16 leaves EIGEN_DEFAULT_ALIGN_BYTES, which the allocator uses, at 64). Eigen's functions are
// templates and inline functions: weak symbols with default visibility, of which the dynamic linker keeps ONE definition
// per process, the first in load order. So a program that uses Eigen itself, built with other flags than MKF, supplies
// MKF's out-of-line Eigen functions while MKF keeps its inlined ones: OMFEM (default flags) loading MKF through
// libmvb++.so had MKF's parallel current split allocate its strand system with OMFEM's plain-malloc
// Eigen::internal::aligned_malloc and release it with MKF's inlined free(block[-1]): free(0x30), SIGSEGV.
//
// With hidden visibility MKF's Eigen instantiations are no dynamic symbols: nothing outside can replace them and MKF
// replaces nobody's. -fvisibility-inlines-hidden would not do: Eigen's out-of-class member templates
// (PartialPivLU::compute, SimplicialCholeskyBase::factorize_preordered, ...) are not inline and stay exported.
//
// The pragma also hides whatever ELSE is first declared inside it, so every system header Eigen includes is included
// here first (a C function first declared hidden does not link: "hidden symbol `__assert_fail' isn't defined").
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cerrno>
#include <cfloat>
#include <climits>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <functional>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <new>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif
#if defined(_OPENMP)
#include <omp.h>
#endif

// An MKF type built from an Eigen one (a functor deriving from Eigen::DenseFunctor, ...) is declared MKF_EIGEN_PRIVATE:
// a type cannot be more visible than its base.
#if defined(__GNUC__) && !defined(_WIN32)
#define MKF_EIGEN_PRIVATE __attribute__((visibility("hidden")))
#pragma GCC visibility push(hidden)
#else
#define MKF_EIGEN_PRIVATE
#endif
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>
#if defined(__GNUC__) && !defined(_WIN32)
#pragma GCC visibility pop
#endif
