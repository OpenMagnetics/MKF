#pragma once
#include <array>
#include <complex>
#include <cstddef>
#include <functional>
#include <vector>

namespace OpenMagnetics {

class BlockEntries;  // the kernel entries a block evaluated (HierarchicalMatrix.cpp)

// A real matrix whose entries come from a kernel between two point sets (rows and columns), stored as a hierarchical
// matrix: the rows and the columns are clustered by position, a pair of clusters far apart against their sizes is
// smooth and stored in low rank (adaptive cross approximation, to `tolerance` relative to the block), and every other
// pair is stored dense. A low-rank block is checked against the kernel on sampled entries and kept dense when it misses
// them. Entries the kernel cannot express smoothly (a kernel that changes form between some row and column pairs) are
// added as sparse corrections, applied with the product. A kernel that changes form with where a point lies (a toroid's
// bore and the outside of its ring, a multi-window core's windows) names that form per point (PointForm): a cluster that
// holds several forms is approximated with no other, only split. The cross approximation follows its pivots and never
// sees a part of a block they cannot reach, and the sampled entries need not land in it: two forms in one block were
// accepted up to 95 % wrong (ABT #2498).
class HierarchicalMatrix {
  public:
    using Entry = std::function<double(size_t row, size_t column)>;

    // Several entries per (row, column) at once, written to values[0..count).
    using MultiEntry = std::function<void(size_t row, size_t column, double* values)>;

    // The kernel's form at a point; none: one form everywhere.
    using PointForm = std::function<size_t(const std::array<double, 2>& point)>;

    HierarchicalMatrix() = default;
    HierarchicalMatrix(const std::vector<std::array<double, 2>>& rowPoints, const std::vector<std::array<double, 2>>& columnPoints, const Entry& entry,
                       double tolerance, const PointForm& form = {});

    // `count` matrices over the same points from a kernel that gives all their entries at once (the two components of a
    // field): they share the clustering and the block partition, and a dense block evaluates the kernel once for all.
    static std::vector<HierarchicalMatrix> build_set(const std::vector<std::array<double, 2>>& rowPoints, const std::vector<std::array<double, 2>>& columnPoints,
                                                     const MultiEntry& entry, size_t count, double tolerance, const PointForm& form = {});

    size_t rows() const { return _rows; }
    size_t columns() const { return _columns; }
    // y += A x.
    void multiply_add(const std::vector<std::complex<double>>& x, std::vector<std::complex<double>>& y) const;
    // y[v] += A x[v] for several vectors, the matrix read once for all of them; each result is the one multiply_add gives.
    void multiply_add_batch(const std::vector<const std::vector<std::complex<double>>*>& x, const std::vector<std::vector<std::complex<double>>*>& y) const;
    void add_sparse(size_t row, size_t column, double value);
    // Values stored (dense entries, low-rank factors and sparse corrections): the memory the matrix holds.
    size_t stored_values() const;

  private:
    struct Cluster {
        size_t begin = 0;
        size_t end = 0;
        double minimumX = 0;
        double maximumX = 0;
        double minimumY = 0;
        double maximumY = 0;
        int first = -1;
        int second = -1;
        bool mixedForms = false;  // holds points of several forms of the kernel (PointForm): never a leaf
    };
    struct Block {
        size_t rowCluster = 0;
        size_t columnCluster = 0;
        size_t rank = 0;      // low rank: U (rows x rank) V^T (rank x columns); 0 with `dense`
        bool lowRank = false;
        std::vector<double> u;      // column-major, rows x rank
        std::vector<double> v;      // column-major, columns x rank
        std::vector<double> dense;  // row-major, rows x columns
    };
    static std::vector<Cluster> build_clusters(const std::vector<std::array<double, 2>>& points, std::vector<size_t>& order, const PointForm& form);
    void build_blocks(size_t rowCluster, size_t columnCluster, const Entry& entry, const BlockEntries* enclosing);
    static void build_set_blocks(std::vector<HierarchicalMatrix>& matrices, size_t rowCluster, size_t columnCluster, const MultiEntry& entry,
                                 const BlockEntries* enclosing);
    bool admissible(size_t rowCluster, size_t columnCluster) const;
    bool split_rows(size_t rowCluster, size_t columnCluster) const;
    // The block's entry at (i, j), in block coordinates (i-th row and j-th column of the block in cluster order).
    using LocalEntry = std::function<double(size_t i, size_t j)>;
    bool approximate(Block& block, const LocalEntry& entry) const;

    size_t _rows = 0;
    size_t _columns = 0;
    double _tolerance = 0;
    std::vector<size_t> _rowOrder;     // cluster order -> row
    std::vector<size_t> _columnOrder;  // cluster order -> column
    std::vector<Cluster> _rowClusters;
    std::vector<Cluster> _columnClusters;
    std::vector<Block> _blocks;
    std::vector<size_t> _sparseRows;
    std::vector<size_t> _sparseColumns;
    std::vector<double> _sparseValues;
};

// Restarted GMRES, right-preconditioned, for A x = b: stops when |b - A x| <= tolerance |b| and throws when it does not
// get there within `maximumIterations`.
std::vector<std::complex<double>> solve_gmres(const std::function<void(const std::vector<std::complex<double>>&, std::vector<std::complex<double>>&)>& apply,
                                              const std::function<void(const std::vector<std::complex<double>>&, std::vector<std::complex<double>>&)>& precondition,
                                              const std::vector<std::complex<double>>& rightHandSide, double tolerance, size_t restart,
                                              size_t maximumIterations, size_t* iterations = nullptr);

using GmresPrecondition = std::function<void(const std::vector<std::complex<double>>&, std::vector<std::complex<double>>&)>;
// outputs[i] = A_s inputs[i] for s = systems[i]; the outputs come zeroed.
using GmresBatchApply = std::function<void(const std::vector<size_t>& systems, const std::vector<const std::vector<std::complex<double>>*>& inputs,
                                           const std::vector<std::vector<std::complex<double>>*>& outputs)>;
// The same GMRES on several systems A_s x_s = b_s at once, each with its own preconditioner, restarts and convergence: they
// share only the calls to `apply`, which takes every product the running systems need in one go. Each solution is the
// one solve_gmres gives for its system alone.
std::vector<std::vector<std::complex<double>>> solve_gmres_batch(const GmresBatchApply& apply, const std::vector<GmresPrecondition>& precondition,
                                                                 const std::vector<std::vector<std::complex<double>>>& rightHandSides, double tolerance,
                                                                 size_t restart, size_t maximumIterations, std::vector<size_t>* iterations = nullptr);

} // namespace OpenMagnetics
