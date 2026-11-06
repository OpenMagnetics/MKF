#include "support/HierarchicalMatrix.h"
#include "support/Exceptions.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>

namespace OpenMagnetics {

namespace {
// Clusters of at most this many points are leaves.
constexpr size_t leafSize = 32;
// A pair of clusters is far when the larger of their diagonals is at most this times the distance between their boxes.
// A far pair that does not compress to the tolerance is split anyway (and its sampled entries guard the ones that do).
constexpr double admissibility = 2.0;
// Entries checked per low-rank block against the kernel.
constexpr size_t checkedEntries = 12;
} // namespace

std::vector<HierarchicalMatrix::Cluster> HierarchicalMatrix::build_clusters(const std::vector<std::array<double, 2>>& points, std::vector<size_t>& order) {
    order.resize(points.size());
    std::iota(order.begin(), order.end(), 0);
    std::vector<Cluster> clusters;
    if (points.empty()) {
        return clusters;
    }
    std::function<int(size_t, size_t)> build = [&](size_t begin, size_t end) -> int {
        Cluster cluster;
        cluster.begin = begin;
        cluster.end = end;
        cluster.minimumX = cluster.maximumX = points[order[begin]][0];
        cluster.minimumY = cluster.maximumY = points[order[begin]][1];
        for (size_t index = begin; index < end; ++index) {
            const auto& point = points[order[index]];
            if (!std::isfinite(point[0]) || !std::isfinite(point[1])) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "Hierarchical matrix: a point has no finite position");
            }
            cluster.minimumX = std::min(cluster.minimumX, point[0]);
            cluster.maximumX = std::max(cluster.maximumX, point[0]);
            cluster.minimumY = std::min(cluster.minimumY, point[1]);
            cluster.maximumY = std::max(cluster.maximumY, point[1]);
        }
        const int index = static_cast<int>(clusters.size());
        clusters.push_back(cluster);
        const double extentX = cluster.maximumX - cluster.minimumX;
        const double extentY = cluster.maximumY - cluster.minimumY;
        if (end - begin > leafSize && (extentX > 0 || extentY > 0)) {
            const size_t axis = extentX >= extentY ? 0 : 1;
            const size_t middle = begin + (end - begin) / 2;
            std::nth_element(order.begin() + static_cast<std::ptrdiff_t>(begin), order.begin() + static_cast<std::ptrdiff_t>(middle),
                             order.begin() + static_cast<std::ptrdiff_t>(end),
                             [&](size_t a, size_t b) { return points[a][axis] < points[b][axis] || (points[a][axis] == points[b][axis] && a < b); });
            const int first = build(begin, middle);
            const int second = build(middle, end);
            clusters[static_cast<size_t>(index)].first = first;
            clusters[static_cast<size_t>(index)].second = second;
        }
        return index;
    };
    build(0, points.size());
    return clusters;
}

// The kernel's entries of one block, kept as they are evaluated (in block coordinates, `count` values per entry), so no
// entry is evaluated twice. `keep` false keeps nothing (every get evaluates).
class BlockEntries {
  public:
    BlockEntries(size_t rows, size_t columns, size_t count, bool keep) : _columns(columns), _count(count), _keep(keep) {
        _scratch.resize(count);
        // A large block keeps the entries it evaluates by index (an approximation evaluates a few of its rows and columns).
        _dense = _keep && rows * columns <= maximumDenseEntries;
        if (_dense) {
            _known.assign(rows * columns, false);
            _values.resize(rows * columns * count);
        }
    }
    template <typename Evaluate> double get(size_t i, size_t j, size_t index, const Evaluate& evaluate) {
        if (!_keep) {
            evaluate(_scratch.data());
            return _scratch[index];
        }
        const size_t entry = i * _columns + j;
        if (_dense) {
            if (!_known[entry]) {
                evaluate(_values.data() + entry * _count);
                _known[entry] = true;
            }
            return _values[entry * _count + index];
        }
        auto found = _sparse.find(entry);
        if (found == _sparse.end()) {
            evaluate(_scratch.data());
            found = _sparse.emplace(entry, _scratch).first;
        }
        return found->second[index];
    }

  private:
    static constexpr size_t maximumDenseEntries = size_t(1) << 20;
    size_t _columns;
    size_t _count;
    bool _keep;
    bool _dense = false;
    std::vector<bool> _known;
    std::vector<double> _values;
    std::vector<double> _scratch;
    std::unordered_map<size_t, std::vector<double>> _sparse;
};

HierarchicalMatrix::HierarchicalMatrix(const std::vector<std::array<double, 2>>& rowPoints, const std::vector<std::array<double, 2>>& columnPoints,
                                       const Entry& entry, double tolerance)
    : _rows(rowPoints.size()), _columns(columnPoints.size()), _tolerance(tolerance) {
    if (!(tolerance > 0) || !(tolerance < 1)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Hierarchical matrix: the tolerance must lie between 0 and 1");
    }
    _rowClusters = build_clusters(rowPoints, _rowOrder);
    _columnClusters = build_clusters(columnPoints, _columnOrder);
    if (_rows > 0 && _columns > 0) {
        build_blocks(0, 0, entry);
    }
}

bool HierarchicalMatrix::admissible(size_t rowCluster, size_t columnCluster) const {
    const Cluster& rows = _rowClusters[rowCluster];
    const Cluster& columns = _columnClusters[columnCluster];
    const double rowDiagonal = std::hypot(rows.maximumX - rows.minimumX, rows.maximumY - rows.minimumY);
    const double columnDiagonal = std::hypot(columns.maximumX - columns.minimumX, columns.maximumY - columns.minimumY);
    const double gapX = std::max({0.0, rows.minimumX - columns.maximumX, columns.minimumX - rows.maximumX});
    const double gapY = std::max({0.0, rows.minimumY - columns.maximumY, columns.minimumY - rows.maximumY});
    const double distance = std::hypot(gapX, gapY);
    return distance > 0 && std::max(rowDiagonal, columnDiagonal) <= admissibility * distance;
}

// Which of the two clusters to split: the larger (a leaf is not split). Only called when one of them is not a leaf.
bool HierarchicalMatrix::split_rows(size_t rowCluster, size_t columnCluster) const {
    const Cluster& rows = _rowClusters[rowCluster];
    const Cluster& columns = _columnClusters[columnCluster];
    const double rowDiagonal = std::hypot(rows.maximumX - rows.minimumX, rows.maximumY - rows.minimumY);
    const double columnDiagonal = std::hypot(columns.maximumX - columns.minimumX, columns.maximumY - columns.minimumY);
    return rows.first >= 0 && (columns.first < 0 || rowDiagonal >= columnDiagonal);
}

void HierarchicalMatrix::build_blocks(size_t rowCluster, size_t columnCluster, const Entry& entry) {
    const Cluster& rows = _rowClusters[rowCluster];
    const Cluster& columns = _columnClusters[columnCluster];
    const size_t numberRows = rows.end - rows.begin;
    const size_t numberColumns = columns.end - columns.begin;
    const bool leaves = rows.first < 0 && columns.first < 0;
    // A pair of leaves that does not compress is stored dense: the entries its approximation evaluated are kept for that.
    BlockEntries entries(numberRows, numberColumns, 1, leaves);
    auto block_entry = [&](size_t i, size_t j) {
        return entries.get(i, j, 0, [&](double* values) { values[0] = entry(_rowOrder[rows.begin + i], _columnOrder[columns.begin + j]); });
    };
    if (admissible(rowCluster, columnCluster)) {
        Block block;
        block.rowCluster = rowCluster;
        block.columnCluster = columnCluster;
        if (approximate(block, block_entry)) {
            _blocks.push_back(std::move(block));
            return;
        }
    }
    if (leaves) {
        Block block;
        block.rowCluster = rowCluster;
        block.columnCluster = columnCluster;
        block.dense.resize(numberRows * numberColumns);
        for (size_t i = 0; i < numberRows; ++i) {
            for (size_t j = 0; j < numberColumns; ++j) {
                block.dense[i * numberColumns + j] = block_entry(i, j);
            }
        }
        _blocks.push_back(std::move(block));
        return;
    }
    if (split_rows(rowCluster, columnCluster)) {
        const int first = rows.first;
        const int second = rows.second;
        build_blocks(static_cast<size_t>(first), columnCluster, entry);
        build_blocks(static_cast<size_t>(second), columnCluster, entry);
    }
    else {
        const int first = columns.first;
        const int second = columns.second;
        build_blocks(rowCluster, static_cast<size_t>(first), entry);
        build_blocks(rowCluster, static_cast<size_t>(second), entry);
    }
}

std::vector<HierarchicalMatrix> HierarchicalMatrix::build_set(const std::vector<std::array<double, 2>>& rowPoints, const std::vector<std::array<double, 2>>& columnPoints,
                                                              const MultiEntry& entry, size_t count, double tolerance) {
    if (count == 0) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Hierarchical matrix: a set needs at least one matrix");
    }
    if (!(tolerance > 0) || !(tolerance < 1)) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Hierarchical matrix: the tolerance must lie between 0 and 1");
    }
    HierarchicalMatrix shape;
    shape._rows = rowPoints.size();
    shape._columns = columnPoints.size();
    shape._tolerance = tolerance;
    shape._rowClusters = build_clusters(rowPoints, shape._rowOrder);
    shape._columnClusters = build_clusters(columnPoints, shape._columnOrder);
    std::vector<HierarchicalMatrix> matrices(count, shape);
    if (shape._rows > 0 && shape._columns > 0) {
        build_set_blocks(matrices, 0, 0, entry);
    }
    return matrices;
}

void HierarchicalMatrix::build_set_blocks(std::vector<HierarchicalMatrix>& matrices, size_t rowCluster, size_t columnCluster, const MultiEntry& entry) {
    const HierarchicalMatrix& shape = matrices.front();
    const Cluster& rows = shape._rowClusters[rowCluster];
    const Cluster& columns = shape._columnClusters[columnCluster];
    const size_t count = matrices.size();
    const size_t numberRows = rows.end - rows.begin;
    const size_t numberColumns = columns.end - columns.begin;
    const bool leaves = rows.first < 0 && columns.first < 0;
    // The kernel gives every matrix's entry at once: the entries one approximation evaluates serve the others, and a pair
    // of leaves stored dense reuses them all.
    const bool admissiblePair = shape.admissible(rowCluster, columnCluster);
    BlockEntries entries(numberRows, numberColumns, count, leaves || (admissiblePair && count > 1));
    auto evaluate = [&](size_t i, size_t j) {
        return [&, i, j](double* values) { entry(shape._rowOrder[rows.begin + i], shape._columnOrder[columns.begin + j], values); };
    };
    if (admissiblePair) {
        // Every matrix in low rank, or none: the set keeps one partition.
        std::vector<Block> blocks(count);
        bool allLowRank = true;
        for (size_t index = 0; index < count && allLowRank; ++index) {
            blocks[index].rowCluster = rowCluster;
            blocks[index].columnCluster = columnCluster;
            allLowRank = matrices[index].approximate(blocks[index], [&](size_t i, size_t j) { return entries.get(i, j, index, evaluate(i, j)); });
        }
        if (allLowRank) {
            for (size_t index = 0; index < count; ++index) {
                matrices[index]._blocks.push_back(std::move(blocks[index]));
            }
            return;
        }
    }
    if (leaves) {
        std::vector<Block> blocks(count);
        for (auto& block : blocks) {
            block.rowCluster = rowCluster;
            block.columnCluster = columnCluster;
            block.dense.resize(numberRows * numberColumns);
        }
        for (size_t i = 0; i < numberRows; ++i) {
            for (size_t j = 0; j < numberColumns; ++j) {
                for (size_t index = 0; index < count; ++index) {
                    blocks[index].dense[i * numberColumns + j] = entries.get(i, j, index, evaluate(i, j));
                }
            }
        }
        for (size_t index = 0; index < count; ++index) {
            matrices[index]._blocks.push_back(std::move(blocks[index]));
        }
        return;
    }
    if (shape.split_rows(rowCluster, columnCluster)) {
        build_set_blocks(matrices, static_cast<size_t>(rows.first), columnCluster, entry);
        build_set_blocks(matrices, static_cast<size_t>(rows.second), columnCluster, entry);
    }
    else {
        build_set_blocks(matrices, rowCluster, static_cast<size_t>(columns.first), entry);
        build_set_blocks(matrices, rowCluster, static_cast<size_t>(columns.second), entry);
    }
}

// Adaptive cross approximation with partial pivoting. Returns false when the block is not of low rank to the tolerance
// (its rank would reach half its smaller side) or when the approximation misses sampled entries of the kernel.
bool HierarchicalMatrix::approximate(Block& block, const LocalEntry& entry) const {
    const Cluster& rows = _rowClusters[block.rowCluster];
    const Cluster& columns = _columnClusters[block.columnCluster];
    const size_t m = rows.end - rows.begin;
    const size_t n = columns.end - columns.begin;
    const size_t maximumRank = std::min(m, n) / 2;
    if (maximumRank == 0) {
        return false;
    }
    std::vector<double> u;
    std::vector<double> v;
    std::vector<bool> usedRow(m, false);
    size_t rank = 0;
    double normSquared = 0;
    size_t pivotRow = 0;
    size_t consecutiveSmall = 0;
    std::vector<double> rowValues(n);
    std::vector<double> columnValues(m);
    while (rank < maximumRank) {
        usedRow[pivotRow] = true;
        for (size_t j = 0; j < n; ++j) {
            double value = entry(pivotRow, j);
            for (size_t k = 0; k < rank; ++k) {
                value -= u[k * m + pivotRow] * v[k * n + j];
            }
            rowValues[j] = value;
        }
        size_t pivotColumn = 0;
        for (size_t j = 1; j < n; ++j) {
            if (std::abs(rowValues[j]) > std::abs(rowValues[pivotColumn])) {
                pivotColumn = j;
            }
        }
        const double pivot = rowValues[pivotColumn];
        if (pivot == 0) {
            // This row is already represented: try the next unused one.
            auto next = std::find(usedRow.begin(), usedRow.end(), false);
            if (next == usedRow.end()) {
                break;
            }
            pivotRow = static_cast<size_t>(next - usedRow.begin());
            continue;
        }
        for (size_t i = 0; i < m; ++i) {
            double value = entry(i, pivotColumn);
            for (size_t k = 0; k < rank; ++k) {
                value -= u[k * m + i] * v[k * n + pivotColumn];
            }
            columnValues[i] = value;
        }
        // New term u v^T, v = row / pivot, u = column.
        double uNorm = 0;
        double vNorm = 0;
        for (size_t i = 0; i < m; ++i) {
            uNorm += columnValues[i] * columnValues[i];
        }
        for (size_t j = 0; j < n; ++j) {
            rowValues[j] /= pivot;
            vNorm += rowValues[j] * rowValues[j];
        }
        double cross = 0;
        for (size_t k = 0; k < rank; ++k) {
            double uDot = 0;
            double vDot = 0;
            for (size_t i = 0; i < m; ++i) {
                uDot += columnValues[i] * u[k * m + i];
            }
            for (size_t j = 0; j < n; ++j) {
                vDot += rowValues[j] * v[k * n + j];
            }
            cross += uDot * vDot;
        }
        u.insert(u.end(), columnValues.begin(), columnValues.end());
        v.insert(v.end(), rowValues.begin(), rowValues.end());
        ++rank;
        normSquared += uNorm * vNorm + 2 * cross;
        const double termNorm = std::sqrt(uNorm * vNorm);
        // Two terms in a row below the tolerance: one small term can come from a pivot that happens to be poor.
        if (termNorm <= _tolerance * std::sqrt(std::max(normSquared, 0.0))) {
            if (++consecutiveSmall >= 2) {
                break;
            }
        }
        else {
            consecutiveSmall = 0;
        }
        // Next pivot row: the largest of the new column among the unused rows.
        double best = -1;
        bool found = false;
        for (size_t i = 0; i < m; ++i) {
            if (!usedRow[i] && std::abs(columnValues[i]) > best) {
                best = std::abs(columnValues[i]);
                pivotRow = i;
                found = true;
            }
        }
        if (!found) {
            break;
        }
    }
    if (rank >= maximumRank) {
        return false;
    }
    // Check sampled entries against the kernel: ACA estimates its error from the terms it found, and a part of the block
    // its pivots never touched would be missed silently.
    const double entryScale = std::sqrt(std::max(normSquared, 0.0) / static_cast<double>(m * n));
    std::mt19937_64 generator(static_cast<uint64_t>(block.rowCluster) * 1000003u + block.columnCluster);
    for (size_t sample = 0; sample < checkedEntries; ++sample) {
        const size_t i = generator() % m;
        const size_t j = generator() % n;
        double approximation = 0;
        for (size_t k = 0; k < rank; ++k) {
            approximation += u[k * m + i] * v[k * n + j];
        }
        const double exact = entry(i, j);
        if (std::abs(exact - approximation) > 100 * _tolerance * std::max(entryScale, std::abs(exact))) {
            return false;
        }
    }
    block.lowRank = true;
    block.rank = rank;
    block.u = std::move(u);
    block.v = std::move(v);
    return true;
}

void HierarchicalMatrix::add_sparse(size_t row, size_t column, double value) {
    if (row >= _rows || column >= _columns) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Hierarchical matrix: a sparse correction lies outside the matrix");
    }
    _sparseRows.push_back(row);
    _sparseColumns.push_back(column);
    _sparseValues.push_back(value);
}

void HierarchicalMatrix::multiply_add(const std::vector<std::complex<double>>& x, std::vector<std::complex<double>>& y) const {
    multiply_add_batch({&x}, {&y});
}

void HierarchicalMatrix::multiply_add_batch(const std::vector<const std::vector<std::complex<double>>*>& x, const std::vector<std::vector<std::complex<double>>*>& y) const {
    const size_t count = x.size();
    if (y.size() != count) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "Hierarchical matrix: as many results as vectors are needed");
    }
    for (size_t vector = 0; vector < count; ++vector) {
        if (x[vector]->size() != _columns || y[vector]->size() != _rows) {
            throw InvalidInputException(ErrorCode::INVALID_INPUT, "Hierarchical matrix: the vectors do not match the matrix");
        }
    }
    if (count == 0) {
        return;
    }
    // In cluster order, real and imaginary parts apart, the vectors interleaved in lanes of `width` (entry j of every vector
    // side by side, the last lane padded with zeros): every block then reads and writes contiguous runs, is read once for
    // all the vectors, and works on whole lanes, so a vector gets the same operations in the same order whichever lane,
    // and whatever company, it has.
    constexpr size_t width = 4;
    const size_t lanes = (count + width - 1) / width;
    const size_t stride = lanes * width;
    std::vector<double> xReal(_columns * stride, 0.0);
    std::vector<double> xImaginary(_columns * stride, 0.0);
    for (size_t j = 0; j < _columns; ++j) {
        for (size_t vector = 0; vector < count; ++vector) {
            const auto value = (*x[vector])[_columnOrder[j]];
            xReal[j * stride + vector] = value.real();
            xImaginary[j * stride + vector] = value.imag();
        }
    }
    std::vector<double> yReal(_rows * stride, 0.0);
    std::vector<double> yImaginary(_rows * stride, 0.0);
    std::vector<double> projectedReal;
    std::vector<double> projectedImaginary;
    for (const auto& block : _blocks) {
        const Cluster& rows = _rowClusters[block.rowCluster];
        const Cluster& columns = _columnClusters[block.columnCluster];
        const size_t m = rows.end - rows.begin;
        const size_t n = columns.end - columns.begin;
        if (block.lowRank) {
            projectedReal.assign(block.rank * stride, 0.0);
            projectedImaginary.assign(block.rank * stride, 0.0);
            for (size_t lane = 0; lane < lanes; ++lane) {
                const double* inReal = xReal.data() + columns.begin * stride + lane * width;
                const double* inImaginary = xImaginary.data() + columns.begin * stride + lane * width;
                for (size_t k = 0; k < block.rank; ++k) {
                    const double* vk = block.v.data() + k * n;
                    double sumReal[width] = {};
                    double sumImaginary[width] = {};
                    for (size_t j = 0; j < n; ++j) {
                        const double value = vk[j];
                        for (size_t l = 0; l < width; ++l) {
                            sumReal[l] += value * inReal[j * stride + l];
                            sumImaginary[l] += value * inImaginary[j * stride + l];
                        }
                    }
                    for (size_t l = 0; l < width; ++l) {
                        projectedReal[k * stride + lane * width + l] = sumReal[l];
                        projectedImaginary[k * stride + lane * width + l] = sumImaginary[l];
                    }
                }
                double* outReal = yReal.data() + rows.begin * stride + lane * width;
                double* outImaginary = yImaginary.data() + rows.begin * stride + lane * width;
                for (size_t k = 0; k < block.rank; ++k) {
                    const double* uk = block.u.data() + k * m;
                    double pReal[width];
                    double pImaginary[width];
                    for (size_t l = 0; l < width; ++l) {
                        pReal[l] = projectedReal[k * stride + lane * width + l];
                        pImaginary[l] = projectedImaginary[k * stride + lane * width + l];
                    }
                    for (size_t i = 0; i < m; ++i) {
                        const double value = uk[i];
                        for (size_t l = 0; l < width; ++l) {
                            outReal[i * stride + l] += value * pReal[l];
                            outImaginary[i * stride + l] += value * pImaginary[l];
                        }
                    }
                }
            }
        }
        else {
            for (size_t lane = 0; lane < lanes; ++lane) {
                const double* inReal = xReal.data() + columns.begin * stride + lane * width;
                const double* inImaginary = xImaginary.data() + columns.begin * stride + lane * width;
                double* outReal = yReal.data() + rows.begin * stride + lane * width;
                double* outImaginary = yImaginary.data() + rows.begin * stride + lane * width;
                for (size_t i = 0; i < m; ++i) {
                    const double* row = block.dense.data() + i * n;
                    double sumReal[width] = {};
                    double sumImaginary[width] = {};
                    for (size_t j = 0; j < n; ++j) {
                        const double value = row[j];
                        for (size_t l = 0; l < width; ++l) {
                            sumReal[l] += value * inReal[j * stride + l];
                            sumImaginary[l] += value * inImaginary[j * stride + l];
                        }
                    }
                    for (size_t l = 0; l < width; ++l) {
                        outReal[i * stride + l] += sumReal[l];
                        outImaginary[i * stride + l] += sumImaginary[l];
                    }
                }
            }
        }
    }
    for (size_t vector = 0; vector < count; ++vector) {
        auto& out = *y[vector];
        const auto& in = *x[vector];
        for (size_t i = 0; i < _rows; ++i) {
            out[_rowOrder[i]] += std::complex<double>(yReal[i * stride + vector], yImaginary[i * stride + vector]);
        }
        for (size_t index = 0; index < _sparseValues.size(); ++index) {
            out[_sparseRows[index]] += _sparseValues[index] * in[_sparseColumns[index]];
        }
    }
}

size_t HierarchicalMatrix::stored_values() const {
    size_t values = _sparseValues.size();
    for (const auto& block : _blocks) {
        values += block.dense.size() + block.u.size() + block.v.size();
    }
    return values;
}

std::vector<std::complex<double>> solve_gmres(const std::function<void(const std::vector<std::complex<double>>&, std::vector<std::complex<double>>&)>& apply,
                                              const std::function<void(const std::vector<std::complex<double>>&, std::vector<std::complex<double>>&)>& precondition,
                                              const std::vector<std::complex<double>>& rightHandSide, double tolerance, size_t restart,
                                              size_t maximumIterations, size_t* iterations) {
    std::vector<size_t> iterationsPerSystem;
    auto solutions = solve_gmres_batch(
        [&](const std::vector<size_t>&, const std::vector<const std::vector<std::complex<double>>*>& inputs, const std::vector<std::vector<std::complex<double>>*>& outputs) {
            apply(*inputs[0], *outputs[0]);
        },
        {precondition}, {rightHandSide}, tolerance, restart, maximumIterations, &iterationsPerSystem);
    if (iterations) {
        *iterations = iterationsPerSystem[0];
    }
    return std::move(solutions[0]);
}

std::vector<std::vector<std::complex<double>>> solve_gmres_batch(const GmresBatchApply& apply, const std::vector<GmresPrecondition>& precondition,
                                                                 const std::vector<std::vector<std::complex<double>>>& rightHandSides, double tolerance,
                                                                 size_t restart, size_t maximumIterations, std::vector<size_t>* iterations) {
    using Complex = std::complex<double>;
    auto norm = [](const std::vector<Complex>& vector) {
        double sum = 0;
        for (const auto& value : vector) {
            sum += std::norm(value);
        }
        return std::sqrt(sum);
    };
    const size_t count = rightHandSides.size();
    if (precondition.size() != count) {
        throw InvalidInputException(ErrorCode::INVALID_INPUT, "GMRES: one preconditioner per system is needed");
    }
    // Each system runs the restarted GMRES on its own; they only take their products with A together. A system checks its
    // residual, then runs Arnoldi steps (each one product), then, at the end of a cycle, takes its true residual (one product).
    enum class Phase { Check, Arnoldi, Residual, Done };
    struct System {
        Phase phase = Phase::Check;
        std::vector<Complex> solution;
        std::vector<Complex> residual;
        std::vector<Complex> work;
        std::vector<Complex> preconditioned;
        double rightHandSideNorm = 0;
        size_t done = 0;
        size_t k = 0;
        std::vector<std::vector<Complex>> basis;
        std::vector<std::vector<Complex>> hessenberg;  // column k: h(0..k+1, k)
        std::vector<Complex> cosines;
        std::vector<Complex> sines;
        std::vector<Complex> g;
    };
    std::vector<System> systems(count);
    for (size_t index = 0; index < count; ++index) {
        auto& system = systems[index];
        const size_t n = rightHandSides[index].size();
        system.solution.assign(n, 0.0);
        system.residual = rightHandSides[index];  // x = 0
        system.work.assign(n, 0.0);
        system.preconditioned.assign(n, 0.0);
        system.rightHandSideNorm = norm(rightHandSides[index]);
        if (system.rightHandSideNorm == 0) {
            system.phase = Phase::Done;
        }
    }
    std::vector<size_t> requested;
    std::vector<const std::vector<Complex>*> inputs;
    std::vector<std::vector<Complex>*> outputs;
    while (true) {
        requested.clear();
        inputs.clear();
        outputs.clear();
        for (size_t index = 0; index < count; ++index) {
            auto& system = systems[index];
            if (system.phase == Phase::Check) {
                const double residualNorm = norm(system.residual);
                if (!std::isfinite(residualNorm)) {
                    throw NaNResultException("GMRES: the residual is not finite");
                }
                if (residualNorm <= tolerance * system.rightHandSideNorm) {
                    system.phase = Phase::Done;
                    continue;
                }
                if (system.done >= maximumIterations) {
                    throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "GMRES: no convergence to " + std::to_string(tolerance) + " in " +
                                               std::to_string(maximumIterations) + " iterations (residual " +
                                               std::to_string(residualNorm / system.rightHandSideNorm) + ")");
                }
                // Arnoldi on A M^-1, from the current residual.
                system.basis.clear();
                system.basis.reserve(restart + 1);
                system.basis.push_back(system.residual);
                for (auto& value : system.basis.back()) {
                    value /= residualNorm;
                }
                system.hessenberg.clear();
                system.cosines.clear();
                system.sines.clear();
                system.g = {Complex(residualNorm, 0)};
                system.k = 0;
                system.phase = Phase::Arnoldi;
            }
            if (system.phase == Phase::Arnoldi) {
                precondition[index](system.basis[system.k], system.preconditioned);
                std::fill(system.work.begin(), system.work.end(), Complex(0));
                requested.push_back(index);
                inputs.push_back(&system.preconditioned);
                outputs.push_back(&system.work);
            }
            else if (system.phase == Phase::Residual) {
                std::fill(system.work.begin(), system.work.end(), Complex(0));
                requested.push_back(index);
                inputs.push_back(&system.solution);
                outputs.push_back(&system.work);
            }
        }
        if (requested.empty()) {
            break;
        }
        apply(requested, inputs, outputs);
        for (size_t index : requested) {
            auto& system = systems[index];
            const size_t n = system.work.size();
            auto& work = system.work;
            if (system.phase == Phase::Residual) {
                // The true residual, so convergence is judged on b - A x itself.
                for (size_t entry = 0; entry < n; ++entry) {
                    system.residual[entry] = rightHandSides[index][entry] - work[entry];
                }
                system.phase = Phase::Check;
                continue;
            }
            const size_t k = system.k;
            auto& basis = system.basis;
            auto& cosines = system.cosines;
            auto& sines = system.sines;
            auto& g = system.g;
            std::vector<Complex> column(k + 2, 0.0);
            for (size_t i = 0; i <= k; ++i) {  // modified Gram-Schmidt
                Complex dot = 0;
                for (size_t entry = 0; entry < n; ++entry) {
                    dot += std::conj(basis[i][entry]) * work[entry];
                }
                column[i] = dot;
                for (size_t entry = 0; entry < n; ++entry) {
                    work[entry] -= dot * basis[i][entry];
                }
            }
            const double nextNorm = norm(work);
            column[k + 1] = nextNorm;
            // Apply the previous rotations, then make this one.
            for (size_t i = 0; i < k; ++i) {
                const Complex top = std::conj(cosines[i]) * column[i] + std::conj(sines[i]) * column[i + 1];
                const Complex bottom = -sines[i] * column[i] + cosines[i] * column[i + 1];
                column[i] = top;
                column[i + 1] = bottom;
            }
            const double magnitude = std::sqrt(std::norm(column[k]) + std::norm(column[k + 1]));
            if (magnitude == 0) {
                throw CalculationException(ErrorCode::CALCULATION_INVALID_RESULT, "GMRES: breakdown (the system is singular)");
            }
            const Complex cosine = column[k] / magnitude;
            const Complex sine = column[k + 1] / magnitude;
            cosines.push_back(cosine);
            sines.push_back(sine);
            column[k] = magnitude;
            column[k + 1] = 0;
            g.push_back(-sine * g[k]);
            g[k] = std::conj(cosine) * g[k];
            system.hessenberg.push_back(std::move(column));
            const bool converged = std::abs(g[k + 1]) <= tolerance * system.rightHandSideNorm;
            if (nextNorm > 0 && !converged) {
                basis.push_back(work);
                for (auto& value : basis.back()) {
                    value /= nextNorm;
                }
            }
            ++system.k;
            ++system.done;
            if (!converged && nextNorm != 0 && system.k < restart && system.done < maximumIterations) {
                continue;  // the next Arnoldi step
            }
            // End of the cycle: back substitution for y, then x += M^-1 (V y).
            const size_t steps = system.k;
            std::vector<Complex> y(steps, 0.0);
            for (size_t i = steps; i-- > 0;) {
                Complex sum = g[i];
                for (size_t j = i + 1; j < steps; ++j) {
                    sum -= system.hessenberg[j][i] * y[j];
                }
                y[i] = sum / system.hessenberg[i][i];
            }
            std::fill(work.begin(), work.end(), Complex(0));
            for (size_t i = 0; i < steps; ++i) {
                for (size_t entry = 0; entry < n; ++entry) {
                    work[entry] += y[i] * basis[i][entry];
                }
            }
            precondition[index](work, system.preconditioned);
            for (size_t entry = 0; entry < n; ++entry) {
                system.solution[entry] += system.preconditioned[entry];
            }
            system.phase = Phase::Residual;
        }
    }
    if (iterations) {
        iterations->resize(count);
        for (size_t index = 0; index < count; ++index) {
            (*iterations)[index] = systems[index].done;
        }
    }
    std::vector<std::vector<Complex>> solutions(count);
    for (size_t index = 0; index < count; ++index) {
        solutions[index] = std::move(systems[index].solution);
    }
    return solutions;
}

} // namespace OpenMagnetics
