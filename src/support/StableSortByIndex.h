#pragma once

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <utility>
#include <vector>

namespace OpenMagnetics {

// std::stable_sort of a vector of heavy items, done on their indices.
//
// The generated MAS classes (CoreMaterial, Magnetic, Mas, ...) declare a destructor, so they
// have no implicit move operations: every "move" std::stable_sort makes is a deep copy (a
// CoreMaterial carries its whole permeability and loss tables). Sorting the pairs directly
// costs O(n log n) of those copies. Sorting the indices with the same comparator and then
// gathering the items costs n copies and gives exactly the same order: std::stable_sort's
// result depends only on the comparator's answers, and those are the same answers about the
// same items. (ABT #1449: the adviser spent ~10% of a CLLC resonant-inductor run in these
// sorts.)
template <typename Item, typename Compare>
void stable_sort_by_index(std::vector<Item>& items, Compare compare) {
    std::vector<size_t> order(items.size());
    std::iota(order.begin(), order.end(), size_t(0));
    std::stable_sort(order.begin(), order.end(), [&items, &compare](size_t left, size_t right) {
        return compare(items[left], items[right]);
    });
    bool alreadyInOrder = true;
    for (size_t position = 0; position < order.size(); ++position) {
        if (order[position] != position) {
            alreadyInOrder = false;
            break;
        }
    }
    if (alreadyInOrder) {
        return;
    }
    std::vector<Item> sorted;
    sorted.reserve(items.size());
    for (size_t index : order) {
        sorted.push_back(std::move(items[index]));
    }
    items.swap(sorted);
}

} // namespace OpenMagnetics
