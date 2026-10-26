#!/usr/bin/env python3
"""Post-generation step for MKF's build of MAS.hpp (quicktype C++ output).

quicktype returns every std::optional member BY VALUE:

    std::optional<std::vector<double>> get_time() const { return time; }

so every presence test or read of an optional waveform, processed description, harmonics record, ... copies it
(a 16k-sample Waveform is ~128 KB per call). This step adds, next to each such getter and without touching it, a
const-reference accessor that returns the same member without the copy:

    const std::optional<std::vector<double>> & get_time_ref() const { return time; }

Nothing is removed or changed: every existing getter, setter and the JSON (de)serialisers stay as generated, so
code that does not use the new accessors compiles and behaves exactly as before. The schemas and the MAS repository
are not touched; this only edits MKF's generated build/MAS/MAS.hpp.

usage: mas_const_ref_getters.py <MAS.hpp>   (edits the file in place; idempotent)
"""
import re
import sys

GETTER = re.compile(r"^(?P<indent>\s*)std::optional<(?P<type>.+)> get_(?P<name>\w+)\(\) const \{ return (?P<member>\w+); \}$")


def main(path):
    with open(path, encoding="utf-8") as f:
        lines = f.read().split("\n")
    existing = {m.group(1) for line in lines for m in [re.search(r" get_(\w+)_ref\(\) const", line)] if m}
    out = []
    added = 0
    for line in lines:
        out.append(line)
        m = GETTER.match(line)
        if not m or m.group("name") in existing:
            continue
        out.append(f"{m.group('indent')}const std::optional<{m.group('type')}> & get_{m.group('name')}_ref() const "
                   f"{{ return {m.group('member')}; }}")
        added += 1
    text = "\n".join(out)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"mas_const_ref_getters: {added} const& accessors added to {path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
