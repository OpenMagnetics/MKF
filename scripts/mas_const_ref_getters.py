#!/usr/bin/env python3
"""Post-generation step for MKF's build of MAS.hpp (quicktype C++ output).

quicktype returns every std::optional member BY VALUE:

    std::optional<std::vector<double>> get_time() const { return time; }

so every presence test or read of an optional waveform, processed description, harmonics record, ... copies it
(a 16k-sample Waveform is ~128 KB per call). This step adds, next to each such getter and without touching it, a
const-reference accessor that returns the same member without the copy:

    const std::optional<std::vector<double>> & get_time_ref() const { return time; }

Every existing getter and the JSON (de)serialisers stay as generated, so code that does not use the new accessors
compiles and behaves exactly as before. The setters, `void set_x(T value) { this->x = value; }`, take their argument
by value and copy-assigned it a second time; they are rewritten to `this->x = std::move(value);`, which stores the
same value without the second copy (a caller passing an lvalue still copies once, into the argument). The schemas and the MAS repository
are not touched; this only edits MKF's generated build/MAS/MAS.hpp.

usage: mas_const_ref_getters.py <MAS.hpp>   (edits the file in place; idempotent)
"""
import re
import sys

SETTER_ASSIGNMENT = re.compile(r"this->(?P<member>\w+) = value; \}$")
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
    # Setters take their argument by value and then copy-assigned it; they now move it (same value, one copy less).
    moved = 0
    for index, line in enumerate(out):
        rewritten, count = SETTER_ASSIGNMENT.subn(r"this->\g<member> = std::move(value); }", line)
        if count:
            out[index] = rewritten
            moved += count
    # quicktype declares `virtual ~X() = default;` in every class, which suppresses the implicit move constructor and
    # move assignment: std::move of a generated object, and every setter that moves its argument, then COPIES it (a
    # 16k-sample waveform per call). The four special members are declared defaulted next to that destructor: copies
    # stay exactly what they are (member-wise), and a move now moves (member-wise), which is what the code that calls
    # std::move asked for. Nothing else changes; a class whose members cannot all be moved keeps copying (a defaulted
    # move that is deleted is ignored by overload resolution).
    DESTRUCTOR = re.compile(r"^(?P<indent>\s*)virtual ~(?P<name>\w+)\(\) = default;$")
    declared = {m.group(1) for line in out for m in [re.search(r"^\s*(\w+)\(\1 &&\) = default;", line)] if m}
    special = 0
    with_special = []
    for line in out:
        with_special.append(line)
        m = DESTRUCTOR.match(line)
        if not m or m.group("name") in declared:
            continue
        name, indent = m.group("name"), m.group("indent")
        with_special += [f"{indent}{name}(const {name} &) = default;", f"{indent}{name}({name} &&) = default;",
                         f"{indent}{name} & operator=(const {name} &) = default;",
                         f"{indent}{name} & operator=({name} &&) = default;"]
        special += 1
    out = with_special
    text = "\n".join(out)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"mas_const_ref_getters: {added} const& accessors added, {moved} setters made to move, "
          f"{special} classes given defaulted copy and move members, in {path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
