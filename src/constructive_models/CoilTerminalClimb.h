#pragma once
// ABT #2274 item E (Alf, "MKF plans what is built"; 2026-11-13 "all parallels terminals should be continuous"):
// the climb a RAMPED terminal lead (ConnectionRoute::rampLength: its stub is shorter than the 2R its two bends
// need) is actually BUILT with. MVB++ never draws the cosine S that rampLength describes: it draws the corner
// between the lead and its turn as the terminal fillet of src/TerminalFillet.cpp (filletHelixCorners, default
// roll): a planar bend of the bend radius in the helix's natural plane plus a tangent biarc to the lead, with a
// one-OD stub too short for its own corner consumed into it. This header is that construction, moved into MKF
// step for step (same tangent-length walk, same phi list, same snap of the bend plane, same acceptance), in the
// MVB++ concentric frame (column axis y, connection face -z, P(az) = (r cos az, y, -r sin az)).
// What a lane must keep clear of is then measured on the built copper (terminal_climb_span), not on a reserved
// ramp nobody draws.
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace OpenMagnetics::terminal_climb {

struct V {
    double x = 0, y = 0, z = 0;
    V operator+(const V& o) const { return {x + o.x, y + o.y, z + o.z}; }
    V operator-(const V& o) const { return {x - o.x, y - o.y, z - o.z}; }
    V operator*(double s) const { return {x * s, y * s, z * s}; }
    double dot(const V& o) const { return x * o.x + y * o.y + z * o.z; }
    V cross(const V& o) const { return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x}; }
    double norm() const { return std::sqrt(dot(*this)); }
    V unit() const { const double n = norm(); return n > 0 ? *this * (1.0 / n) : *this; }
};
inline double angle(const V& a, const V& b) {
    // gp_Vec::Angle: the angle between two vectors, in [0, pi].
    const double c = a.cross(b).norm(), d = a.dot(b);
    return std::atan2(c, d);
}
// p rotated by `ang` about the axis through `c` along unit `ax` (right-handed).
inline V rotate(const V& p, const V& c, const V& ax, double ang) {
    const V v = p - c;
    const double co = std::cos(ang), si = std::sin(ang);
    return c + v * co + ax.cross(v) * si + ax * (ax.dot(v) * (1 - co));
}

struct Arc {
    bool straight = false;
    V start, end;       // centreline end points, in travel order
    V c, axis, v0;      // centre, unit axis, start - centre
    double sweep = 0, radius = std::numeric_limits<double>::infinity();
    V tEnd;
    V at(double f) const { return straight ? start + (end - start) * f : rotate(start, c, axis, sweep * f); }
};

// TerminalFillet.cpp arcFromTangent: the circular arc that leaves P along unit t and reaches Q.
inline Arc arc_from_tangent(const V& P, const V& t, const V& Q) {
    Arc out;
    out.start = P;
    out.end = Q;
    const V q = Q - P;
    V n = q - t * q.dot(t);
    if (n.norm() < 1e-15 * std::max(q.norm(), 1e-12)) {
        out.straight = true;
        out.tEnd = t;
        return out;
    }
    const double R = q.dot(q) / (2.0 * n.norm());
    n = n.unit();
    out.c = P + n * R;
    out.axis = t.cross(n).unit();
    out.v0 = P - out.c;
    out.sweep = angle(P - out.c, Q - out.c);
    out.radius = R;
    out.tEnd = out.axis.cross(Q - out.c).unit();
    if ((rotate(P, out.c, out.axis, out.sweep) - Q).norm() > 1e-9 * std::max(R, 1e-9)) {
        throw std::logic_error("terminal_climb: a fillet arc misses its landing");
    }
    return out;
}

constexpr double kMaxFilletSweep = 120.0 * std::numbers::pi / 180.0;   // TerminalFillet.cpp: a fillet never loops
inline bool loops(const Arc& a) { return !a.straight && a.sweep > kMaxFilletSweep; }

// TerminalFillet.cpp biarc: equal tangent lengths.
inline std::optional<std::array<Arc, 2>> biarc(const V& P1, const V& t1, const V& P2, const V& t2) {
    const V d = P2 - P1;
    const V s = t1 + t2;
    const double c = t1.dot(t2);
    const double qa = 2.0 * c - 2.0, qb = -2.0 * d.dot(s), qc = d.dot(d);
    double a = -1.0;
    if (std::abs(qa) < 1e-12) {
        if (std::abs(qb) < 1e-300) return std::nullopt;
        a = -qc / qb;
    }
    else {
        const double disc = qb * qb - 4.0 * qa * qc;
        if (disc < 0) return std::nullopt;
        const double r1 = (-qb - std::sqrt(disc)) / (2.0 * qa), r2 = (-qb + std::sqrt(disc)) / (2.0 * qa);
        a = std::numeric_limits<double>::max();
        for (double r : {r1, r2})
            if (r > 0 && r < a) a = r;
        if (a == std::numeric_limits<double>::max()) return std::nullopt;
    }
    const V A = P1 + t1 * a, B = P2 - t2 * a;
    const V J = (A + B) * 0.5;
    const Arc a1 = arc_from_tangent(P1, t1, J);
    const V tJ = B - A;
    if (tJ.norm() < 1e-15) return std::nullopt;
    const Arc a2 = arc_from_tangent(J, tJ.unit(), P2);
    if (angle(a2.tEnd, t2) > 1e-9) return std::nullopt;
    if (loops(a1) || loops(a2)) return std::nullopt;
    return std::array<Arc, 2>{a1, a2};
}

// TerminalFillet.cpp triarcHelixPlane.
inline std::optional<std::vector<Arc>> triarc(const V& P1, const V& t1, const V& P2, const V& t2, bool helixAtP2,
                                              const V& nHelix, double R, double phi) {
    const V PH = helixAtP2 ? P2 : P1;
    const V tH = helixAtP2 ? t2 : t1;
    const V PL = helixAtP2 ? P1 : P2;
    if (nHelix.norm() < 1e-30 || tH.cross(nHelix).norm() < 1e-30) return std::nullopt;
    const double side = (PL - PH).dot(nHelix) >= 0 ? 1.0 : -1.0;
    const V n = nHelix * side;
    const V C = PH + n * R;
    const V ax = tH.cross(n).unit();
    const V origin{0, 0, 0};
    if (helixAtP2) {
        const V J = rotate(P2, C, ax, -phi);
        const V tJ = rotate(tH, origin, ax, -phi).unit();
        const auto b = biarc(P1, t1, J, tJ);
        if (!b) return std::nullopt;
        const Arc a3 = arc_from_tangent(J, tJ, P2);
        if (angle(a3.tEnd, t2) > 1e-9) return std::nullopt;
        return std::vector<Arc>{(*b)[0], (*b)[1], a3};
    }
    const V J = rotate(P1, C, ax, phi);
    const V tJ = rotate(tH, origin, ax, phi).unit();
    const Arc a1 = arc_from_tangent(P1, t1, J);
    if (angle(a1.tEnd, tJ) > 1e-9) return std::nullopt;
    const auto b = biarc(J, tJ, P2, t2);
    if (!b) return std::nullopt;
    return std::vector<Arc>{a1, (*b)[0], (*b)[1]};
}

// The attach turn's helix as MKF lays it at the lane (terminal_turn_height_along: station - advance * along / L).
struct Helix {
    double radius;        // the attach turn's centreline radius (round column)
    double heightAtLane;  // y of the turn at az(lane)
    double dydaz;         // axial travel per radian of az (signed)
    int sgn;              // +1: the wire travels towards increasing az
    double cz = 0;        // the circle's centre in z: a dragback ride pushes the connection face out (TerminalTurnPath::rideOwnFace)
    double y(double az, double azLane) const { return heightAtLane + dydaz * (az - azLane); }
    V p(double az, double azLane) const { return {radius * std::cos(az), y(az, azLane), cz - radius * std::sin(az)}; }
    V t(double az) const {   // unit TRAVEL tangent
        return (V{-radius * std::sin(az), dydaz, -radius * std::cos(az)} * double(sgn)).unit();
    }
    V principalNormal(double az) const {   // analytic: the tangent-orthogonal part of d2P/daz2
        const V dd{-radius * std::cos(az), 0, radius * std::sin(az)};
        const V tt = t(az);
        return (dd - tt * dd.dot(tt)).unit();
    }
};

struct Climb {
    std::vector<Arc> arcs;     // travel order
    double tangentLength = 0;  // L taken from the lead leg and (as arc length) from the helix
    double phiDegrees = 0;
    bool stubConsumed = false;
};

// The terminal corner of a lead whose run reaches its lane at `rowY` and whose axial stub climbs from the row to
// the turn at the lane, as MVB++ builds it at the default roll. entrance: lead then helix (travel), exit: helix then
// lead. runLength: the radial run's length (its usable part is what lies beyond its far corner's mitre, as
// TerminalFillet's legInfo). wireRadius: the radius TerminalFillet is given for that mitre. Throws when no fillet
// of radius >= bend exists (the build would throw too).
inline Climb plan(bool entrance, const Helix& helix, double laneX, double rowY, double runLength, double bend,
                  double wireRadius, const std::string& who) {
    constexpr double pi = std::numbers::pi;
    if (!(std::abs(laneX) < helix.radius)) {
        throw std::logic_error("terminal_climb: the lane of " + who + " lies beyond its turn's radius");
    }
    const double azLane = std::acos(laneX / helix.radius);
    const V attach = helix.p(azLane, azLane);
    const double stub = std::abs(helix.heightAtLane - rowY);
    const V up{0, helix.heightAtLane >= rowY ? 1.0 : -1.0, 0};   // from the row towards the turn
    const V inward{0, 0, 1};                                        // the run, travelling towards the turn
    const V corner{laneX, rowY, attach.z};                          // run meets stub
    // Legs in TRAVEL order, adjacent to the joint first: the stub, then the run.
    struct Leg { V t; double len; double usable; V anchor; };
    auto usableOf = [&](double len, const V& tk, const V& tOther) {
        const double phiC = angle(tk, tOther);
        return (phiC > 1e-9 && phiC < pi - 1e-9) ? std::max(0.0, len - (wireRadius * std::tan(0.5 * phiC) + wireRadius)) : len;
    };
    const V tStub = entrance ? up : up * -1.0;
    const V tRun = entrance ? inward : inward * -1.0;
    Leg stubLeg{tStub, stub, usableOf(stub, tStub, tRun), entrance ? attach : attach};
    Leg runLeg{tRun, runLength, usableOf(runLength, tRun, V{0, 0, 0}), corner};
    // TerminalFillet: the lead-side anchor is the leg's end at the joint (entrance seg.b, exit seg.a).
    stubLeg.anchor = attach;
    runLeg.anchor = corner;
    const V tJoint = helix.t(azLane);
    auto thetaOf = [&](const V& tLeg) { return entrance ? angle(tLeg, tJoint) : angle(tJoint, tLeg); };
    Leg leg = stubLeg;
    bool consumed = false;
    if (!(stub > 1e-12)) {
        leg = runLeg;
    }
    else {
        const double L0 = std::max(bend * std::tan(0.5 * thetaOf(stubLeg.t)), 1e-9);
        if (0.9 * stubLeg.usable < L0) {
            leg = runLeg;
            consumed = true;
        }
    }
    const double theta = thetaOf(leg.t);
    if (theta < 1e-9) {
        throw std::logic_error("terminal_climb: the terminal corner of " + who + " is already tangent");
    }
    static const double kPhiDeg[] = {5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 85};
    const double dsdaz = V{-helix.radius * std::sin(azLane), helix.dydaz, -helix.radius * std::cos(azLane)}.norm();
    double L = std::max(bend * std::tan(0.5 * theta), 1e-9);
    for (int it = 0; it < 24; ++it, L *= 1.2) {
        if (L > 0.9 * leg.usable) break;
        const double dAz = L / std::max(dsdaz, 1e-18);
        const double azCut = entrance ? azLane + helix.sgn * dAz : azLane - helix.sgn * dAz;
        const V H = helix.p(azCut, azLane), tH = helix.t(azCut);
        V P1, t1, P2, t2;
        if (entrance) { P1 = leg.anchor - leg.t * L; t1 = leg.t; P2 = H; t2 = tH; }
        else          { P1 = H; t1 = tH; P2 = leg.anchor + leg.t * L; t2 = leg.t; }
        {   // TerminalFillet's planar-corner shortcut applies only when the corner is exactly coplanar.
            const V nPlane = t1.cross(t2);
            if (nPlane.norm() > 1e-12 && std::abs((P2 - P1).dot(nPlane.unit())) < 1e-12) {
                throw std::logic_error("terminal_climb: the terminal corner of " + who +
                                       " is planar (TerminalFillet's single-arc corner), which this plan does not model");
            }
        }
        V nH = helix.principalNormal(azCut);
        {   // the snap to the helix's natural plane nearest the lead (roll NaN)
            const V& PHx = entrance ? P2 : P1;
            const V& PLx = entrance ? P1 : P2;
            V toLead = PLx - PHx;
            toLead = toLead - tH * toLead.dot(tH);
            if (toLead.norm() > 1e-12) {
                toLead = toLead.unit();
                const V B = tH.cross(nH).unit();
                const double cN = toLead.dot(nH), cB = toLead.dot(B);
                nH = std::abs(cN) >= std::abs(cB) ? nH * (cN >= 0 ? 1.0 : -1.0) : B * (cB >= 0 ? 1.0 : -1.0);
            }
        }
        for (double phiDeg : kPhiDeg) {
            auto f = triarc(P1, t1, P2, t2, entrance, nH, bend, phiDeg * pi / 180.0);
            if (!f) continue;
            double minRadius = std::numeric_limits<double>::infinity();
            for (const auto& a : *f) minRadius = std::min(minRadius, a.radius);
            if (!(minRadius >= bend * (1.0 - 1e-9))) continue;
            Climb out;
            out.arcs = std::move(*f);
            out.tangentLength = L;
            out.phiDegrees = phiDeg;
            out.stubConsumed = consumed;
            return out;
        }
    }
    throw std::logic_error("terminal_climb: no terminal fillet of radius >= " + std::to_string(bend * 1e3) + " mm fits the corner of " + who);
}

// How far beyond its lane (along the face, towards `outwards` = -1 lower x / +1 higher x) a LEVEL lead of coated
// diameter `clearance`'s partner at height `levelY` must stand from the climb so the two centrelines are at least
// `clearance` apart, minus nothing: the separation s* (metres, >= 0) from the lane. The level lead is taken as a
// straight line along z (the run on the face), the climb as its built centreline: the arcs, sampled and refined,
// and the stub when it is not consumed (vertical at the lane: s* >= clearance from it wherever it lies within
// `clearance` of the level height).
inline double level_separation(const Climb& climb, double laneX, int outwards, double levelY, double clearance) {
    double need = 0;
    auto take = [&](const V& p) {
        const double dy = p.y - levelY;
        if (std::abs(dy) < clearance) {
            need = std::max(need, outwards * (p.x - laneX) + std::sqrt(clearance * clearance - dy * dy));
        }
    };
    take(V{laneX, levelY, 0});   // the run itself at the lane
    for (const auto& a : climb.arcs) {
        constexpr int n = 4096;
        int bestK = 0;
        double before = need;
        for (int k = 0; k <= n; ++k) {
            const double prev = need;
            take(a.at(double(k) / n));
            if (need > prev) bestK = k;
        }
        if (need > before) {   // refine about the best sample (golden section on the piecewise-smooth need)
            double lo = std::max(0.0, (bestK - 1.0) / n), hi = std::min(1.0, (bestK + 1.0) / n);
            auto f = [&](double s) {
                const V p = a.at(s);
                const double dy = p.y - levelY;
                return std::abs(dy) < clearance ? outwards * (p.x - laneX) + std::sqrt(clearance * clearance - dy * dy)
                                                : -std::numeric_limits<double>::infinity();
            };
            const double g = (std::sqrt(5.0) - 1) / 2;
            double c = hi - g * (hi - lo), d = lo + g * (hi - lo);
            for (int it = 0; it < 100; ++it) {
                if (f(c) > f(d)) hi = d; else lo = c;
                c = hi - g * (hi - lo);
                d = lo + g * (hi - lo);
            }
            need = std::max(need, f((lo + hi) / 2));
        }
    }
    return need;
}

}  // namespace OpenMagnetics::terminal_climb
