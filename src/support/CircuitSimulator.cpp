#include <vector>
#include <string>
#include <cstdio>
#include <cmath>

// Taken from https://gist.github.com/signaldust/74ce12ae67bf21a8141f9b4a19ce8365

// set to 1 to make LU factorization show pivots
#define VERBOSE_LU 0

// gMin for diodes etc..
static const double gMin = 1e-12;

// voltage tolerance
static const double vTolerance = 5e-5;

// thermal voltage for diodes/transistors
static const double vThermal = 0.026;

static const unsigned maxIter = 200;

//
// General overview
// ----------------
//
// Circuits are built from nodes and Components, where nodes are
// simply positive integers (with 0 designating ground).
//
// Every Component has one or more pins connecting to the circuit
// nodes as well as zero or more internal nets.
//
// While we map pins directly to nets here, the separation would
// be useful if the solver implemented stuff like net-reordering.
//
// MNACell represents a single entry in the solution matrix,
// where we store constants and time-step dependent constants
// separately, plus collect pointers to dynamic variables.
//
// We track enough information here that we only need to stamp once.
//
struct MNACell
{
    double g;       // simple values (eg. resistor conductance)
    double gtimed;  // time-scaled values (eg. capacitor conductance)

    // pointers to dynamic variables, added in once per solve
    std::vector<double*>    gdyn;

    double  lu, prelu;  // lu-solver values and matrix pre-LU cache

    std::string txt;    // text version of MNA for pretty-printing

    void clear()
    {
        g = 0;
        gtimed = 0;
        txt = "";
    }

    void initLU(double stepScale)
    {
        prelu = g + gtimed * stepScale;
    }

    // restore matrix state and update dynamic values
    void updatePre()
    {
        lu = prelu;
        for(size_t i = 0; i < gdyn.size(); ++i)
        {
            lu += 0[gdyn[i]];
        }
    }
};

// this is for keeping track of node information
// for the purposes of more intelligent plotting
struct MNANodeInfo
{
    enum Type
    {
        tVoltage,
        tCurrent,

        tCount
    };

    Type    type;   // one auto-range per unit-type
    double  scale;  // scale factor (eg. charge to voltage)

    std::string name;   // node name for display
};

// Stores A and b for A*x - b = 0, where x is the solution.
//
// A is stored as a vector of rows, for easy in-place pivots
//
struct MNASystem
{
    typedef std::vector<MNACell>    MNAVector;
    typedef std::vector<MNAVector>  MNAMatrix;

    // node names - for output
    std::vector<MNANodeInfo>    nodes;

    MNAMatrix   A;
    MNAVector   b;

    double      time;

    void setSize(size_t n)
    {
        A.resize(n);
        b.resize(n);

        nodes.resize(n);

        for(unsigned i = 0; i < n; ++i)
        {
            b[i].clear();
            A[i].resize(n);

            char buf[16];
            sprintf(buf, "v%d", i);
            nodes[i].name = buf;
            nodes[i].type = MNANodeInfo::tVoltage;
            nodes[i].scale = 1;

            for(unsigned j = 0; j < n; ++j)
            {
                A[i][j].clear();
            }
        }

        time = 0;
    }

    void stampTimed(double g, int r, int c, const std::string & txt)
    {
        A[r][c].gtimed += g;
        A[r][c].txt += txt;
    }

    void stampStatic(double g, int r, int c, const std::string & txt)
    {
        A[r][c].g += g;
        A[r][c].txt += txt;
    }
};

struct IComponent
{
    virtual ~IComponent() {}

    // return the number of pins for this component
    virtual int pinCount() = 0;

    // return a pointer to array of pin locations
    // NOTE: these will eventually be GUI locations to be unified
    virtual int * getPinLocs() = 0;

    // setup pins and calculate the size of the full netlist
    // the Component<> will handle this automatically
    //
    //  - netSize is the current size of the netlist
    //  - pins is an array of circuits nodes
    //
    virtual void setupNets(int & netSize, int & states, int * pins) = 0;

    // stamp constants into the matrix
    virtual void stamp(MNASystem & m) = 0;

    // this is for allocating state variables
    virtual void setupStates([[maybe_unused]] int & states) {}

    // update state variables, only tagged nodes
    // this is intended for fixed-time compatible
    // testing to make sure we can code-gen stuff
    virtual void update([[maybe_unused]] MNASystem & m) {}

    // return true if we're done - will keep iterating
    // until all the components are happy
    virtual bool newton([[maybe_unused]] MNASystem & m) { return true; }

    // time-step change, for caps to fix their state-variables
    virtual void scaleTime([[maybe_unused]] double told_per_new) {}
};

template <int nPins = 0, int nInternalNets = 0>
struct Component : IComponent
{
    static const int nNets = nPins + nInternalNets;

    int pinLoc[nPins];
    int nets[nNets];

    int pinCount() { return nPins; }

    int * getPinLocs() { return pinLoc; }

    void setupNets(int & netSize, int & states, int * pins)
    {
        for(int i = 0; i < nPins; ++i)
        {
            nets[i] = pins[i];
        }

        for(int i = 0; i < nInternalNets; ++i)
        {
            nets[nPins + i] = netSize++;
        }

        setupStates(states);
    }
};

static const int unitValueOffset = 4;
static const int unitValueMax = 8;
static const char * unitValueSuffixes[unitValueMax+1] = {
    "p", "n", "u", "m", "", "k", "M", "G"
};

static void formatUnitValue(char * buf, double v, const char * unit)
{
    int suff = unitValueOffset + int(log(v) / log(10.)) / 3;

    if(v < 1) suff -= 1;

    if(suff < 0) suff = 0;
    if(suff > unitValueMax) suff = unitValueMax;

    double vr = v / pow(10., 3*double(suff - unitValueOffset));

    sprintf(buf, "%.0f%s%s", vr, unitValueSuffixes[suff], unit);
}

struct Resistor : Component<2>
{
    double  r;

    Resistor(double r, int l0, int l1) : r(r)
    {
        pinLoc[0] = l0;
        pinLoc[1] = l1;
    }

    void stamp(MNASystem & m)
    {
        char txt[16];
        txt[0] = 'R';
        formatUnitValue(txt+1, r, "");
        double g = 1. / r;
        m.stampStatic(+g, nets[0], nets[0], std::string("+") + txt);
        m.stampStatic(-g, nets[0], nets[1], std::string("-") + txt);
        m.stampStatic(-g, nets[1], nets[0], std::string("-") + txt);
        m.stampStatic(+g, nets[1], nets[1], std::string("+") + txt);
    }
};

struct Capacitor : Component<2, 1>
{
    double c;
    double stateVar;
    double voltage;

    Capacitor(double c, int l0, int l1) : c(c)
    {
        pinLoc[0] = l0;
        pinLoc[1] = l1;

        stateVar = 0;
        voltage = 0;
    }

    void stamp(MNASystem & m)
    {
        char buf[16];
        formatUnitValue(buf, c, "F");

        // we can use a trick here, to get the capacitor to
        // work on it's own line with direct trapezoidal:
        //
        // | -g*t  +g*t  +t | v+
        // | +g*t  -g*t  -t | v-
        // | +2*g  -2*g  -1 | state
        //
        // the logic with this is that for constant timestep:
        //
        //  i1 = g*v1 - s0   , s0 = g*v0 + i0
        //  s1 = 2*g*v1 - s0 <-> s0 = 2*g*v1 - s1
        //
        // then if we substitute back:
        //  i1 = g*v1 - (2*g*v1 - s1)
        //     = s1 - g*v1
        //
        // this way we just need to copy the new state to the
        // next timestep and there's no actual integration needed
        //
        // the "half time-step" error here means that our state
        // is 2*c*v - i/t but we fix this for display in update
        // and correct the current-part on time-step changes

        // trapezoidal needs another factor of two for the g
        // since c*(v1 - v0) = (i1 + i0)/(2*t), where t = 1/T
        double g = 2*c;

        m.stampTimed(+1, nets[0], nets[2], "+t");
        m.stampTimed(-1, nets[1], nets[2], "-t");

        m.stampTimed(-g, nets[0], nets[0], std::string("-t*") + buf);
        m.stampTimed(+g, nets[0], nets[1], std::string("+t*") + buf);
        m.stampTimed(+g, nets[1], nets[0], std::string("+t*") + buf);
        m.stampTimed(-g, nets[1], nets[1], std::string("-t*") + buf);

        m.stampStatic(+2*g, nets[2], nets[0], std::string("+2*") + buf);
        m.stampStatic(-2*g, nets[2], nets[1], std::string("-2*") + buf);

        m.stampStatic(-1, nets[2], nets[2], "-1");

        // see the comment about v:C[%d] below
        sprintf(buf, "q:C:%d,%d", pinLoc[0], pinLoc[1]);
        m.b[nets[2]].gdyn.push_back(&stateVar);
        m.b[nets[2]].txt = buf;

        // this isn't quite right as state stores 2*c*v - i/t
        // however, we'll fix this in updateFull() for display
        sprintf(buf, "v:C:%d,%d", pinLoc[0], pinLoc[1]);
        m.nodes[nets[2]].name = buf;
        m.nodes[nets[2]].scale = 1 / c;
    }

    void update(MNASystem & m)
    {
        stateVar = m.b[nets[2]].lu;

        // solve legit voltage from the pins
        voltage = m.b[nets[0]].lu - m.b[nets[1]].lu;

        // then we can store this for display here
        // since this value won't be used at this point
        m.b[nets[2]].lu = c*voltage;
    }

    void scaleTime(double told_per_new)
    {
        // the state is 2*c*voltage - i/t0
        // so we subtract out the voltage, scale current
        // and then add the voltage back to get new state
        //
        // note that this also works if the old rate is infinite
        // (ie. t0=0) when going from DC analysis to transient
        //
        double qq = 2*c*voltage;
        stateVar = qq + (stateVar - qq)*told_per_new;
    }
};

struct Voltage : Component<2, 1>
{
    double v;

    Voltage(double v, int l0, int l1) : v(v)
    {
        pinLoc[0] = l0;
        pinLoc[1] = l1;
    }

    void stamp(MNASystem & m)
    {
        m.stampStatic(-1, nets[0], nets[2], "-1");
        m.stampStatic(+1, nets[1], nets[2], "+1");

        m.stampStatic(+1, nets[2], nets[0], "+1");
        m.stampStatic(-1, nets[2], nets[1], "-1");

        char buf[16];
        sprintf(buf, "%+.2gV", v);
        m.b[nets[2]].g = v;
        m.b[nets[2]].txt = buf;

        sprintf(buf, "i:V(%+.2g):%d,%d", v, pinLoc[0], pinLoc[1]);
        m.nodes[nets[2]].name = buf;
        m.nodes[nets[2]].type = MNANodeInfo::tCurrent;
    }
};

// probe a differential voltage
// also forces this voltage to actually get solved :)
struct Probe : Component<2, 1>
{
    Probe(int l0, int l1)
    {
        pinLoc[0] = l0;
        pinLoc[1] = l1;
    }

    void stamp(MNASystem & m)
    {
        // vp + vn - vd = 0
        m.stampStatic(+1, nets[2], nets[0], "+1");
        m.stampStatic(-1, nets[2], nets[1], "-1");
        m.stampStatic(-1, nets[2], nets[2], "-1");

        m.nodes[nets[2]].name = "v:probe";
    }

    void update([[maybe_unused]] MNASystem & m)
    {
        // we could do output here :)
    }
};

// function voltage generator
struct Function : Component<2,1>
{
    typedef double (*FuncPtr)(double t);

    FuncPtr fn;
    double  v;

    Function(FuncPtr fn, int l0, int l1) : fn(fn)
    {
        pinLoc[0] = l0;
        pinLoc[1] = l1;

        v = fn(0);
    }

    void stamp(MNASystem & m)
    {
        // this is identical to voltage source
        // except voltage is dynanic
        m.stampStatic(-1, nets[0], nets[2], "-1");
        m.stampStatic(+1, nets[1], nets[2], "+1");

        m.stampStatic(+1, nets[2], nets[0], "+1");
        m.stampStatic(-1, nets[2], nets[1], "-1");

        char buf[16];
        m.b[nets[2]].gdyn.push_back(&v);
        sprintf(buf, "Vfn:%d,%d", pinLoc[0], pinLoc[1]);
        m.b[nets[2]].txt = buf;

        sprintf(buf, "i:Vfn:%d,%d", pinLoc[0], pinLoc[1]);
        m.nodes[nets[2]].name = buf;
        m.nodes[nets[2]].type = MNANodeInfo::tCurrent;
    }

    void update(MNASystem & m)
    {
        v = fn(m.time);
    }

};

// POD-struct for PN-junction data, for diodes and BJTs
//
struct JunctionPN
{
    // variables
    double geq, ieq, veq;

    // parameters
    double is, nvt, rnvt, vcrit;
};

void initJunctionPN(JunctionPN & pn, double is, double n)
{
    pn.is = is;
    pn.nvt = n * vThermal;
    pn.rnvt = 1 / pn.nvt;
    pn.vcrit = pn.nvt * log(pn.nvt / (pn.is * sqrt(2.)));
}

// linearize junction at the specified voltage
//
// ideally we could handle series resistance here as well
// to avoid putting it on a separate node, but not sure how
// to make that work as it looks like we'd need Lambert-W then
void linearizeJunctionPN(JunctionPN & pn, double v)
{
    double e = pn.is * exp(v * pn.rnvt);
    double i = e - pn.is + gMin * v;
    double g = e * pn.rnvt + gMin;

    pn.geq = g;
    pn.ieq = v*g - i;
    pn.veq = v;
}

// returns true if junction is good enough
bool newtonJunctionPN(JunctionPN & pn, double v)
{
    double dv = v - pn.veq;
    if(fabs(dv) < vTolerance) return true;

    // check critical voltage and adjust voltage if over
    if(v > pn.vcrit)
    {
        // this formula comes from Qucs documentation
        v = pn.veq + pn.nvt*log((std::max)(pn.is, 1+dv*pn.rnvt));
    }

    linearizeJunctionPN(pn, v);

    return false;
}

struct Diode : Component<2, 2>
{
    JunctionPN  pn;

    // should make these parameters
    double rs;

    // l0 -->|-- l1 -- parameters default to approx 1N4148
    Diode(int l0, int l1,
        double rs = 10., double is = 35e-12, double n = 1.24)
        : rs(rs)
    {
        pinLoc[0] = l0;
        pinLoc[1] = l1;

        initJunctionPN(pn, is, n);

        // FIXME: move init to some restart routine?

        // initial condition v = 0
        linearizeJunctionPN(pn, 0);
    }

    bool newton(MNASystem & m)
    {
        return newtonJunctionPN(pn, m.b[nets[2]].lu);
    }

    void stamp(MNASystem & m)
    {
        // Diode could be built with 3 extra nodes:
        //
        // |  .  .    .       . +1 | V+
        // |  .  .    .       . -1 | V-
        // |  .  .  grs    -grs -1 | v:D
        // |  .  . -grs grs+geq  . | v:pn = ieq
        // | -1 +1   +1       .  . | i:pn
        //
        // Here grs is the 1/rs series conductance.
        //
        // This gives us the junction voltage (v:pn) and
        // current (i:pn) and the composite voltage (v:D).
        //
        // The i:pn row is an ideal transformer connecting
        // the floating diode to the ground-referenced v:D
        // where we connect the series resistance to v:pn
        // that solves the diode equation with Newton.
        //
        // We can then add the 3rd row to the bottom 2 with
        // multipliers 1 and -rs = -1/grs and drop it:
        //
        // |  .  .   . +1 | V+
        // |  .  .   . -1 | V-
        // |  .  . geq -1 | v:pn = ieq
        // | -1 +1  +1 rs | i:pn
        //
        // Note that only the v:pn row here is non-linear.
        //
        // We could even do away without the separate row for
        // the current, which would lead to the following:
        //
        // | +grs -grs     -grs |
        // | -grs +grs     +grs |
        // | -grs +grs +grs+geq | = ieq
        //
        // In practice we keep the current row since it's
        // nice to have it as an output anyway.
        //
        m.stampStatic(-1, nets[3], nets[0], "-1");
        m.stampStatic(+1, nets[3], nets[1], "+1");
        m.stampStatic(+1, nets[3], nets[2], "+1");

        m.stampStatic(+1, nets[0], nets[3], "+1");
        m.stampStatic(-1, nets[1], nets[3], "-1");
        m.stampStatic(-1, nets[2], nets[3], "-1");

        m.stampStatic(rs, nets[3], nets[3], "rs:pn");

        m.A[nets[2]][nets[2]].gdyn.push_back(&pn.geq);
        m.A[nets[2]][nets[2]].txt = "gm:D";
        m.b[nets[2]].gdyn.push_back(&pn.ieq);

        char buf[16];
        sprintf(buf, "i0:D:%d,%d", pinLoc[0], pinLoc[1]);
        m.b[nets[2]].txt = buf;

        sprintf(buf, "v:D:%d,%d", pinLoc[0], pinLoc[1]);
        m.nodes[nets[2]].name = buf;

        sprintf(buf, "i:D:%d,%d", pinLoc[0], pinLoc[1]);
        m.nodes[nets[3]].name = buf;
        m.nodes[nets[3]].type = MNANodeInfo::tCurrent;

    }
};

struct BJT : Component<3, 4>
{
    // emitter and collector junctions
    JunctionPN  pnC, pnE;

    // forward and reverse alpha
    double af, ar, rsbc, rsbe;

    bool pnp;

    BJT(int b, int c, int e, bool pnp = false) : pnp(pnp)
    {
        pinLoc[0] = b;
        pinLoc[1] = c;
        pinLoc[2] = e;

        // this attempts a 2n3904-style small-signal
        // transistor, although the values are a bit
        // arbitrarily set to "something reasonable"

        // forward and reverse beta
        double bf = 200;
        double br = 20;

        // forward and reverse alpha
        af = bf / (1 + bf);
        ar = br / (1 + br);

        // these are just rb+re and rb+rc
        // this is not necessarily the best way to
        // do anything, but having junction series
        // resistances helps handle degenerate cases
        rsbc = 5.8376+0.0001;
        rsbe = 5.8376+2.65711;

        //
        // the basic rule is that:
        //  af * ise = ar * isc = is
        //
        // FIXME: with non-equal ideality factors
        // we can get non-sensical results, why?
        //
        double is = 6.734e-15;
        double n = 1.24;
        initJunctionPN(pnE, is / af, n);
        initJunctionPN(pnC, is / ar, n);

        linearizeJunctionPN(pnE, 0);
        linearizeJunctionPN(pnC, 0);
    }

    bool newton(MNASystem & m)
    {
        return newtonJunctionPN(pnC, m.b[nets[3]].lu)
             & newtonJunctionPN(pnE, m.b[nets[4]].lu);
    }

    void stamp(MNASystem & m)
    {
        // The basic idea here is the same as with diodes
        // except we do it once for each junction.
        //
        // With the transfer currents sourced from the
        // diode currents, NPN then looks like this:
        //
        // 0 |  .  .  .  .  . 1-ar 1-af | vB
        // 1 |  .  .  .  .  .   -1  +af | vC
        // 2 |  .  .  .  .  .  +ar   -1 | vE
        // 3 |  .  .  . gc  .   -1    . | v:Qbc  = ic
        // 4 |  .  .  .  . ge    .   -1 | v:Qbe  = ie
        // 5 | -1 +1  . +1  . rsbc    . | i:Qbc
        // 6 | -1  . +1  . +1    . rsbe | i:Qbe
        //     ------------------------
        //      0  1  2  3  4    5    6
        //
        // For PNP version, we simply flip the junctions
        // by changing signs of (3,5),(5,3) and (4,6),(6,4).
        //
        // Also just like diodes, we have junction series
        // resistances, rather than terminal resistances.
        //
        // This works just as well, but should be kept
        // in mind when fitting particular transistors.
        //

        // diode currents to external base
        m.stampStatic(1-ar, nets[0], nets[5], "1-ar");
        m.stampStatic(1-af, nets[0], nets[6], "1-af");

        // diode currents to external collector and emitter
        m.stampStatic(-1, nets[1], nets[5], "-1");
        m.stampStatic(-1, nets[2], nets[6], "-1");

        // series resistances
        m.stampStatic(rsbc, nets[5], nets[5], "rsbc");
        m.stampStatic(rsbe, nets[6], nets[6], "rsbe");

        // current - junction connections
        // for the PNP case we flip the signs of these
        // to flip the diode junctions wrt. the above
        if(pnp)
        {
            m.stampStatic(-1, nets[5], nets[3], "-1");
            m.stampStatic(+1, nets[3], nets[5], "+1");

            m.stampStatic(-1, nets[6], nets[4], "-1");
            m.stampStatic(+1, nets[4], nets[6], "+1");

        }
        else
        {
            m.stampStatic(+1, nets[5], nets[3], "+1");
            m.stampStatic(-1, nets[3], nets[5], "-1");

            m.stampStatic(+1, nets[6], nets[4], "+1");
            m.stampStatic(-1, nets[4], nets[6], "-1");
        }

        // external voltages to collector current
        m.stampStatic(-1, nets[5], nets[0], "-1");
        m.stampStatic(+1, nets[5], nets[1], "+1");

        // external voltages to emitter current
        m.stampStatic(-1, nets[6], nets[0], "-1");
        m.stampStatic(+1, nets[6], nets[2], "+1");

        // source transfer currents to external pins
        m.stampStatic(+ar, nets[2], nets[5], "+ar");
        m.stampStatic(+af, nets[1], nets[6], "+af");

        char buf[16];

        // dynamic variables
        m.A[nets[3]][nets[3]].gdyn.push_back(&pnC.geq);
        m.A[nets[3]][nets[3]].txt = "gm:Qbc";
        m.b[nets[3]].gdyn.push_back(&pnC.ieq);
        sprintf(buf, "i0:Q:%d,%d,%d:cb", pinLoc[0], pinLoc[1], pinLoc[2]);
        m.b[nets[3]].txt = buf;

        m.A[nets[4]][nets[4]].gdyn.push_back(&pnE.geq);
        m.A[nets[4]][nets[4]].txt = "gm:Qbe";
        m.b[nets[4]].gdyn.push_back(&pnE.ieq);
        sprintf(buf, "i0:Q:%d,%d,%d:eb", pinLoc[0], pinLoc[1], pinLoc[2]);
        m.b[nets[4]].txt = buf;

        sprintf(buf, "v:Q:%d,%d,%d:%s",
            pinLoc[0], pinLoc[1], pinLoc[2], pnp ? "cb" : "bc");
        m.nodes[nets[3]].name = buf;

        sprintf(buf, "v:Q:%d,%d,%d:%s",
            pinLoc[0], pinLoc[1], pinLoc[2], pnp ? "eb" : "be");
        m.nodes[nets[4]].name = buf;

        sprintf(buf, "i:Q:%d,%d,%d:bc", pinLoc[0], pinLoc[1], pinLoc[2]);
        m.nodes[nets[5]].name = buf;
        m.nodes[nets[5]].type = MNANodeInfo::tCurrent;
        m.nodes[nets[5]].scale = 1 - ar;

        sprintf(buf, "i:Q:%d,%d,%d:be", pinLoc[0], pinLoc[1], pinLoc[2]);
        m.nodes[nets[6]].name = buf;
        m.nodes[nets[6]].type = MNANodeInfo::tCurrent;
        m.nodes[nets[6]].scale = 1 - af;
    }
};


struct NetList
{
    typedef std::vector<IComponent*> ComponentList;

    NetList(int nodes) : nets(nodes), states(0)
    {

    }

    void addComponent(IComponent * c)
    {
        // this is a bit "temporary" for now
        c->setupNets(nets, states, c->getPinLocs());
        components.push_back(c);
    }

    void buildSystem()
    {
        system.setSize(nets);
        for(size_t i = 0; i < components.size(); ++i)
        {
            components[i]->stamp(system);
        }
        printf("Prepare for DC analysis..\n");
        setStepScale(0);
        tStep = 0;
    }

    void dumpMatrix()
    {
        std::vector<int> maxWidth(nets);

        for(int i = 0; i < nets; ++i) maxWidth[i] = 1;
        int nnMax = 1;

        for(int i = 0; i < nets; ++i)
        {
            nnMax = std::max(nnMax, (int)system.nodes[i].name.size());
            for(int j = 0; j < nets; ++j)
            {
                maxWidth[j] = std::max(maxWidth[j],
                    (int)system.A[i][j].txt.size());
            }
        }


        char buf[1024];
        for(int i = 0; i < nets; ++i)
        {
            int off = sprintf(buf, "%2d: | ", i);
            for(int j = 0; j < nets; ++j)
            {
                off += sprintf(buf+off,
                    " %*s ", maxWidth[j],
                    system.A[i][j].txt.size()
                    ? system.A[i][j].txt.c_str()
                    : ((system.A[i][j].lu==0) ? "." : "#"));
            }
            sprintf(buf+off, " | %-*s = %s",
                nnMax, system.nodes[i].name.c_str(),
                system.b[i].txt.size()
                ? system.b[i].txt.c_str() : (!i ? "ground" : "0"));

            puts(buf);
        }
    }

    void setTimeStep(double tStepSize)
    {
        for(size_t i = 0; i < components.size(); ++i)
        {
            components[i]->scaleTime(tStep / tStepSize);
        }

        tStep = tStepSize;
        double stepScale = 1. / tStep;
        printf("timeStep changed to %.2g (%.2g Hz)\n", tStep, stepScale);
        setStepScale(stepScale);
    }

    void simulateTick()
    {
        size_t iter;
        for(iter = 0; iter < maxIter; ++iter)
        {
            // restore matrix state and add dynamic values
            updatePre();
            luFactor();
            luForward();
            luSolve();

            if(newton()) break;
        }

        system.time += tStep;

        update();

        printf(" %02.4f |", system.time);
        int fillPost = 0;
        for(int i = 1; i < nets; ++i)
        {
            printf("\t%+.4e", system.b[i].lu * system.nodes[i].scale);
            for(int j = 1; j < nets; ++j)
            {
                if(system.A[i][j].lu != 0) ++fillPost;
            }
        }
        printf("\t %lu iters, LU density: %.1f%%\n",
            iter, 100 * fillPost / ((nets-1.f)*(nets-1.f)));
    }

    void printHeaders()
    {
        printf("\n  time: |  ");
        for(int i = 1; i < nets; ++i)
        {
            printf("%16s", system.nodes[i].name.c_str());
        }
        printf("\n\n");
    }

    // plotting and such would want to use this
    const MNASystem & getMNA() { return system; }

protected:
    double  tStep;

    int nets, states;
    ComponentList   components;

    MNASystem       system;

    void update()
    {
        for(size_t i = 0; i < components.size(); ++i)
        {
            components[i]->update(system);
        }
    }

    // return true if we're done
    bool newton()
    {
        bool done = 1;
        for(size_t i = 0; i < components.size(); ++i)
        {
            done &= components[i]->newton(system);
        }
        return done;
    }

    void initLU(double stepScale)
    {
        for(int i = 0; i < nets; ++i)
        {
            system.b[i].initLU(stepScale);
            for(int j = 0; j < nets; ++j)
            {
                system.A[i][j].initLU(stepScale);
            }
        }
    }

    void setStepScale(double stepScale)
    {
        // initialize matrix for LU and save it to cache
        initLU(stepScale);

        int fill = 0;
        for(int i = 1; i < nets; ++i)
        {
            for(int j = 1; j < nets; ++j)
            {
                if(system.A[i][j].prelu != 0
                    || system.A[i][j].gdyn.size()) ++fill;
            }
        }
        printf("MNA density %.1f%%\n", 100 * fill / ((nets-1.)*(nets-1.)));
    }

    void updatePre()
    {
        for(int i = 0; i < nets; ++i)
        {
            system.b[i].updatePre();
            for(int j = 0; j < nets; ++j)
            {
                system.A[i][j].updatePre();
            }
        }
    }

    void luFactor()
    {
        int p;
        for(p = 1; p < nets; ++p)
        {
            // FIND PIVOT
            {
                int pr = p;
                for(int r = p; r < nets; ++r)
                {
                    if(fabs(system.A[r][p].lu)
                    > fabs(system.A[pr][p].lu))
                    {
                        pr = r;
                    }
                }
                // swap if necessary
                if(pr != p)
                {
                    std::swap(system.A[p], system.A[pr]);
                    std::swap(system.b[p], system.b[pr]);
                }
                #if VERBOSE_LU
                printf("pivot %d (from %d): %+.2e\n",
                p, pr, system.A[p][p].lu);
                #endif
            }
            if(0 == system.A[p][p].lu)
            {
                printf("Failed to find a pivot!!");
                return;
            }

            // take reciprocal for D entry
            system.A[p][p].lu = 1 / system.A[p][p].lu;

            // perform reduction on rows below
            for(int r = p+1; r < nets; ++r)
            {
                if(system.A[r][p].lu == 0) continue;

                system.A[r][p].lu *= system.A[p][p].lu;
                for(int c = p+1; c < nets; ++c)
                {
                    if(system.A[p][c].lu == 0) continue;

                    system.A[r][c].lu -=
                    system.A[p][c].lu * system.A[r][p].lu;
                }

            }
        }
    }

    // this does forward substitution for the solution vector
    int luForward()
    {
        int p;
        for(p = 1; p < nets; ++p)
        {
            // perform reduction on rows below
            for(int r = p+1; r < nets; ++r)
            {
                if(system.A[r][p].lu == 0) continue;
                if(system.b[p].lu == 0) continue;

                system.b[r].lu -= system.b[p].lu * system.A[r][p].lu;
            }
        }
        return p;
    }

    // solves nodes backwards from limit-1
    // if solveAll is true, solves all the nodes
    // otherwise if solveNoIter is true, solves until !wantUpdate
    // if both flags are false, solves until !wantIter
    int luSolve()
    {
        for(int r = nets; --r;)
        {
            //printf("solve node %d\n", r);
            for(int s = r+1; s < nets; ++s)
            {
                system.b[r].lu -= system.b[s].lu * system.A[r][s].lu;
            }

            system.b[r].lu *= system.A[r][r].lu;
        }
        return 1;
    }
};


double fnGen(double t)
{
    if(t < 0.0001) return 0;
    return (fmod(2000*t, 1) > (.5+.4*sin(2*acos(-1)*100*t))) ? .25 : -.25;
}

// int main()
// {
//     /*

//         BJT test circuit

//         this is kinda good at catching problems :)

//     */
//     NetList * net = new NetList(8);

//     // node 1 is +12V
//     net->addComponent(new Voltage(+12, 1, 0));

//     // bias for base
//     net->addComponent(new Resistor(100e3, 2, 1));
//     net->addComponent(new Resistor(11e3, 2, 0));

//     // input cap and function
//     net->addComponent(new Capacitor(.1e-6, 2, 5));
//     net->addComponent(new Resistor(1e3, 5, 6));
//     net->addComponent(new Function(fnGen, 6, 0));

//     // collector resistance
//     net->addComponent(new Resistor(10e3, 1, 3));

//     // emitter resistance
//     net->addComponent(new Resistor(680, 4, 0));

//     // bjt on b:2,c:3,e:4
//     net->addComponent(new BJT(2,3,4));

//     // output emitter follower (direct rail collector?)
//     net->addComponent(new BJT(3, 1, 7));
//     net->addComponent(new Resistor(100e3, 7, 0));

//     net->addComponent(new Probe(7, 0));

//     net->buildSystem();

//     // get DC solution (optional)
//     net->simulateTick();
//     net->setTimeStep(1e-4);
//     net->printHeaders();

//     for(int i = 0; i < 25; ++i)
//     {
//         net->simulateTick();
//     }
//     net->printHeaders();

//     net->dumpMatrix();
// };