#include "support/CircuitSimulator.cpp"

#include <UnitTest++.h>
#include <filesystem>
#include <cfloat>
#include <limits>
#include <fstream>
#include <iostream>
#include <magic_enum.hpp>
#include <vector>
#include <typeinfo>

SUITE(CircuitSimulator) {

    // TEST(Test_Circuit_Simulator) {
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
    // }

    // TEST(Test_Circuit_Simulator_R) {
    //     NetList * net = new NetList(4);

    //     // node 1 is +12V
    //     net->addComponent(new Voltage(+12, 1, 0));

    //     // bias for base
    //     net->addComponent(new Resistor(100e3, 1, 0));

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
    // }

    TEST(Test_Circuit_Simulator_R) {
        NetList * net = new NetList(3);

        // node 1 is +12V
        net->addComponent(new Voltage(+12, 2, 0));

        // bias for base
        net->addComponent(new Resistor(100e3, 1, 0));
        net->addComponent(new Resistor(200e3, 2, 1));

        net->buildSystem();

        // get DC solution (optional)
        net->simulateTick();
        net->printHeaders();

        net->dumpMatrix();
    }

    double sinVoltage(double time) {
        double angle = 2 * std::numbers::pi * time;
        return sin(angle);
    }

    TEST(Test_Circuit_Simulator_C) {
        NetList * net = new NetList(3);
        // node 1 is +12V
        net->addComponent(new Function(sinVoltage, 2, 0));

        // bias for base
        net->addComponent(new Capacitor(10e-12, 1, 0));
        net->addComponent(new Capacitor(10e-12, 2, 1));

        net->buildSystem();

        // get DC solution (optional)
        net->simulateTick();
        net->setTimeStep(1e-4);
        net->printHeaders();

        for(int i = 0; i < 100; ++i)
        {
            net->simulateTick();
        }
        net->printHeaders();

        net->dumpMatrix();
    }
}