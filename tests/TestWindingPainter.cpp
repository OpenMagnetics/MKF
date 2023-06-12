#include "WindingPainter.h"
#include "json.hpp"
#include "TestingUtils.h"

#include <UnitTest++.h>

#include <fstream>
#include <svg.hpp>
#include <string>


SUITE(WindingPainter) {
    std::string filePath = __FILE__;
    auto outputFilePath = filePath.substr(0, filePath.rfind("/")).append("/../output/");

    TEST(Test_Painter_Pq_Core_Distributed_Gap) {
        std::vector<uint64_t> numberTurns = {42};
        std::vector<uint64_t> numberParallels = {3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Distributed_Gap.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        auto svg = windingPainter.paint_core(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 1);
        CHECK(svg->get_children<SVG::Polygon>().size() == 4);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Distributed_Gap_Many) {
        std::vector<uint64_t> numberTurns = {42};
        std::vector<uint64_t> numberParallels = {3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.001, 9);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Distributed_Gap_Many.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        auto svg = windingPainter.paint_core(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 1);
        CHECK(svg->get_children<SVG::Polygon>().size() == 10);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap) {
        std::vector<uint64_t> numberTurns = {42};
        std::vector<uint64_t> numberParallels = {3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.003);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        auto svg = windingPainter.paint_core(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 1);
        CHECK(svg->get_children<SVG::Polygon>().size() == 2);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_U_Core_Distributed_Gap) {
        std::vector<uint64_t> numberTurns = {42};
        std::vector<uint64_t> numberParallels = {3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "U 10/8/3";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.001, 3);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_U_Core_Distributed_Gap.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        auto svg = windingPainter.paint_core(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 1);
        CHECK(svg->get_children<SVG::Polygon>().size() == 4);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_U_Core_Grinded_Gap) {
        std::vector<uint64_t> numberTurns = {42};
        std::vector<uint64_t> numberParallels = {3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "U 10/8/3";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.003);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_U_Core_Grinded_Gap.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        auto svg = windingPainter.paint_core(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 1);
        CHECK(svg->get_children<SVG::Polygon>().size() == 2);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Bobbin) {
        std::vector<uint64_t> numberTurns = {42};
        std::vector<uint64_t> numberParallels = {3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        auto svg = windingPainter.paint_bobbin(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 2);
        CHECK(svg->get_children<SVG::Polygon>().size() == 5);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Section) {
        std::vector<uint64_t> numberTurns = {42};
        std::vector<uint64_t> numberParallels = {3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Sections.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 2);
        CHECK(svg->get_children<SVG::Polygon>().size() == 6);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Bobbin_And_Section) {
        std::vector<uint64_t> numberTurns = {42};
        std::vector<uint64_t> numberParallels = {3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_And_Section.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Bobbin_And_Sections) {
        std::vector<uint64_t> numberTurns = {42, 42};
        std::vector<uint64_t> numberParallels = {3, 3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_And_Sections.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 9);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Epx_Core_Grinded_Gap) {
        std::vector<uint64_t> numberTurns = {42, 42};
        std::vector<uint64_t> numberParallels = {3, 3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "EPX 9/9";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Epx_Core_Grinded_Gap.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Epx_Core_Spacer_Gap) {
        std::vector<uint64_t> numberTurns = {42, 42};
        std::vector<uint64_t> numberParallels = {3, 3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "EPX 9/9";
        std::string coreMaterial = "3C97";
        auto gapping = json::array();
        auto basicSpacerGap = json();
        basicSpacerGap["type"] = "additive";
        basicSpacerGap["length"] = 0.0003;
        gapping.push_back(basicSpacerGap);
        gapping.push_back(basicSpacerGap);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Epx_Core_Spacer_Gap.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_P_Core_Grinded_Gap) {
        std::vector<uint64_t> numberTurns = {42, 42};
        std::vector<uint64_t> numberParallels = {3, 3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "P 3.3/2.6";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_P_Core_Grinded_Gap.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }


    TEST(Test_Painter_U80_Core_Grinded_Gap) {
        std::vector<uint64_t> numberTurns = {42, 42};
        std::vector<uint64_t> numberParallels = {3, 3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "U 80/65/32";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_U80_Core_Grinded_Gap.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Ep_Core_Grinded_Gap) {
        std::vector<uint64_t> numberTurns = {42, 42};
        std::vector<uint64_t> numberParallels = {3, 3};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "EP 10";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Ep_Core_Grinded_Gap.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_All_Cores) {

        for (std::string& shapeName : OpenMagnetics::get_shape_names()){
            if (shapeName.contains("PQI") || shapeName.contains("R ") || shapeName.contains("T ") || shapeName.contains("UI ")) {
                continue;
            }
            std::vector<uint64_t> numberTurns = {42, 42};
            std::vector<uint64_t> numberParallels = {3, 3};
            uint64_t interleavingLevel = 2;
            uint64_t numberStacks = 1;
            std::string coreShape = shapeName;
            std::string coreMaterial = "3C97";
            auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

            auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
            auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

            std::replace( shapeName.begin(), shapeName.end(), '.', '_');
            std::replace( shapeName.begin(), shapeName.end(), '/', '_');
            auto outFile = outputFilePath;
            outFile.append("Test_Painter_Core_" + shapeName + ".svg");
            OpenMagnetics::WindingPainter windingPainter(outFile);
            OpenMagnetics::Magnetic magnetic;
            magnetic.set_core(core);
            magnetic.set_winding(winding);

            windingPainter.paint_core(magnetic);
            windingPainter.paint_bobbin(magnetic);
            auto svg = windingPainter.paint_winding_sections(magnetic);

            CHECK(svg->get_children<SVG::Group>().size() == 3);
            CHECK(svg->get_children<SVG::Polygon>().size() == 7);
            CHECK(std::filesystem::exists(outFile));
        }
    }


    TEST(Test_Painter_Pq_Core_Grinded_Gap_Layers_No_Interleaving) {
        std::vector<uint64_t> numberTurns = {42, 42};
        std::vector<uint64_t> numberParallels = {1, 1};
        uint64_t interleavingLevel = 1;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Layers_No_Interleaving.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_layers(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap_Turns_No_Interleaving) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {2, 2};
        uint64_t interleavingLevel = 1;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_No_Interleaving.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 3);
        CHECK(svg->get_children<SVG::Circle>().size() == 140);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {4, 4};
        uint64_t interleavingLevel = 3;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 3);
        CHECK(svg->get_children<SVG::Circle>().size() == 280);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Top_Alignment) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {4, 4};
        uint64_t interleavingLevel = 3;
        uint64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Top_Alignment.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 3);
        CHECK(svg->get_children<SVG::Circle>().size() == 280);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Bottom_Alignment) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {4, 4};
        uint64_t interleavingLevel = 3;
        uint64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Bottom_Alignment.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 3);
        CHECK(svg->get_children<SVG::Circle>().size() == 280);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Spread_Alignment) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {4, 4};
        uint64_t interleavingLevel = 3;
        uint64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Grinded_Gap_Turns_Interleaving_Spread_Alignment.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 3);
        CHECK(svg->get_children<SVG::Circle>().size() == 280);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Vertical_Sections) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {1, 1};
        uint64_t interleavingLevel = 3;
        uint64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 9);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Vertical_Sections_Vectical_Layers) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {3, 3};
        uint64_t interleavingLevel = 3;
        uint64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Vectical_Layers.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_layers(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 45);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Vertical_Sections_Horizontal_Layers) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {1, 1};
        uint64_t interleavingLevel = 3;
        uint64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_layers(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 9);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Vertical_Sections_Horizontal_Layers_Spread_Turns) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {4, 4};
        uint64_t interleavingLevel = 3;
        uint64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Spread_Turns.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 3);
        CHECK(svg->get_children<SVG::Circle>().size() == 280);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Vertical_Sections_Horizontal_Layers_Inner_Turns) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {4, 4};
        uint64_t interleavingLevel = 3;
        uint64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Inner_Turns.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 3);
        CHECK(svg->get_children<SVG::Circle>().size() == 280);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Vertical_Sections_Horizontal_Layers_Outer_Turns) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {4, 4};
        uint64_t interleavingLevel = 3;
        uint64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Outer_Turns.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 3);
        CHECK(svg->get_children<SVG::Circle>().size() == 280);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Vertical_Sections_Horizontal_Layers_Centered_Turns) {
        std::vector<uint64_t> numberTurns = {35, 35};
        std::vector<uint64_t> numberParallels = {4, 4};
        uint64_t interleavingLevel = 3;
        uint64_t numberStacks = 1;
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        std::string coreShape = "PQ 35/30";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Vertical_Sections_Horizontal_Layers_Centered_Turns.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 3);
        CHECK(svg->get_children<SVG::Circle>().size() == 280);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Foil_Centered) {
        std::vector<uint64_t> numberTurns = {4};
        std::vector<uint64_t> numberParallels = {1};
        uint64_t interleavingLevel = 1;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_outer_height(0.014);
        wire.set_nominal_value_outer_width(0.0002);
        wire.set_type("foil");
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Foil_Centered.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Foil_Top) {
        std::vector<uint64_t> numberTurns = {4};
        std::vector<uint64_t> numberParallels = {1};
        uint64_t interleavingLevel = 1;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        OpenMagnetics::WireWrapper wire;
        wire.set_nominal_value_outer_height(0.014);
        wire.set_nominal_value_outer_width(0.0002);
        wire.set_type("foil");
        auto wires = std::vector<OpenMagnetics::WireWrapper>({wire});

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment, wires);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Foil_Top.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Centered) {
        std::vector<uint64_t> numberTurns = {23, 23};
        std::vector<uint64_t> numberParallels = {2, 2};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Centered.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Centered) {
        std::vector<uint64_t> numberTurns = {23, 23};
        std::vector<uint64_t> numberParallels = {2, 2};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Centered.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Top) {
        std::vector<uint64_t> numberTurns = {23, 23};
        std::vector<uint64_t> numberParallels = {2, 2};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Top.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Inner) {
        std::vector<uint64_t> numberTurns = {23, 23};
        std::vector<uint64_t> numberParallels = {2, 2};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::INNER_OR_TOP;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Inner.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Outer) {
        std::vector<uint64_t> numberTurns = {23, 23};
        std::vector<uint64_t> numberParallels = {2, 2};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Outer.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Bottom) {
        std::vector<uint64_t> numberTurns = {23, 23};
        std::vector<uint64_t> numberParallels = {2, 2};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::OUTER_OR_BOTTOM;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Bottom.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Spread) {
        std::vector<uint64_t> numberTurns = {23, 23};
        std::vector<uint64_t> numberParallels = {2, 2};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Spread.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Spread_Two_Sections) {
        std::vector<uint64_t> numberTurns = {23, 23};
        std::vector<uint64_t> numberParallels = {2, 2};
        uint64_t interleavingLevel = 1;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Spread_Two_Sections.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 5);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Vertical_Spread_One_Section) {
        std::vector<uint64_t> numberTurns = {23};
        std::vector<uint64_t> numberParallels = {2};
        uint64_t interleavingLevel = 1;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Vertical_Spread_One_Section.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 4);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Spread) {
        std::vector<uint64_t> numberTurns = {23, 23};
        std::vector<uint64_t> numberParallels = {2, 2};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 7);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_Two_Sections) {
        std::vector<uint64_t> numberTurns = {23, 23};
        std::vector<uint64_t> numberParallels = {2, 2};
        uint64_t interleavingLevel = 1;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_Two_Sections.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 5);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_One_Section) {
        std::vector<uint64_t> numberTurns = {23};
        std::vector<uint64_t> numberParallels = {2};
        uint64_t interleavingLevel = 1;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 26/25";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_grinded_gap(0.0001);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::SPREAD;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Delimit_Coil_Sections_Horizontal_Spread_One_Section.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);
        // auto svg = windingPainter.paint_winding_layers(magnetic);
        svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 4);
        CHECK(svg->get_children<SVG::Polygon>().size() == 4);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Bobbin_Vertical_Sections_And_Insulation) {
        std::vector<uint64_t> numberTurns = {42, 42};
        std::vector<uint64_t> numberParallels = {3, 3};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        double voltagePeakToPeak = 20000;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operation_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        winding.set_inputs(inputs);
        winding.wind();
        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Vertical_Sections_And_Insulation.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 12);
        CHECK(std::filesystem::exists(outFile));
    }


    TEST(Test_Painter_Pq_Core_Bobbin_Horizontal_Sections_And_Insulation) {
        std::vector<uint64_t> numberTurns = {42, 42};
        std::vector<uint64_t> numberParallels = {3, 3};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);
        OpenMagnetics::WindingOrientation sectionOrientation = OpenMagnetics::WindingOrientation::HORIZONTAL;
        OpenMagnetics::WindingOrientation layersOrientation = OpenMagnetics::WindingOrientation::VERTICAL;
        OpenMagnetics::CoilAlignment sectionsAlignment = OpenMagnetics::CoilAlignment::CENTERED;
        OpenMagnetics::CoilAlignment turnsAlignment = OpenMagnetics::CoilAlignment::CENTERED;

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel, sectionOrientation, layersOrientation, turnsAlignment, sectionsAlignment);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        double voltagePeakToPeak = 20000;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operation_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        winding.set_inputs(inputs);
        winding.wind();
        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Horizontal_Sections_And_Insulation.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_sections(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 12);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Bobbin_Layers_And_Insulation) {
        std::vector<uint64_t> numberTurns = {42, 42};
        std::vector<uint64_t> numberParallels = {3, 3};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        double voltagePeakToPeak = 20000;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operation_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        winding.set_inputs(inputs);
        winding.wind();
        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Layers_And_Insulation.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_layers(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 19);
        CHECK(std::filesystem::exists(outFile));
    }

    TEST(Test_Painter_Pq_Core_Bobbin_Turns_And_Insulation) {
        std::vector<uint64_t> numberTurns = {42, 42};
        std::vector<uint64_t> numberParallels = {3, 3};
        std::vector<double> turnsRatios = {double(numberTurns[0]) / numberTurns[1]};
        uint64_t interleavingLevel = 2;
        uint64_t numberStacks = 1;
        std::string coreShape = "PQ 40/40";
        std::string coreMaterial = "3C97";
        auto gapping = OpenMagneticsTesting::get_distributed_gap(0.003, 3);

        auto winding = OpenMagneticsTesting::get_quick_winding(numberTurns, numberParallels, coreShape, interleavingLevel);
        auto core = OpenMagneticsTesting::get_core(coreShape, gapping, numberStacks, coreMaterial);

        double voltagePeakToPeak = 20000;
        auto inputs = OpenMagnetics::InputsWrapper::create_quick_operation_point(125000, 0.001, 25, OpenMagnetics::WaveformLabel::SINUSOIDAL, voltagePeakToPeak, 0.5, 0, turnsRatios);
        winding.set_inputs(inputs);
        winding.wind();
        winding.delimit_and_compact();

        auto outFile = outputFilePath;
        outFile.append("Test_Painter_Pq_Core_Bobbin_Turns_And_Insulation.svg");
        OpenMagnetics::WindingPainter windingPainter(outFile);
        OpenMagnetics::Magnetic magnetic;
        magnetic.set_core(core);
        magnetic.set_winding(winding);

        windingPainter.paint_core(magnetic);
        windingPainter.paint_bobbin(magnetic);
        auto svg = windingPainter.paint_winding_turns(magnetic);

        CHECK(svg->get_children<SVG::Group>().size() == 3);
        CHECK(svg->get_children<SVG::Polygon>().size() == 11);
        CHECK(std::filesystem::exists(outFile));
    }
}