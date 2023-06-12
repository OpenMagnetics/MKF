#pragma once

#include "CoreWrapper.h"
#include "WindingWrapper.h"
#include "Utils.h"
#include "svg.hpp"
#include <MAS.hpp>

namespace OpenMagnetics {

class WindingPainter : public Winding {
    private:
        std::filesystem::path _filepath = ".";
        std::filesystem::path _filename = "_.svg";
    public:
        SVG::SVG* _root;
        WindingPainter(std::filesystem::path filepath){
            _filename = filepath.filename();
            _filepath = filepath.remove_filename();
            _root = new SVG::SVG();
            _root->style(".ferrite").set_attr("fill", "#7b7c7d");
            _root->style(".bobbin").set_attr("fill", "#1b1b1b");
        };
        virtual ~WindingPainter() = default;

    SVG::SVG* paint_core(Magnetic magnetic);
    SVG::SVG* paint_two_piece_set_core(CoreWrapper core);
    SVG::SVG* paint_bobbin(Magnetic magnetic);
    SVG::SVG* paint_two_piece_set_bobbin(Magnetic magnetic);
    SVG::SVG* paint_winding_sections(Magnetic magnetic);
    SVG::SVG* paint_two_piece_set_winding_sections(Magnetic magnetic);
    SVG::SVG* paint_winding_layers(Magnetic magnetic);
    SVG::SVG* paint_two_piece_set_winding_layers(Magnetic magnetic);
    SVG::SVG* paint_winding_turns(Magnetic magnetic);
    SVG::SVG* paint_two_piece_set_winding_turns(Magnetic magnetic);

};
} // namespace OpenMagnetics
