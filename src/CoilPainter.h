#pragma once

#include "CoreWrapper.h"
#include "CoilWrapper.h"
#include "Utils.h"
#include "svg.hpp"
#include <MAS.hpp>

namespace OpenMagnetics {

class CoilPainter : public Coil {
    private:
        std::filesystem::path _filepath = ".";
        std::filesystem::path _filename = "_.svg";
        double _opacity = 1;
    public:
        SVG::SVG* _root;
        CoilPainter(std::filesystem::path filepath){
            _filename = filepath.filename();
            _filepath = filepath.remove_filename();
            _root = new SVG::SVG();
            _root->style(".ferrite").set_attr("fill", "#7b7c7d");
            _root->style(".bobbin").set_attr("fill", "#1b1b1b");
            _root->style(".copper").set_attr("fill", "#B87333");
        };
        virtual ~CoilPainter() = default;

    void set_opacity(double opacity) {
        _opacity = opacity;
    }
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
