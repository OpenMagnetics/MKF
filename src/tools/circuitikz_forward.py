#!/usr/bin/env python3
"""
CircuiTikz Single Switch Forward Converter Schematic Generator

This tool generates a Single Switch Forward converter circuit diagram using CircuiTikz (LaTeX)
and exports it to PDF, PNG, or SVG formats.

Usage:
    python circuitikz_forward.py [--output OUTPUT] [--format FORMAT] [--no-values]

Requirements:
    - pdflatex with circuitikz package
    - pdf2svg (for SVG output)
    - ImageMagick convert (for PNG output)
"""

import argparse
import subprocess
import tempfile
import shutil
from pathlib import Path
from dataclasses import dataclass
from typing import Optional


@dataclass
class ForwardConverterParams:
    """Parameters for the Forward converter circuit"""
    input_voltage: float = 400.0         # V (rectified mains)
    output_voltage: float = 12.0         # V
    output_current: float = 5.0          # A
    switching_frequency: float = 100e3   # Hz
    magnetizing_inductance: float = 1e-3 # H
    output_inductance: float = 10e-6     # H
    input_capacitance: float = 47e-6     # F
    output_capacitance: float = 470e-6   # F
    turns_ratio: float = 10.0            # Np/Ns
    duty_cycle: float = 0.45


def format_value_simple(value: float, unit: str) -> str:
    """Format a value with SI prefix for display (simple text, no math mode)"""
    prefixes = [
        (1e-12, 'p'), (1e-9, 'n'), (1e-6, 'u'), (1e-3, 'm'),
        (1, ''), (1e3, 'k'), (1e6, 'M'), (1e9, 'G')
    ]
    
    if unit == r"\Omega":
        unit_str = r"$\Omega$"
    else:
        unit_str = unit
    
    for scale, prefix in reversed(prefixes):
        if abs(value) >= scale:
            formatted = value / scale
            if formatted == int(formatted):
                return f"{int(formatted)}{prefix}{unit_str}"
            return f"{formatted:.1f}{prefix}{unit_str}"
    
    return f"{value}{unit_str}"


def generate_forward_circuitikz(params: Optional[ForwardConverterParams] = None, 
                                 show_values: bool = True) -> str:
    """
    Generate CircuiTikz LaTeX code for a Single Switch Forward converter
    
    Based on TI SLYU036A: The reset winding (tertiary) is used to demagnetize the core.
    During off-time, reset winding conducts through reset diode to return magnetizing 
    energy to the input.
    """
    if params is None:
        params = ForwardConverterParams()
    
    load_r = params.output_voltage / params.output_current
    
    if show_values:
        vin_label = f"{params.input_voltage:.0f}V"
        cin_label = format_value_simple(params.input_capacitance, "F")
        lout_label = format_value_simple(params.output_inductance, "H")
        cout_label = format_value_simple(params.output_capacitance, "F")
        rload_label = format_value_simple(load_r, r"\Omega")
    else:
        vin_label = "$V_{in}$"
        cin_label = "$C_{in}$"
        lout_label = "$L_{out}$"
        cout_label = "$C_{out}$"
        rload_label = "$R_{load}$"
    
    latex = r"""\documentclass[border=10pt]{standalone}
\usepackage[siunitx,RPvoltages]{circuitikz}
\usepackage{amsmath}

\begin{document}
\begin{circuitikz}[american, scale=1.0]

% === PRIMARY SIDE ===
% Input voltage source
\draw (0, 0) to[V, l=$V_{in}$] (0, 5);

% Input capacitor  
\draw (0, 5) -- (1.5, 5);
\draw (1.5, 5) to[C, l^="""
    latex += cin_label
    latex += r"""] (1.5, 0);
\draw (1.5, 0) -- (0, 0);

% Primary winding (Np turns, dot at top)
\draw (1.5, 5) -- (3, 5);
\draw (3, 5) to[L, l_=$N_p$] (3, 3);
\node[circ, fill=black, scale=0.5] at (2.7, 4.8) {};  % Dot marking

% Reset winding (Nr turns, same polarity as primary, wound in same direction)
% Connected from drain (sw_node) through diode back to Vin
\draw (3, 3) -- (3.8, 3);  % Connection point at drain
\draw (3.8, 3) -- (3.8, 5);
\draw (3.8, 5) to[L, l^=$N_r$] (3.8, 7);
\node[circ, fill=black, scale=0.5] at (4.1, 6.8) {};  % Dot marking (same side as Np)
\draw (3.8, 7) to[D, l^=$D_r$] (1.5, 7);
\draw (1.5, 7) -- (1.5, 5);

% MOSFET switch (low-side)
\draw (3, 3) to[Tnmos, n=mosfet, mirror] (3, 0);
\node[left, font=\footnotesize] at (mosfet.G) {PWM};
\draw (3, 0) -- (1.5, 0);

% Transformer core (two parallel lines)
\draw[thick] (4.3, 3.2) -- (4.3, 4.8);
\draw[thick] (4.5, 3.2) -- (4.5, 4.8);

% === SECONDARY SIDE ===
% Secondary winding (Ns turns, opposite dot position for forward operation)
\draw (5, 5) to[L, l^=$N_s$] (5, 3);
\node[circ, fill=black, scale=0.5] at (5.3, 4.8) {};  % Dot at top (same orientation)

% Rectifier diode D1 (conducts during on-time)
\draw (5, 5) to[D, l^=$D_1$] (7, 5);

% Freewheeling diode D2 (conducts during off-time, from GND to Lout input)
\draw (5, 3) -- (7, 3);
\draw (7, 3) to[D, l_=$D_2$] (7, 5);

% Output inductor
\draw (7, 5) to[L, l^="""
    latex += lout_label
    latex += r"""] (9.5, 5);

% Output capacitor
\draw (9.5, 5) -- (11, 5);
\draw (11, 5) to[C, l_="""
    latex += cout_label
    latex += r"""] (11, 0);

% Load resistor
\draw (11, 5) -- (12.5, 5);
\draw (12.5, 5) to[R, l^="""
    latex += rload_label
    latex += r"""] (12.5, 0);
\draw (12.5, 0) -- (11, 0);

% Secondary ground connection
\draw (5, 3) -- (5, 0) -- (11, 0);

% Primary ground
\draw (3, 0) node[ground] {};

% Isolation barrier (dashed line)
\draw[dashed, gray] (4.4, -0.5) -- (4.4, 7.5);
\node[gray, font=\footnotesize] at (4.4, -0.8) {Isolation};

% Output voltage label
\draw (11, 5) node[above] {$V_{out}$};

% Side labels
\node[font=\footnotesize] at (1.5, -0.8) {Primary};
\node[font=\footnotesize] at (9, -0.8) {Secondary};

\end{circuitikz}
\end{document}
"""
    return latex


def compile_latex(latex_content: str, output_path: Path, output_format: str = "pdf") -> bool:
    """Compile LaTeX content to the specified output format"""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        tex_file = tmpdir / "forward_circuit.tex"
        pdf_file = tmpdir / "forward_circuit.pdf"
        
        tex_file.write_text(latex_content)
        
        try:
            result = subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", 
                 str(tex_file)],
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode != 0:
                print(f"LaTeX compilation failed:\n{result.stdout}\n{result.stderr}")
                return False
                
        except FileNotFoundError:
            print("Error: pdflatex not found. Please install a LaTeX distribution.")
            return False
        except subprocess.TimeoutExpired:
            print("Error: LaTeX compilation timed out.")
            return False
        
        if not pdf_file.exists():
            print("Error: PDF file was not created.")
            return False
        
        if output_format == "pdf":
            shutil.copy(pdf_file, output_path)
        elif output_format == "png":
            try:
                subprocess.run(
                    ["convert", "-density", "300", str(pdf_file), 
                     "-quality", "100", str(output_path)],
                    check=True, capture_output=True
                )
            except FileNotFoundError:
                print("Error: ImageMagick 'convert' not found.")
                return False
        elif output_format == "svg":
            try:
                subprocess.run(
                    ["pdf2svg", str(pdf_file), str(output_path)],
                    check=True, capture_output=True
                )
            except FileNotFoundError:
                print("Error: pdf2svg not found.")
                return False
        
        return True


def save_latex_only(latex_content: str, output_path: Path) -> bool:
    try:
        output_path.write_text(latex_content)
        return True
    except Exception as e:
        print(f"Error saving LaTeX file: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Generate Single Switch Forward converter circuit diagram using CircuiTikz"
    )
    parser.add_argument("-o", "--output", type=str, default="forward_circuit",
                        help="Output filename (without extension)")
    parser.add_argument("-f", "--format", type=str, choices=["pdf", "png", "svg", "tex"],
                        default="pdf", help="Output format (default: pdf)")
    parser.add_argument("--no-values", action="store_true",
                        help="Don't show component values on the schematic")
    parser.add_argument("--vin", type=float, default=400.0, help="Input voltage (V)")
    parser.add_argument("--vout", type=float, default=12.0, help="Output voltage (V)")
    parser.add_argument("--iout", type=float, default=5.0, help="Output current (A)")
    
    args = parser.parse_args()
    
    params = ForwardConverterParams(
        input_voltage=args.vin,
        output_voltage=args.vout,
        output_current=args.iout,
    )
    
    show_values = not args.no_values
    latex_content = generate_forward_circuitikz(params, show_values)
    
    output_path = Path(args.output)
    if not output_path.suffix:
        output_path = output_path.with_suffix(f".{args.format}")
    
    if args.format == "tex":
        success = save_latex_only(latex_content, output_path)
    else:
        success = compile_latex(latex_content, output_path, args.format)
    
    if success:
        print(f"Successfully generated: {output_path}")
    else:
        print("Failed to generate circuit diagram.")
        tex_path = output_path.with_suffix(".tex")
        save_latex_only(latex_content, tex_path)
        print(f"LaTeX source saved to: {tex_path}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
