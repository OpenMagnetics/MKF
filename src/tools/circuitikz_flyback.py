#!/usr/bin/env python3
"""
CircuiTikz Flyback Converter Schematic Generator

This tool generates a Flyback converter circuit diagram using CircuiTikz (LaTeX)
and exports it to PDF, PNG, or SVG formats.

Usage:
    python circuitikz_flyback.py [--output OUTPUT] [--format FORMAT] [--show-values]

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
class FlybackConverterParams:
    """Parameters for the Flyback converter circuit"""
    input_voltage: float = 400.0         # V (rectified mains)
    output_voltage: float = 12.0         # V
    output_current: float = 1.0          # A
    switching_frequency: float = 100e3   # Hz
    magnetizing_inductance: float = 500e-6  # H
    input_capacitance: float = 47e-6     # F
    output_capacitance: float = 470e-6   # F
    turns_ratio: float = 10.0            # Np/Ns
    duty_cycle: float = 0.4


def format_value_simple(value: float, unit: str) -> str:
    """Format a value with SI prefix for display (simple text, no math mode)"""
    prefixes = [
        (1e-12, 'p'), (1e-9, 'n'), (1e-6, 'u'), (1e-3, 'm'),
        (1, ''), (1e3, 'k'), (1e6, 'M'), (1e9, 'G')
    ]
    
    # Handle ohms specially - needs math mode for the symbol
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


def generate_flyback_circuitikz(params: Optional[FlybackConverterParams] = None, 
                                 show_values: bool = True) -> str:
    """
    Generate CircuiTikz LaTeX code for a Flyback converter
    
    Args:
        params: Optional circuit parameters
        show_values: If True, show component values on the schematic
        
    Returns:
        Complete LaTeX document as string
    """
    if params is None:
        params = FlybackConverterParams()
    
    # Calculate load resistance
    load_r = params.output_voltage / params.output_current
    
    # Component labels
    if show_values:
        vin_label = f"{params.input_voltage:.0f}V"
        cin_label = format_value_simple(params.input_capacitance, "F")
        cout_label = format_value_simple(params.output_capacitance, "F")
        rload_label = format_value_simple(load_r, r"\Omega")
    else:
        vin_label = "$V_{in}$"
        cin_label = "$C_{in}$"
        cout_label = "$C_{out}$"
        rload_label = "$R_{load}$"
    
    latex = r"""\documentclass[border=10pt]{standalone}
\usepackage[siunitx,RPvoltages]{circuitikz}
\usepackage{amsmath}

\begin{document}
\begin{circuitikz}[american, scale=1.0]

% Define coordinates - Primary side
\coordinate (vin_plus) at (0, 4);
\coordinate (vin_minus) at (0, 0);
\coordinate (pri_top) at (3.5, 4);
\coordinate (pri_bot) at (3.5, 2);
\coordinate (sw_node) at (3.5, 2);
\coordinate (gnd_pri) at (3.5, 0);

% Define coordinates - Secondary side
\coordinate (sec_top) at (6, 4);
\coordinate (sec_bot) at (6, 2);
\coordinate (diode_out) at (8, 4);
\coordinate (output) at (10, 4);
\coordinate (gnd_sec) at (10, 0);

% Input voltage source
\draw (vin_plus) to[V, l="""
    latex += vin_label
    latex += r"""] (vin_minus);

% Input capacitor
\draw (vin_plus) -- (1.8, 4);
\draw (1.8, 4) to[C, l^="""
    latex += cin_label
    latex += r"""] (1.8, 0);
\draw (1.8, 0) -- (vin_minus);

% Primary winding (with dot)
\draw (1.8, 4) -- (pri_top);
\draw (pri_top) to[L, l_=$L_p$] (pri_bot);
\node[circ, fill=black, scale=0.5] at (3.2, 3.8) {};  % Dot marking

% Low-side MOSFET switch
\draw (pri_bot) to[Tnmos, n=mosfet, mirror] (gnd_pri);
\node[left, font=\footnotesize] at (mosfet.G) {PWM};

% RCD Clamp circuit - positioned above and to the side
\draw (sw_node) -- (4.5, 2);
\draw (4.5, 2) to[D, l_=$D_c$] (4.5, 4);
\draw (4.5, 4) -- (4.5, 5.5);
\draw (4.5, 5.5) to[C, l^=$C_c$] (2.5, 5.5);
\draw (2.5, 5.5) to[R, l^=$R_c$] (2.5, 4);
\draw (2.5, 4) -- (1.8, 4);

% Secondary winding (with dot - opposite polarity for flyback)
\draw (sec_top) to[L, l^=$L_s$] (sec_bot);
\node[circ, fill=black, scale=0.5] at (6.3, 2.2) {};  % Dot marking (bottom for flyback)

% Transformer core symbol (two parallel lines)
\draw[thick] (5.4, 2.3) -- (5.4, 3.7);
\draw[thick] (5.6, 2.3) -- (5.6, 3.7);

% Secondary rectifier diode
\draw (sec_top) to[D, l^=$D$] (diode_out);

% Output capacitor
\draw (diode_out) -- (output);
\draw (output) to[C, l_="""
    latex += cout_label
    latex += r"""] (gnd_sec);

% Load resistor
\draw (output) -- (11.5, 4);
\draw (11.5, 4) to[R, l^="""
    latex += rload_label
    latex += r"""] (11.5, 0);
\draw (11.5, 0) -- (gnd_sec);

% Ground connections
\draw (vin_minus) -- (gnd_pri);
\draw (sec_bot) -- (6, 0) -- (gnd_sec);

% Ground symbols
\draw (3.5, 0) node[ground] {};

% Isolation barrier (dashed line)
\draw[dashed, gray] (5.5, -0.5) -- (5.5, 6);
\node[gray, font=\footnotesize] at (5.5, -0.8) {Isolation};

% Labels
\draw (output) node[above] {$V_{out}$};

% Primary/Secondary labels
\node[font=\footnotesize] at (2, -0.8) {Primary};
\node[font=\footnotesize] at (9, -0.8) {Secondary};

\end{circuitikz}
\end{document}
"""
    return latex


def compile_latex(latex_content: str, output_path: Path, output_format: str = "pdf") -> bool:
    """
    Compile LaTeX content to the specified output format
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        tex_file = tmpdir / "flyback_circuit.tex"
        pdf_file = tmpdir / "flyback_circuit.pdf"
        
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
        description="Generate Flyback converter circuit diagram using CircuiTikz"
    )
    parser.add_argument("-o", "--output", type=str, default="flyback_circuit",
                        help="Output filename (without extension)")
    parser.add_argument("-f", "--format", type=str, choices=["pdf", "png", "svg", "tex"],
                        default="pdf", help="Output format (default: pdf)")
    parser.add_argument("--no-values", action="store_true",
                        help="Don't show component values on the schematic")
    parser.add_argument("--vin", type=float, default=400.0, help="Input voltage (V)")
    parser.add_argument("--vout", type=float, default=12.0, help="Output voltage (V)")
    parser.add_argument("--iout", type=float, default=1.0, help="Output current (A)")
    
    args = parser.parse_args()
    
    params = FlybackConverterParams(
        input_voltage=args.vin,
        output_voltage=args.vout,
        output_current=args.iout,
    )
    
    show_values = not args.no_values
    latex_content = generate_flyback_circuitikz(params, show_values)
    
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
