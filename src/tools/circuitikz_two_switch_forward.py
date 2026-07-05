#!/usr/bin/env python3
"""
CircuiTikz Two Switch Forward Converter Schematic Generator

This tool generates a Two Switch Forward converter circuit diagram using CircuiTikz (LaTeX)
and exports it to PDF, PNG, or SVG formats.

The Two Switch Forward topology uses two switches (high-side and low-side) with two
freewheeling diodes for core reset. The switch voltage stress is limited to Vin
(compared to 2*Vin for single switch forward).

Usage:
    python circuitikz_two_switch_forward.py [--output OUTPUT] [--format FORMAT] [--no-values]

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
class TwoSwitchForwardParams:
    """Parameters for the Two Switch Forward converter circuit"""
    input_voltage: float = 400.0         # V (rectified mains or high voltage DC)
    output_voltage: float = 24.0         # V
    output_current: float = 5.0          # A
    switching_frequency: float = 100e3   # Hz
    magnetizing_inductance: float = 2e-3 # H
    output_inductance: float = 50e-6     # H
    input_capacitance: float = 47e-6     # F
    output_capacitance: float = 470e-6   # F
    turns_ratio: float = 15.0            # Np/Ns
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


def generate_two_switch_forward_circuitikz(params: Optional[TwoSwitchForwardParams] = None, 
                                            show_values: bool = True) -> str:
    """
    Generate CircuiTikz LaTeX code for a Two Switch Forward converter
    
    Topology:
    - High-side MOSFET Q1 between Vin+ and transformer primary top
    - Low-side MOSFET Q2 between transformer primary bottom and GND
    - Reset diode D1 from primary top to GND (conducts during off, clamps to 0V)
    - Reset diode D2 from Vin+ to primary bottom (conducts during off, clamps to Vin)
    - Switches see max Vin (not 2*Vin like single-switch)
    - Secondary with rectifier, freewheeling diode, and output inductor
    """
    if params is None:
        params = TwoSwitchForwardParams()
    
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
\draw (0, 0) to[V, l=$V_{in}$] (0, 6);

% Input capacitor  
\draw (0, 6) -- (1.5, 6);
\draw (1.5, 6) to[C, l^="""
    latex += cin_label
    latex += r"""] (1.5, 0);
\draw (1.5, 0) -- (0, 0);

% High-side MOSFET Q1 (between Vin and primary top)
\draw (1.5, 6) -- (3, 6);
\draw (3, 6) to[Tnmos, n=q1, l=$Q_1$] (3, 4.5);
\node[left, font=\footnotesize] at (q1.G) {PWM};

% Primary winding (Np turns, dot at top)
\draw (3, 4.5) to[L, l_=$N_p$] (3, 2);
\node[circ, fill=black, scale=0.5] at (2.7, 4.3) {};  % Dot marking

% Low-side MOSFET Q2 (between primary bottom and GND)
\draw (3, 2) to[Tnmos, n=q2, l=$Q_2$, mirror] (3, 0);
\node[left, font=\footnotesize] at (q2.G) {PWM};
\draw (3, 0) -- (1.5, 0);

% Reset diode D1 (from primary top to GND - clamps pri_high to 0V during off)
\draw (3, 4.5) -- (5, 4.5);
\draw (5, 4.5) to[D, l^=$D_{r1}$] (5, 0);
\draw (5, 0) -- (3, 0);

% Reset diode D2 (from Vin to primary bottom - clamps pri_low to Vin during off)
\draw (1.5, 6) -- (1.5, 7);
\draw (1.5, 7) -- (5.5, 7);
\draw (5.5, 7) to[D, l^=$D_{r2}$] (5.5, 2);
\draw (5.5, 2) -- (3, 2);

% Transformer core (two parallel lines)
\draw[thick] (6.3, 2.5) -- (6.3, 4);
\draw[thick] (6.5, 2.5) -- (6.5, 4);

% === SECONDARY SIDE ===
% Secondary winding (Ns turns)
\draw (7, 4.5) to[L, l^=$N_s$] (7, 2);
\node[circ, fill=black, scale=0.5] at (7.3, 4.3) {};  % Dot at top (same orientation as primary)

% Rectifier diode D3 (conducts during on-time)
\draw (7, 4.5) to[D, l^=$D_1$] (9, 4.5);

% Freewheeling diode D4 (conducts during off-time)
\draw (7, 2) -- (9, 2);
\draw (9, 2) to[D, l_=$D_2$] (9, 4.5);

% Output inductor
\draw (9, 4.5) to[L, l^="""
    latex += lout_label
    latex += r"""] (11.5, 4.5);

% Output capacitor
\draw (11.5, 4.5) -- (13, 4.5);
\draw (13, 4.5) to[C, l_="""
    latex += cout_label
    latex += r"""] (13, 0);

% Load resistor
\draw (13, 4.5) -- (14.5, 4.5);
\draw (14.5, 4.5) to[R, l^="""
    latex += rload_label
    latex += r"""] (14.5, 0);
\draw (14.5, 0) -- (13, 0);

% Secondary ground connection
\draw (7, 2) -- (7, 0) -- (13, 0);

% Primary ground
\draw (3, 0) node[ground] {};

% Isolation barrier (dashed line)
\draw[dashed, gray] (6.4, -0.5) -- (6.4, 7.5);
\node[gray, font=\footnotesize] at (6.4, -0.8) {Isolation};

% Output voltage label
\draw (13, 4.5) node[above] {$V_{out}$};

% Side labels
\node[font=\footnotesize] at (2.5, -0.8) {Primary};
\node[font=\footnotesize] at (11, -0.8) {Secondary};

% Node labels for key points
\node[font=\scriptsize, blue] at (3, 4.8) {pri\_high};
\node[font=\scriptsize, blue] at (3.3, 1.7) {pri\_low};

\end{circuitikz}
\end{document}
"""
    return latex


def compile_latex(latex_content: str, output_path: Path, output_format: str = "pdf") -> bool:
    """Compile LaTeX content to the specified output format"""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        tex_file = tmpdir / "two_switch_forward_circuit.tex"
        pdf_file = tmpdir / "two_switch_forward_circuit.pdf"
        
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
            
            if not pdf_file.exists():
                print("PDF file was not created")
                return False
            
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            if output_format.lower() == "pdf":
                shutil.copy(pdf_file, output_path)
            elif output_format.lower() == "svg":
                svg_file = tmpdir / "two_switch_forward_circuit.svg"
                result = subprocess.run(
                    ["pdf2svg", str(pdf_file), str(svg_file)],
                    capture_output=True, text=True, timeout=30
                )
                if result.returncode != 0:
                    print(f"pdf2svg conversion failed:\n{result.stderr}")
                    return False
                shutil.copy(svg_file, output_path)
            elif output_format.lower() == "png":
                result = subprocess.run(
                    ["convert", "-density", "300", str(pdf_file), 
                     "-quality", "100", str(output_path)],
                    capture_output=True, text=True, timeout=30
                )
                if result.returncode != 0:
                    print(f"ImageMagick conversion failed:\n{result.stderr}")
                    return False
            else:
                print(f"Unknown output format: {output_format}")
                return False
            
            print(f"Successfully created: {output_path}")
            return True
            
        except subprocess.TimeoutExpired:
            print("Compilation timed out")
            return False
        except FileNotFoundError as e:
            print(f"Required tool not found: {e}")
            return False


def main():
    parser = argparse.ArgumentParser(
        description="Generate Two Switch Forward converter schematic using CircuiTikz"
    )
    parser.add_argument(
        "--output", "-o",
        default="two_switch_forward_schematic.pdf",
        help="Output file path (default: two_switch_forward_schematic.pdf)"
    )
    parser.add_argument(
        "--format", "-f",
        choices=["pdf", "svg", "png"],
        default=None,
        help="Output format (auto-detected from extension if not specified)"
    )
    parser.add_argument(
        "--no-values",
        action="store_true",
        help="Hide component values (show only symbolic names)"
    )
    parser.add_argument(
        "--latex-only",
        action="store_true",
        help="Only output the LaTeX code to stdout (don't compile)"
    )
    parser.add_argument(
        "--vin", type=float, default=400.0,
        help="Input voltage in V (default: 400)"
    )
    parser.add_argument(
        "--vout", type=float, default=24.0,
        help="Output voltage in V (default: 24)"
    )
    parser.add_argument(
        "--iout", type=float, default=5.0,
        help="Output current in A (default: 5)"
    )
    parser.add_argument(
        "--fsw", type=float, default=100e3,
        help="Switching frequency in Hz (default: 100000)"
    )
    
    args = parser.parse_args()
    
    params = TwoSwitchForwardParams(
        input_voltage=args.vin,
        output_voltage=args.vout,
        output_current=args.iout,
        switching_frequency=args.fsw
    )
    
    latex = generate_two_switch_forward_circuitikz(params, not args.no_values)
    
    if args.latex_only:
        print(latex)
        return
    
    output_path = Path(args.output)
    if args.format:
        output_format = args.format
    else:
        output_format = output_path.suffix.lstrip(".").lower() or "pdf"
    
    success = compile_latex(latex, output_path, output_format)
    exit(0 if success else 1)


if __name__ == "__main__":
    main()
