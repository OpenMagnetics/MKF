#!/usr/bin/env python3
"""
CircuiTikz Boost Converter Schematic Generator

This tool generates a Boost converter circuit diagram using CircuiTikz (LaTeX)
and exports it to PDF, PNG, or SVG formats.

Usage:
    python circuitikz_boost.py [--output OUTPUT] [--format FORMAT] [--show-values]

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
class BoostConverterParams:
    """Parameters for the Boost converter circuit"""
    input_voltage: float = 12.0          # V
    output_voltage: float = 24.0         # V
    output_current: float = 1.0          # A
    switching_frequency: float = 100e3   # Hz
    inductance: float = 100e-6           # H
    input_capacitance: float = 10e-6     # F
    output_capacitance: float = 100e-6   # F
    duty_cycle: float = 0.5
    rds_on: float = 0.05                 # Ohm
    diode_vf: float = 0.5                # V


def format_value(value: float, unit: str) -> str:
    """Format a value with SI prefix for display (with math mode)"""
    prefixes = [
        (1e-12, 'p'), (1e-9, 'n'), (1e-6, '\\mu '), (1e-3, 'm'),
        (1, ''), (1e3, 'k'), (1e6, 'M'), (1e9, 'G')
    ]
    
    for scale, prefix in reversed(prefixes):
        if abs(value) >= scale:
            formatted = value / scale
            if formatted == int(formatted):
                return f"${int(formatted)}\\,\\mathrm{{{prefix}{unit}}}$"
            return f"${formatted:.1f}\\,\\mathrm{{{prefix}{unit}}}$"
    
    return f"${value}\\,\\mathrm{{{unit}}}$"


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


def generate_boost_circuitikz(params: Optional[BoostConverterParams] = None, 
                               show_values: bool = True) -> str:
    """
    Generate CircuiTikz LaTeX code for a Boost converter
    
    Args:
        params: Optional circuit parameters
        show_values: If True, show component values on the schematic
        
    Returns:
        Complete LaTeX document as string
    """
    if params is None:
        params = BoostConverterParams()
    
    # Calculate load resistance
    load_r = params.output_voltage / params.output_current
    
    # Component labels (use simple labels for voltage source to avoid parsing issues)
    if show_values:
        vin_label = f"{params.input_voltage:.0f}V"
        cin_label = format_value_simple(params.input_capacitance, "F")
        l_label = format_value_simple(params.inductance, "H")
        cout_label = format_value_simple(params.output_capacitance, "F")
        rload_label = format_value_simple(load_r, r"\Omega")
    else:
        vin_label = "$V_{in}$"
        cin_label = "$C_{in}$"
        l_label = "$L$"
        cout_label = "$C_{out}$"
        rload_label = "$R_{load}$"
    
    latex = r"""\documentclass[border=10pt]{standalone}
\usepackage[siunitx,RPvoltages]{circuitikz}
\usepackage{amsmath}

\begin{document}
\begin{circuitikz}[american, scale=1.0]

% Define coordinates
\coordinate (input) at (0, 0);
\coordinate (vin_plus) at (0, 3);
\coordinate (vin_minus) at (0, 0);
\coordinate (l_in) at (2, 3);
\coordinate (sw_node) at (5, 3);
\coordinate (diode_out) at (7, 3);
\coordinate (output) at (9, 3);
\coordinate (gnd_sw) at (5, 0);
\coordinate (gnd_out) at (9, 0);

% Input voltage source
\draw (vin_minus) to[V, l="""
    latex += vin_label
    latex += r""", invert] (vin_plus);

% Input capacitor
\draw (vin_plus) -- (1, 3);
\draw (1, 3) to[C, l="""
    latex += cin_label
    latex += r"""] (1, 0);
\draw (1, 0) -- (vin_minus);

% Inductor
\draw (1, 3) -- (l_in);
\draw (l_in) to[L, l_="""
    latex += l_label
    latex += r"""] (sw_node);

% Low-side MOSFET switch
\draw (sw_node) to[Tnmos, n=mosfet, mirror] (gnd_sw);
\node[left, font=\footnotesize] at (mosfet.G) {PWM};

% Boost diode
\draw (sw_node) to[D, l=$D$] (diode_out);

% Output capacitor
\draw (diode_out) -- (output);
\draw (output) to[C, l_="""
    latex += cout_label
    latex += r"""] (gnd_out);

% Load resistor
\draw (output) -- (10.5, 3);
\draw (10.5, 3) to[R, l^="""
    latex += rload_label
    latex += r"""] (10.5, 0);
\draw (10.5, 0) -- (gnd_out);

% Ground connections
\draw (vin_minus) -- (gnd_sw) -- (gnd_out);

% Ground symbol
\draw (5, 0) node[ground] {};

% Output voltage label
\draw (output) node[above] {$V_{out}$};

% Input label  
\draw (vin_plus) node[left] {$+$};
\draw (vin_minus) node[left] {$-$};

% Switch node label
\draw (sw_node) node[above left] {$V_{sw}$};

\end{circuitikz}
\end{document}
"""
    return latex


def generate_boost_circuitikz_detailed(params: Optional[BoostConverterParams] = None,
                                        show_values: bool = True) -> str:
    """
    Generate a more detailed CircuiTikz schematic with additional annotations
    
    Args:
        params: Optional circuit parameters
        show_values: If True, show component values on the schematic
        
    Returns:
        Complete LaTeX document as string
    """
    if params is None:
        params = BoostConverterParams()
    
    load_r = params.output_voltage / params.output_current
    
    latex = r"""\documentclass[border=10pt]{standalone}
\usepackage[siunitx,RPvoltages]{circuitikz}
\usepackage{amsmath}
\usepackage{xcolor}

\begin{document}
\begin{circuitikz}[american, scale=1.2, transform shape]

% Title
\node[font=\Large\bfseries] at (5, 5) {Boost Converter};

% Define styles
\ctikzset{
    resistors/scale=0.8,
    capacitors/scale=0.8,
    inductors/scale=0.8
}

% Input Stage (left)
\draw (0, 0) to[V, v=$V_{in}$, invert] (0, 3);
\draw (0, 3) -- (1.5, 3);
\draw (1.5, 3) to[C, l_=$C_{in}$, *-*] (1.5, 0);
\draw (0, 0) -- (1.5, 0);

% Inductor
\draw (1.5, 3) to[L, l=$L$, i>^=$i_L$, -*] (4.5, 3);

% Switch node
\node[circ] at (4.5, 3) {};
\node[above] at (4.5, 3.2) {\small $V_{sw}$};

% MOSFET Switch (low-side)
\draw (4.5, 3) to[Tnmos, n=Q1, mirror] (4.5, 0);
\node[right, font=\footnotesize] at (Q1.G) {Gate};

% Gate driver block
\draw[dashed, rounded corners] (3.2, 0.3) rectangle (5.8, 1.8);
\node[font=\scriptsize] at (4.5, 1.05) {PWM};
\draw[->, thick] (4.5, 1.5) -- (Q1.G);

% Boost Diode
\draw (4.5, 3) to[D, l=$D$, -*] (7, 3);

% Output node
\node[above] at (7, 3.2) {\small $V_{out}$};

% Output Capacitor
\draw (7, 3) to[C, l=$C_{out}$, *-*] (7, 0);

% Load Resistor
\draw (7, 3) -- (9, 3);
\draw (9, 3) to[R, l=$R_{load}$] (9, 0);
\draw (9, 0) -- (7, 0);

% Ground plane
\draw[thick] (0, 0) -- (9, 0);
\draw (4.5, 0) node[ground] {};

% Current flow annotations (ON state)
\draw[->, very thick, blue!70] (0.5, 3.6) -- (4, 3.6);
\node[blue!70, font=\scriptsize, above] at (2.25, 3.6) {ON: $i_L$ charges L};

% Current flow annotations (OFF state)  
\draw[->, very thick, red!70] (5, 3.8) -- (6.5, 3.8);
\node[red!70, font=\scriptsize, above] at (5.75, 3.8) {OFF: L discharges};

% Energy flow annotations
\draw[<->, thick, green!50!black] (0, -0.8) -- (4.5, -0.8);
\node[green!50!black, font=\scriptsize] at (2.25, -1.1) {Energy storage};

\draw[<->, thick, orange] (4.5, -0.8) -- (9, -0.8);
\node[orange, font=\scriptsize] at (6.75, -1.1) {Energy transfer};

% Voltage step-up indication
\draw[<->, thick, purple] (-0.5, 0) -- (-0.5, 3);
\node[purple, font=\scriptsize, rotate=90] at (-0.8, 1.5) {$V_{in}$};

\draw[<->, thick, purple] (9.5, 0) -- (9.5, 3);
\node[purple, font=\scriptsize, rotate=90] at (9.8, 1.5) {$V_{out} > V_{in}$};

% Operating principle box
\draw[rounded corners, fill=gray!10] (0, -2.5) rectangle (9, -1.5);
\node[font=\footnotesize, align=left] at (4.5, -2) {
    Boost: $V_{out} = \frac{V_{in}}{1-D}$ \quad where $D$ = duty cycle
};

\end{circuitikz}
\end{document}
"""
    return latex


def compile_latex(latex_content: str, output_path: Path, output_format: str = "pdf") -> bool:
    """
    Compile LaTeX content to the specified output format
    
    Args:
        latex_content: LaTeX document content
        output_path: Path for the output file
        output_format: One of 'pdf', 'png', 'svg'
        
    Returns:
        True if compilation succeeded
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        tex_file = tmpdir / "boost_circuit.tex"
        pdf_file = tmpdir / "boost_circuit.pdf"
        
        # Write LaTeX file
        tex_file.write_text(latex_content)
        
        # Compile to PDF
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
        
        # Convert to final format
        if output_format == "pdf":
            shutil.copy(pdf_file, output_path)
            
        elif output_format == "png":
            try:
                subprocess.run(
                    ["convert", "-density", "300", str(pdf_file), 
                     "-quality", "100", str(output_path)],
                    check=True,
                    capture_output=True
                )
            except FileNotFoundError:
                print("Error: ImageMagick 'convert' not found. Please install ImageMagick.")
                return False
            except subprocess.CalledProcessError as e:
                print(f"Error converting to PNG: {e}")
                return False
                
        elif output_format == "svg":
            try:
                subprocess.run(
                    ["pdf2svg", str(pdf_file), str(output_path)],
                    check=True,
                    capture_output=True
                )
            except FileNotFoundError:
                print("Error: pdf2svg not found. Please install pdf2svg.")
                return False
            except subprocess.CalledProcessError as e:
                print(f"Error converting to SVG: {e}")
                return False
        else:
            print(f"Error: Unknown output format '{output_format}'")
            return False
        
        return True


def save_latex_only(latex_content: str, output_path: Path) -> bool:
    """
    Save only the LaTeX source file without compilation
    
    Args:
        latex_content: LaTeX document content
        output_path: Path for the output .tex file
        
    Returns:
        True if save succeeded
    """
    try:
        output_path.write_text(latex_content)
        return True
    except Exception as e:
        print(f"Error saving LaTeX file: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Generate Boost converter circuit diagram using CircuiTikz"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="boost_circuit",
        help="Output filename (without extension)"
    )
    parser.add_argument(
        "-f", "--format",
        type=str,
        choices=["pdf", "png", "svg", "tex"],
        default="pdf",
        help="Output format (default: pdf)"
    )
    parser.add_argument(
        "--detailed",
        action="store_true",
        help="Generate detailed schematic with annotations"
    )
    parser.add_argument(
        "--no-values",
        action="store_true",
        help="Don't show component values on the schematic"
    )
    parser.add_argument(
        "--vin",
        type=float,
        default=12.0,
        help="Input voltage (V)"
    )
    parser.add_argument(
        "--vout",
        type=float,
        default=24.0,
        help="Output voltage (V)"
    )
    parser.add_argument(
        "--iout",
        type=float,
        default=1.0,
        help="Output current (A)"
    )
    parser.add_argument(
        "--fsw",
        type=float,
        default=100e3,
        help="Switching frequency (Hz)"
    )
    parser.add_argument(
        "-L", "--inductance",
        type=float,
        default=100e-6,
        help="Inductance (H)"
    )
    
    args = parser.parse_args()
    
    # Create parameters
    params = BoostConverterParams(
        input_voltage=args.vin,
        output_voltage=args.vout,
        output_current=args.iout,
        switching_frequency=args.fsw,
        inductance=args.inductance,
        duty_cycle=1 - args.vin / args.vout  # Calculate from Vin/Vout
    )
    
    # Generate LaTeX
    show_values = not args.no_values
    if args.detailed:
        latex_content = generate_boost_circuitikz_detailed(params, show_values)
    else:
        latex_content = generate_boost_circuitikz(params, show_values)
    
    # Determine output path
    output_path = Path(args.output)
    if not output_path.suffix:
        output_path = output_path.with_suffix(f".{args.format}")
    
    # Generate output
    if args.format == "tex":
        success = save_latex_only(latex_content, output_path)
    else:
        success = compile_latex(latex_content, output_path, args.format)
    
    if success:
        print(f"Successfully generated: {output_path}")
    else:
        print("Failed to generate circuit diagram.")
        # Save LaTeX as fallback
        tex_path = output_path.with_suffix(".tex")
        save_latex_only(latex_content, tex_path)
        print(f"LaTeX source saved to: {tex_path}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
