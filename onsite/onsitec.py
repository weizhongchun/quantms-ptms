#!/usr/bin/env python3
"""
OnSite CLI - Integrated command line interface for phosphorylation site localization tools.

This module provides a unified interface for accessing AScore, PhosphoRS, and LucXor algorithms.
"""

import argparse
import sys
import os
from typing import List, Optional

def create_parser():
    """Create the main argument parser for OnSite CLI."""
    parser = argparse.ArgumentParser(
        prog='onsite',
        description='OnSite: Mass spectrometry post-translational modification localization tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Available algorithms:
  ascore      AScore algorithm for phosphorylation site localization
  phosphors   PhosphoRS algorithm for phosphorylation site localization  
  lucxor      LucXor (LuciPHOr2) algorithm for PTM localization

Examples:
  onsite ascore -in spectra.mzML -id identifications.idXML -out results.idXML
  onsite phosphors -in spectra.mzML -id identifications.idXML -out results.idXML
  onsite lucxor -in spectra.mzML -id identifications.idXML -out results.idXML

For algorithm-specific help:
  onsite ascore --help
  onsite phosphors --help
  onsite lucxor --help
        """
    )
    
    subparsers = parser.add_subparsers(
        dest='algorithm',
        help='Algorithm to use for phosphorylation site localization',
        metavar='ALGORITHM'
    )
    
    # AScore subcommand
    ascore_parser = subparsers.add_parser(
        'ascore',
        help='AScore algorithm for phosphorylation site localization',
        description='AScore algorithm for phosphorylation site localization'
    )
    _setup_ascore_parser(ascore_parser)
    
    # PhosphoRS subcommand
    phosphors_parser = subparsers.add_parser(
        'phosphors',
        help='PhosphoRS algorithm for phosphorylation site localization',
        description='PhosphoRS algorithm for phosphorylation site localization'
    )
    _setup_phosphors_parser(phosphors_parser)
    
    # LucXor subcommand
    lucxor_parser = subparsers.add_parser(
        'lucxor',
        help='LucXor (LuciPHOr2) algorithm for PTM localization',
        description='LucXor (LuciPHOr2) algorithm for PTM localization'
    )
    _setup_lucxor_parser(lucxor_parser)
    
    return parser

def _setup_ascore_parser(parser):
    """Setup AScore-specific arguments."""
    # Required arguments
    parser.add_argument('-in', dest='in_file', required=True,
                       help='Input mzML file path')
    parser.add_argument('-id', dest='id_file', required=True,
                       help='Input idXML file path')
    parser.add_argument('-out', dest='out_file', required=True,
                       help='Output idXML file path')
    
    # Optional arguments
    parser.add_argument('-fragment_mass_tolerance', type=float, default=0.05,
                       help='Fragment mass tolerance value (default: 0.05)')
    parser.add_argument('-fragment_mass_unit', choices=['Da', 'ppm'], default='Da',
                       help='Tolerance unit (default: Da)')
    parser.add_argument('--threads', dest='threads', type=int, default=4,
                       help='Number of parallel threads (default: 4)')
    parser.add_argument('-debug', dest='debug', action='store_true',
                       help='Enable debug output and write debug log')
    parser.add_argument('--add_decoys', dest='add_decoys', action='store_true', default=False,
                       help='Include A (PhosphoDecoy) as potential phosphorylation site')

def _setup_phosphors_parser(parser):
    """Setup PhosphoRS-specific arguments."""
    # Required arguments
    parser.add_argument('-in', dest='in_file', required=True,
                       help='Input mzML file path')
    parser.add_argument('-id', dest='id_file', required=True,
                       help='Input idXML file path')
    parser.add_argument('-out', dest='out_file', required=True,
                       help='Output idXML file path')
    
    # Optional arguments
    parser.add_argument('-fragment_mass_tolerance', type=float, default=0.05,
                       help='Fragment mass tolerance value (default: 0.05)')
    parser.add_argument('-fragment_mass_unit', choices=['Da', 'ppm'], default='Da',
                       help='Tolerance unit (default: Da)')
    parser.add_argument('--threads', dest='threads', type=int, default=1,
                       help='Number of parallel processes (default: 1)')
    parser.add_argument('-debug', dest='debug', action='store_true',
                       help='Enable debug output and write debug log')
    parser.add_argument('--add_decoys', dest='add_decoys', action='store_true', default=False,
                       help='Include A (PhosphoDecoy) as potential phosphorylation site')

def _setup_lucxor_parser(parser):
    """Setup LucXor-specific arguments."""
    # Required arguments
    parser.add_argument('-in', '--input_spectrum', required=True,
                       help='Input spectrum file (mzML)')
    parser.add_argument('-id', '--input_id', required=True,
                       help='Input identification file (idXML)')
    parser.add_argument('-out', '--output', required=True,
                       help='Output file (idXML)')
    
    # Optional arguments
    parser.add_argument('--fragment_method', choices=['CID', 'HCD'], default='CID',
                       help='Fragmentation method (default: CID)')
    parser.add_argument('--fragment_mass_tolerance', type=float, default=0.5,
                       help='Tolerance of the peaks in the fragment spectrum (default: 0.5)')
    parser.add_argument('--fragment_error_units', choices=['Da', 'ppm'], default='Da',
                       help='Unit of fragment mass tolerance (default: Da)')
    parser.add_argument('--min_mz', type=float, default=150.0,
                       help='Do not consider peaks below this value (default: 150.0)')
    parser.add_argument('--target_modifications', nargs='+',
                       default=["Phospho (S)", "Phospho (T)", "Phospho (Y)"],
                       help='List of target modifications (default: Phospho (S) Phospho (T) Phospho (Y))')
    parser.add_argument('--neutral_losses', nargs='+',
                       default=["sty -H3PO4 -97.97690"],
                       help='List of neutral losses (default: sty -H3PO4 -97.97690)')
    parser.add_argument('--decoy_mass', type=float, default=79.966331,
                       help='Mass to add for decoy generation (default: 79.966331)')
    parser.add_argument('--decoy_neutral_losses', nargs='+',
                       default=["X -H3PO4 -97.97690"],
                       help='List of decoy neutral losses (default: X -H3PO4 -97.97690)')
    parser.add_argument('--max_charge_state', type=int, default=5,
                       help='Maximum charge state to consider (default: 5)')
    parser.add_argument('--max_peptide_length', type=int, default=40,
                       help='Maximum peptide length (default: 40)')
    parser.add_argument('--max_num_perm', type=int, default=16384,
                       help='Maximum number of permutations (default: 16384)')
    parser.add_argument('--modeling_score_threshold', type=float, default=0.95,
                       help='Minimum score for modeling (default: 0.95)')
    parser.add_argument('--scoring_threshold', type=float, default=0.0,
                       help='Minimum score threshold (default: 0.0)')
    parser.add_argument('--min_num_psms_model', type=int, default=50,
                       help='Minimum number of PSMs for modeling (default: 50)')
    parser.add_argument('--threads', type=int, default=4,
                       help='Number of threads to use (default: 4)')
    parser.add_argument('--rt_tolerance', type=float, default=0.01,
                       help='Retention time tolerance (default: 0.01)')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug mode')
    parser.add_argument('--log_file', type=str,
                       help='Log file path (only used in debug mode)')

def run_ascore(args):
    """Run AScore algorithm."""
    # Import and run phospho_scoring.py
    from onsite.phospho_scoring import main as ascore_main
    
    # Set up sys.argv for phospho_scoring
    original_argv = sys.argv[:]
    sys.argv = ['phospho_scoring.py']
    
    # Add arguments
    sys.argv.extend(['-in', args.in_file])
    sys.argv.extend(['-id', args.id_file])
    sys.argv.extend(['-out', args.out_file])
    sys.argv.extend(['-fragment_mass_tolerance', str(args.fragment_mass_tolerance)])
    sys.argv.extend(['-fragment_mass_unit', args.fragment_mass_unit])
    sys.argv.extend(['-threads', str(args.threads)])
    if args.debug:
        sys.argv.append('-debug')
    if args.add_decoys:
        sys.argv.append('--add_decoys')
    
    try:
        ascore_main()
    finally:
        sys.argv = original_argv

def run_phosphors(args):
    """Run PhosphoRS algorithm."""
    # Import and run phosphors_scoring.py
    from onsite.phosphors_scoring import main as phosphors_main
    
    # Set up sys.argv for phosphors_scoring
    original_argv = sys.argv[:]
    sys.argv = ['phosphors_scoring.py']
    
    # Add arguments
    sys.argv.extend(['-in', args.in_file])
    sys.argv.extend(['-id', args.id_file])
    sys.argv.extend(['-out', args.out_file])
    sys.argv.extend(['-fragment_mass_tolerance', str(args.fragment_mass_tolerance)])
    sys.argv.extend(['-fragment_mass_unit', args.fragment_mass_unit])
    sys.argv.extend(['-threads', str(args.threads)])
    if args.debug:
        sys.argv.append('-debug')
    if args.add_decoys:
        sys.argv.append('--add_decoys')
    
    try:
        phosphors_main()
    finally:
        sys.argv = original_argv

def run_lucxor(args):
    """Run LucXor algorithm."""
    # Import and run lucxor CLI
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from lucxor.cli import main as lucxor_main
    
    # Set up sys.argv for LucXor
    original_argv = sys.argv[:]
    sys.argv = ['lucxor/cli.py']
    
    # Add arguments
    sys.argv.extend(['-in', args.input_spectrum])
    sys.argv.extend(['-id', args.input_id])
    sys.argv.extend(['-out', args.output])
    sys.argv.extend(['--fragment_method', args.fragment_method])
    sys.argv.extend(['--fragment_mass_tolerance', str(args.fragment_mass_tolerance)])
    sys.argv.extend(['--fragment_error_units', args.fragment_error_units])
    sys.argv.extend(['--min_mz', str(args.min_mz)])
    sys.argv.extend(['--target_modifications'] + args.target_modifications)
    sys.argv.extend(['--neutral_losses'] + args.neutral_losses)
    sys.argv.extend(['--decoy_mass', str(args.decoy_mass)])
    sys.argv.extend(['--decoy_neutral_losses'] + args.decoy_neutral_losses)
    sys.argv.extend(['--max_charge_state', str(args.max_charge_state)])
    sys.argv.extend(['--max_peptide_length', str(args.max_peptide_length)])
    sys.argv.extend(['--max_num_perm', str(args.max_num_perm)])
    sys.argv.extend(['--modeling_score_threshold', str(args.modeling_score_threshold)])
    sys.argv.extend(['--scoring_threshold', str(args.scoring_threshold)])
    sys.argv.extend(['--min_num_psms_model', str(args.min_num_psms_model)])
    sys.argv.extend(['--threads', str(args.threads)])
    sys.argv.extend(['--rt_tolerance', str(args.rt_tolerance)])
    if args.debug:
        sys.argv.append('--debug')
    if args.log_file:
        sys.argv.extend(['--log_file', args.log_file])
    
    try:
        lucxor_main()
    finally:
        sys.argv = original_argv

def main():
    """Main entry point for OnSite CLI."""
    parser = create_parser()
    args = parser.parse_args()
    
    if not args.algorithm:
        parser.print_help()
        sys.exit(1)
    
    try:
        if args.algorithm == 'ascore':
            run_ascore(args)
        elif args.algorithm == 'phosphors':
            run_phosphors(args)
        elif args.algorithm == 'lucxor':
            run_lucxor(args)
        else:
            print(f"Unknown algorithm: {args.algorithm}")
            sys.exit(1)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user")
        sys.exit(1)
    except Exception as e:
        print(f"Error running {args.algorithm}: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    main()
