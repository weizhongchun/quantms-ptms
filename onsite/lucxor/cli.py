#!/usr/bin/env python3
"""
Command line interface for pyLuciPHOr2
"""

import click
import os
import sys
import logging
import time
import json
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

import numpy as np
from pyopenms import (
    IdXMLFile,
    MzMLFile,
    MSExperiment,
    PeptideIdentification, 
    ProteinIdentification,
    IDFilter
)

from .psm import PSM
from .peptide import Peptide
from .models import CIDModel, HCDModel
from .constants import (
    NTERM_MOD,
    CTERM_MOD,
    AA_MASSES,
    DEFAULT_CONFIG
)
from .spectrum import Spectrum
from .flr import FLRCalculator
from .parallel import parallel_psm_processing, get_optimal_thread_count

logger = logging.getLogger(__name__)

@click.command()
@click.option('-in', '--input-spectrum', 'input_spectrum', required=True,
              help='Input spectrum file (mzML)', type=click.Path(exists=True))
@click.option('-id', '--input-id', 'input_id', required=True,
              help='Input identification file (idXML)', type=click.Path(exists=True))
@click.option('-out', '--output', 'output', required=True,
              help='Output file (idXML)', type=click.Path())
@click.option('--fragment-method', type=click.Choice(['CID', 'HCD']), default='CID',
              help='Fragmentation method (default: CID)')
@click.option('--fragment-mass-tolerance', type=float, default=0.5,
              help='Tolerance of the peaks in the fragment spectrum (default: 0.5)')
@click.option('--fragment-error-units', type=click.Choice(['Da', 'ppm']), default='Da',
              help='Unit of fragment mass tolerance (default: Da)')
@click.option('--min-mz', type=float, default=150.0,
              help='Do not consider peaks below this value (default: 150.0)')
@click.option('--target-modifications', multiple=True,
              default=["Phospho (S)", "Phospho (T)", "Phospho (Y)"],
              help='List of target modifications (default: Phospho (S) Phospho (T) Phospho (Y))')
@click.option('--neutral-losses', multiple=True,
              default=["sty -H3PO4 -97.97690"],
              help='List of neutral losses (default: sty -H3PO4 -97.97690)')
@click.option('--decoy-mass', type=float, default=79.966331,
              help='Mass to add for decoy generation (default: 79.966331)')
@click.option('--decoy-neutral-losses', multiple=True,
              default=["X -H3PO4 -97.97690"],
              help='List of decoy neutral losses (default: X -H3PO4 -97.97690)')
@click.option('--max-charge-state', type=int, default=5,
              help='Maximum charge state to consider (default: 5)')
@click.option('--max-peptide-length', type=int, default=40,
              help='Maximum peptide length (default: 40)')
@click.option('--max-num-perm', type=int, default=16384,
              help='Maximum number of permutations (default: 16384)')
@click.option('--modeling-score-threshold', type=float, default=0.95,
              help='Minimum score for modeling (default: 0.95)')
@click.option('--scoring-threshold', type=float, default=0.0,
              help='Minimum score threshold (default: 0.0)')
@click.option('--min-num-psms-model', type=int, default=50,
              help='Minimum number of PSMs for modeling (default: 50)')
@click.option('--threads', type=int, default=4,
              help='Number of threads to use (default: 4)')
@click.option('--rt-tolerance', type=float, default=0.01,
              help='Retention time tolerance (default: 0.01)')
@click.option('--debug', is_flag=True,
              help='Enable debug mode')
@click.option('--log-file', type=str,
              help='Log file path (only used in debug mode, default: {output_base}_debug.log)')
def main(input_spectrum, input_id, output, fragment_method, fragment_mass_tolerance,
         fragment_error_units, min_mz, target_modifications, neutral_losses,
         decoy_mass, decoy_neutral_losses, max_charge_state, max_peptide_length,
         max_num_perm, modeling_score_threshold, scoring_threshold,
         min_num_psms_model, threads, rt_tolerance, debug, log_file):
    """
    Modification site localization using pyLuciPHOr2 algorithm.
    
    This tool processes MS/MS spectra and peptide identifications to localize
    post-translational modifications using the LuciPHOr2 algorithm with
    false localization rate (FLR) calculation.
    """
    try:
        # Setup logging first
        setup_logging(debug, log_file, output)
        
        # Create tool instance and run
        tool = PyLuciPHOr2()
        exit_code = tool.run(input_spectrum, input_id, output, fragment_method,
                           fragment_mass_tolerance, fragment_error_units, min_mz,
                           target_modifications, neutral_losses, decoy_mass,
                           decoy_neutral_losses, max_charge_state, max_peptide_length,
                           max_num_perm, modeling_score_threshold, scoring_threshold,
                           min_num_psms_model, threads, rt_tolerance, debug)
        
        sys.exit(exit_code)
        
    except KeyboardInterrupt:
        click.echo("\nOperation cancelled by user")
        sys.exit(1)
    except Exception as e:
        click.echo(f"Error: {str(e)}")
        if debug:
            logger.error(f"Error: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
        sys.exit(1)

def setup_logging(debug, log_file, output):
    """Setup logging configuration"""
    # Configure log format
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    formatter = logging.Formatter(log_format)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG if debug else logging.INFO)
    
    # Clear existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Configure console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG if debug else logging.INFO)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)
    
    # Only configure file handler in debug mode
    if debug:
        # Get output filename (without extension)
        output_base = os.path.splitext(output)[0]
        log_file_path = log_file or f"{output_base}_debug.log"
        
        # Configure file handler
        file_handler = logging.FileHandler(log_file_path, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)
    
    # Set third-party library log levels
    logging.getLogger('numpy').setLevel(logging.WARNING)
    logging.getLogger('scipy').setLevel(logging.WARNING)

class PyLuciPHOr2:
    """Main class for pyLuciPHOr2 command line tool"""
    
    def __init__(self):
        """Initialize PyLuciPHOr2."""
        self.logger = logging.getLogger(__name__)
        self.model = None
        self.psms = []
        self.config = {}
        
    def load_input_files(self, input_id: str, input_spectrum: str) -> Tuple[List[PeptideIdentification], List[ProteinIdentification], MSExperiment]:
        """Load input files"""
        # Load identifications
        pep_ids = []
        prot_ids = []
        IdXMLFile().load(input_id, prot_ids, pep_ids)
        
        if not pep_ids:
            self.logger.warning("No peptide identifications found in input file")
            return [], [], None
            
        # Keep only best hits
        IDFilter().keepNBestHits(pep_ids, 1)
        
        # Load spectra
        exp = MSExperiment()
        MzMLFile().load(input_spectrum, exp)
        
        self.logger.info(f"Loaded {len(pep_ids)} peptide identifications")
        self.logger.info(f"Loaded {exp.size()} spectra")
        
        return pep_ids, prot_ids, exp
        
    def process_peptide(self, pep_id: PeptideIdentification, spectrum, config: Dict) -> Optional[PeptideIdentification]:
        """Process a single peptide identification"""
        try:
            # Create PSM object
            psm = PSM(pep_id, spectrum, config)
            
            # Score permutations
            psm.score_permutations(self.model)
            
            # Create new peptide identification with results
            new_pep_id = PeptideIdentification()
            new_pep_id.setMZ(pep_id.getMZ())
            new_pep_id.setRT(pep_id.getRT())
            new_pep_id.setMetaValue("spectrum_reference", pep_id.getMetaValue("spectrum_reference"))
            
            # Add scored hits
            for hit in psm.hits:
                new_pep_id.insertHit(hit)
            
            new_pep_id.assignRanks()
            return new_pep_id
            
        except Exception as e:
            self.logger.error(f"Error processing peptide: {e}")
            return None

    def run(self, input_spectrum, input_id, output, fragment_method, fragment_mass_tolerance,
            fragment_error_units, min_mz, target_modifications, neutral_losses,
            decoy_mass, decoy_neutral_losses, max_charge_state, max_peptide_length,
            max_num_perm, modeling_score_threshold, scoring_threshold,
            min_num_psms_model, threads, rt_tolerance, debug):
        """
        LuciPHOr2 main workflow:
        1. Read input files and collect all PSMs.
        2. Train CID/HCD model with high-scoring PSMs.
        3. Score all PSMs with the trained model.
        4. Exit with error if insufficient high-scoring PSMs.
        """
        config = DEFAULT_CONFIG.copy()
        
        # Parse target_modifications to handle comma-separated format
        parsed_target_modifications = []
        for mod in target_modifications:
            if ',' in mod:
                # Split comma-separated modifications
                parsed_target_modifications.extend([m.strip() for m in mod.split(',')])
            else:
                parsed_target_modifications.append(mod.strip())
        
        config.update({
            "fragment_method": fragment_method,
            "fragment_mass_tolerance": fragment_mass_tolerance,
            "fragment_error_units": fragment_error_units,
            "min_mz": min_mz,
            "target_modifications": parsed_target_modifications,
            "neutral_losses": list(neutral_losses),
            "decoy_mass": decoy_mass,
            "decoy_neutral_losses": list(decoy_neutral_losses),
            "max_charge_state": max_charge_state,
            "max_peptide_length": max_peptide_length,
            "max_num_perm": max_num_perm,
            "modeling_score_threshold": modeling_score_threshold,
            "scoring_threshold": scoring_threshold,
            "min_num_psms_model": min_num_psms_model,
            "num_threads": threads,
            "rt_tolerance": rt_tolerance
        })
        
        # Start timing
        start_time = time.time()

        self.logger.info("Loading input files...")
        self.logger.debug(f"Debug mode: {debug}")
        self.logger.debug(f"Log level: {logging.getLogger().level}")
        
        pep_ids, prot_ids, exp = self.load_input_files(input_id, input_spectrum)
        
        if not pep_ids:
            self.logger.error("No peptide identifications found")
            return 1
        
        # Collect all PSMs
        self.logger.info("Collecting PSMs...")
        all_psms = []
        
        for pep_id in pep_ids:
            # Find corresponding spectrum
            spectrum = self.find_spectrum_by_mz_rt(exp, pep_id.getMZ(), pep_id.getRT(), rt_tolerance)
            if spectrum:
                psm = PSM(pep_id, spectrum, config)
                all_psms.append(psm)
        
        self.logger.info(f"Collected {len(all_psms)} PSMs")
        
        if len(all_psms) < min_num_psms_model:
            self.logger.error(f"Insufficient PSMs for modeling: {len(all_psms)} < {min_num_psms_model}")
            return 1
        
        # Initialize model
        if fragment_method == "CID":
            self.model = CIDModel()
        else:
            self.model = HCDModel()
        
        # Train model with high-scoring PSMs
        self.logger.info("Training model...")
        high_scoring_psms = [psm for psm in all_psms if psm.score >= modeling_score_threshold]
        
        if len(high_scoring_psms) < min_num_psms_model:
            self.logger.warning(f"Only {len(high_scoring_psms)} high-scoring PSMs found, using all PSMs")
            high_scoring_psms = all_psms
        
        self.model.train(high_scoring_psms)
        
        # Score all PSMs
        self.logger.info("Scoring all PSMs...")
        if threads > 1:
            parallel_psm_processing(all_psms, self.model, threads)
        else:
            for psm in all_psms:
                psm.score_permutations(self.model)
        
        # Calculate FLR
        self.logger.info("Calculating FLR...")
        flr_calculator = FLRCalculator()
        flr_calculator.calculate_flr(all_psms)
        
        # Save results
        self.logger.info("Saving results...")
        self.save_results(output, prot_ids, all_psms)
        
        # Print summary
        elapsed = time.time() - start_time
        self.logger.info(f"Processing completed in {elapsed:.2f} seconds")
        self.logger.info(f"Processed {len(all_psms)} PSMs")
        
        return 0
    
    def find_spectrum_by_mz_rt(self, exp: MSExperiment, mz: float, rt: float, rt_tolerance: float) -> Optional[Spectrum]:
        """Find spectrum by m/z and retention time"""
        best_spectrum = None
        best_mz_diff = float('inf')
        
        for spec in exp:
            if spec.getMSLevel() == 2 and spec.getPrecursors():
                spec_mz = spec.getPrecursors()[0].getMZ()
                spec_rt = spec.getRT()
                
                mz_diff = abs(spec_mz - mz)
                rt_diff = abs(spec_rt - rt)
                
                if rt_diff <= rt_tolerance and mz_diff < best_mz_diff:
                    best_mz_diff = mz_diff
                    best_spectrum = spec
        
        return best_spectrum
    
    def save_results(self, output_file: str, prot_ids: List[ProteinIdentification], psms: List[PSM]):
        """Save results to output file"""
        pep_ids = []
        
        for psm in psms:
            new_pep_id = PeptideIdentification()
            new_pep_id.setMZ(psm.pep_id.getMZ())
            new_pep_id.setRT(psm.pep_id.getRT())
            new_pep_id.setMetaValue("spectrum_reference", psm.pep_id.getMetaValue("spectrum_reference"))
            
            # Add scored hits
            for hit in psm.hits:
                new_pep_id.insertHit(hit)
            
            new_pep_id.assignRanks()
            pep_ids.append(new_pep_id)
        
        # Save to file
        IdXMLFile().store(output_file, prot_ids, pep_ids)
        self.logger.info(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()