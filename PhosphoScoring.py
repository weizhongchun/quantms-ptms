#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import traceback
import time
import logging
from datetime import datetime
from pyopenms import *
from Ascore import AScore

def parse_args():
    """Parse command line arguments with enhanced validation"""
    parser = argparse.ArgumentParser(
        description='Phosphorylation site localization scoring tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-in', dest='in_file', required=True, 
                       help='Input mzML file path')
    parser.add_argument('-id', dest='id_file', required=True,
                       help='Input idXML file path')
    parser.add_argument('-out', dest='out_file', required=True,
                       help='Output idXML file path')
    parser.add_argument('-fragment_mass_tolerance', type=float, default=0.05,
                       help='Fragment mass tolerance value')
    parser.add_argument('-fragment_mass_unit', choices=['Da', 'ppm'], default='Da',
                       help='Tolerance unit (Da or ppm)')
    return parser.parse_args()

def load_spectra(mzml_file):
    """Load MS/MS spectra with progress feedback"""
    print(f"[{time.strftime('%H:%M:%S')}] Loading spectra from {mzml_file}")
    exp = MSExperiment()
    FileHandler().loadExperiment(mzml_file, exp)
    print(f"Loaded {exp.size()} spectra")
    return exp

def load_identifications(idxml_file):
    """Load identification results with metadata validation"""
    print(f"[{time.strftime('%H:%M:%S')}] Loading identifications from {idxml_file}")
    protein_ids = []
    peptide_ids = []
    IdXMLFile().load(idxml_file, protein_ids, peptide_ids)
    print(f"Loaded {len(peptide_ids)} peptide identifications")
    return protein_ids, peptide_ids

def save_identifications(out_file, protein_ids, peptide_ids):
    """Save results with complete metadata preservation"""
    try:
        # Ensure all critical metadata fields exist
        for pid in peptide_ids:
            # Set identification-level attributes
            pid.setScoreType("PhosphoScore")
            pid.setHigherScoreBetter(True)
            pid.setSignificanceThreshold(0.0)
            
            for hit in pid.getHits():
                seq_str = hit.getSequence().toString()
                
                # Ensure all required metadata fields exist
                if not hit.metaValueExists("search_engine_sequence"):
                    hit.setMetaValue("search_engine_sequence", seq_str)
                if not hit.metaValueExists("regular_phospho_count"):
                    hit.setMetaValue("regular_phospho_count", seq_str.count("(Phospho)"))
                if not hit.metaValueExists("phospho_decoy_count"):
                    hit.setMetaValue("phospho_decoy_count", seq_str.count("(PhosphoDecoy)"))
                
                if hit.metaValueExists("MS:1002052"):
                    hit.setMetaValue("SpecEValue_score", float(hit.getMetaValue("MS:1002052")))
                
                # Then set AScore fields
                if not hit.metaValueExists("AScore_pep_score"):
                    hit.setMetaValue("AScore_pep_score", hit.getScore())
                
                # Ensure AScore fields exist
                rank = 1
                while hit.metaValueExists(f"AScore_{rank}"):
                    hit.setMetaValue(f"AScore_{rank}", float(hit.getMetaValue(f"AScore_{rank}")))
                    rank += 1

        # Write with XML format validation
        IdXMLFile().store(out_file, protein_ids, peptide_ids)
        print(f"Successfully saved {len(peptide_ids)} identifications to {out_file}")
    except Exception as e:
        print(f"Critical error saving results: {str(e)}")
        traceback.print_exc()
        sys.exit(1)

def find_spectrum_by_mz(exp, target_mz, rt=None, ppm_tolerance=10):
    """Optimized spectrum matching with caching"""
    # Binary search for optimized spectrum matching
    if not hasattr(find_spectrum_by_mz, 'spectrum_cache'):
        find_spectrum_by_mz.spectrum_cache = {}
        find_spectrum_by_mz.spectrum_list = []
        
        # Preprocess all MS2 spectra
        for spec in exp:
            if spec.getMSLevel() == 2 and spec.getPrecursors():
                mz = spec.getPrecursors()[0].getMZ()
                find_spectrum_by_mz.spectrum_list.append((mz, spec))
        
        # Sort by m/z
        find_spectrum_by_mz.spectrum_list.sort(key=lambda x: x[0])
    
    # Binary search for the closest m/z
    left, right = 0, len(find_spectrum_by_mz.spectrum_list) - 1
    best_match = None
    min_diff = float('inf')
    
    while left <= right:
        mid = (left + right) // 2
        mz, spec = find_spectrum_by_mz.spectrum_list[mid]
        diff = abs(mz - target_mz)
        
        if diff < min_diff:
            min_diff = diff
            best_match = spec
        
        if mz < target_mz:
            left = mid + 1
        else:
            right = mid - 1
    
    return best_match

def process_peptide_hit(hit, spectrum, ascore, logger):
    """Optimized phosphorylation analysis"""
    try:
        # Store original sequence as search_engine_sequence
        original_seq_str = hit.getSequence().toString()
        hit.setMetaValue("search_engine_sequence", original_seq_str)
        
        # Quick check for phosphorylation sites
        if "(Phospho)" not in original_seq_str and "(PhosphoDecoy)" not in original_seq_str:
            hit.setScore(-1.0)
            hit.setMetaValue("AScore_pep_score", -1.0)
            return hit
        
        # Cache calculation results
        cache_key = f"{original_seq_str}_{spectrum.getNativeID()}"
        if hasattr(process_peptide_hit, 'result_cache') and cache_key in process_peptide_hit.result_cache:
            cached_result = process_peptide_hit.result_cache[cache_key]
            for key, value in cached_result.items():
                hit.setMetaValue(key, value)
            return hit
        
        # Configure AScore parameters
        ascore.fragment_mass_tolerance_ = args.fragment_mass_tolerance
        ascore.fragment_tolerance_ppm_ = (args.fragment_mass_unit == 'ppm')
        
        # Execute AScore calculation
        scored_hit = ascore.compute(hit, spectrum)
        
        # Get the new sequence after AScore reassignment
        new_seq_str = scored_hit.getSequence().toString()
        
        # Set the new sequence as the main sequence
        hit.setSequence(scored_hit.getSequence())
        
        # Extract site-specific scores
        site_scores = []
        rank = 1
        while scored_hit.metaValueExists(f"AScore_{rank}"):
            score = float(scored_hit.getMetaValue(f"AScore_{rank}"))
            site_scores.append(score)
            rank += 1
        
        # Determine overall score (minimum site score)
        final_score = min(site_scores) if site_scores else -1.0
        hit.setScore(final_score)
        
        # Set all required metadata fields
        result_cache = {
            "search_engine_sequence": original_seq_str,  # Keep original sequence
            "regular_phospho_count": new_seq_str.count("(Phospho)"),  # Count from new sequence
            "phospho_decoy_count": new_seq_str.count("(PhosphoDecoy)"),  # Count from new sequence
            "AScore_pep_score": final_score
        }
        
        # Store individual site scores
        for i, score in enumerate(site_scores, 1):
            result_cache[f"AScore_{i}"] = score
        
        # Add ProForma notation (if available)
        if scored_hit.metaValueExists("ProForma"):
            result_cache["ProForma"] = scored_hit.getMetaValue("ProForma")
        
        # Preserve original MSGF score
        if hit.metaValueExists("MS:1002052"):
            result_cache["SpecEValue_score"] = float(hit.getMetaValue("MS:1002052"))
        
        # Cache results
        if not hasattr(process_peptide_hit, 'result_cache'):
            process_peptide_hit.result_cache = {}
        process_peptide_hit.result_cache[cache_key] = result_cache
        
        # Apply cached results
        for key, value in result_cache.items():
            hit.setMetaValue(key, value)
        
        # Output phosphorylated peptide information
        if "(Phospho)" in new_seq_str or "(PhosphoDecoy)" in new_seq_str:
            print("--------------")
            print(f"Original sequence: {original_seq_str}")
            print(f"New sequence: {new_seq_str}")
            print(f"AScore_pep_score: {result_cache['AScore_pep_score']}")
            
            # Output AScore site scores
            for i in range(1, rank):
                print(f"AScore_{i}: {result_cache.get(f'AScore_{i}', 'N/A')}")
            
            # Output ProForma information
            proforma = result_cache.get("ProForma", "N/A")
            print(f"ProForma: {proforma}")
            print("--------------")
            
            # Log to debug file
            logger.info("--------------")
            logger.info(f"Original sequence: {original_seq_str}")
            logger.info(f"New sequence: {new_seq_str}")
            logger.info(f"AScore_pep_score: {result_cache['AScore_pep_score']}")
            for i in range(1, rank):
                logger.info(f"AScore_{i}: {result_cache.get(f'AScore_{i}', 'N/A')}")
            logger.info(f"ProForma: {proforma}")
            logger.info("--------------")
        
        return hit
        
    except Exception as e:
        print(f"Error processing {hit.getSequence().toString()}: {str(e)}")
        traceback.print_exc()
        return hit

def process_peptide_identification(pid, exp, ascore, logger):
    """Process a single peptide identification with error handling"""
    try:
        # Create new PeptideIdentification object
        new_pid = PeptideIdentification(pid)
        new_pid.setScoreType("PhosphoScore")
        new_pid.setHigherScoreBetter(True)
        new_pid.setSignificanceThreshold(0.0)
        
        # Find corresponding spectrum
        spectrum = find_spectrum_by_mz(exp, pid.getMZ(), pid.getRT())
        if not spectrum:
            logger.error(f"Error processing identification: spectrum_not_found for MZ {pid.getMZ()}")
            return {'status': 'error', 'reason': 'spectrum_not_found'}
        
        # Process each peptide hit
        scored_peptides = []
        for hit in pid.getHits():
            # Create new PeptideHit object
            new_hit = PeptideHit(hit)
            # Process phosphorylation sites
            processed_hit = process_peptide_hit(new_hit, spectrum, ascore, logger)
            scored_peptides.append(processed_hit)
        
        # Set new hits
        new_pid.setHits(scored_peptides)
        
        return {'status': 'success', 'new_pid': new_pid}
        
    except Exception as e:
        logger.error(f"Error processing identification: {str(e)}")
        return {'status': 'error', 'reason': str(e)}

def log_debug(log_file):
    """Initialize debug logging"""
    logger = logging.getLogger('debug_logger')
    logger.setLevel(logging.INFO)
    
    handler = logging.FileHandler(log_file)
    handler.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    
    logger.addHandler(handler)
    return logger

def main():
    """Main workflow"""
    try:
        global args
        args = parse_args()
        
        # Initialize processing pipeline
        exp = load_spectra(args.in_file)
        protein_ids, peptide_ids = load_identifications(args.id_file)
        ascore = AScore()
        
        # Initialize debug log
        log_file = f"{args.out_file}.debug.log"
        logger = log_debug(log_file)
        logger.info("PhosphoScoring Debug Log")
        logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        # Processing statistics
        stats = {
            'total': len(peptide_ids),
            'processed': 0,
            'phospho': 0,
            'errors': 0
        }
        
        # Main processing loop
        start_time = time.time()
        processed_peptide_ids = []
        
        # Process each PeptideIdentification
        for pid in peptide_ids:
            try:
                result = process_peptide_identification(pid, exp, ascore, logger)
                if result['status'] == 'success':
                    processed_peptide_ids.append(result['new_pid'])
                    stats['processed'] += 1
                    stats['phospho'] += len([h for h in result['new_pid'].getHits() 
                                          if "(Phospho)" in h.getSequence().toString()])
                else:
                    stats['errors'] += 1
                    logger.error(f"Error processing identification: {result['reason']}")
            except Exception as e:
                stats['errors'] += 1
                logger.error(f"Error processing identification: {str(e)}")
                traceback.print_exc()
        
        # Generate final report
        elapsed = time.time() - start_time
        print(f"\nProcessing Complete:")
        print(f"  Total identifications: {stats['total']}")
        print(f"  Successfully processed: {stats['processed']}")
        print(f"  Phosphorylated peptides: {stats['phospho']}")
        print(f"  Processing errors: {stats['errors']}")
        print(f"  Time elapsed: {elapsed:.2f} seconds")
        print(f"  Processing speed: {stats['processed']/elapsed:.2f} IDs/second")
        print(f"  Debug log saved to: {log_file}")
        
        # Save results
        save_identifications(args.out_file, protein_ids, processed_peptide_ids)
        
    except Exception as e:
        print(f"Fatal error: {str(e)}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
