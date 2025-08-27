import itertools
import math
import numpy as np
from pyopenms import (
    AASequence,
    ResidueModification,
    ModificationsDB,
    TheoreticalSpectrumGenerator,
    MSSpectrum,
    MSExperiment,
    PeptideIdentification,
    PeptideHit,
    Precursor,
    Constants
)
import sys # For checking float limits

# --- Configuration ---
MODIFICATION_NAME = "Phospho"
POTENTIAL_SITES = set(['S', 'T', 'Y'])
FRAGMENT_TOLERANCE = 0.05 # More typical tolerance for PhosphoRS scoring (Da)
FRAGMENT_METHOD_PPM = False # PhosphoRS typically uses Da tolerance
ADD_PRECURSOR_PEAK = False
ADD_ION_TYPES = ("b", "y")
MAX_ION_CHARGE = 2 # Adjust based on typical fragmentation

# --- Constants ---
LOG10_ZERO_REPLACEMENT = -100.0 # Value to use for log10(0) or log10(very small number)
MIN_PROBABILITY = 1e-10 # Minimum probability to avoid log10(0)

# --- Helper: Calculate Occurrence Probability 'p' ---
def get_occurrence_probability(exp_spectrum: MSSpectrum, tolerance_da: float) -> float:
    """
    Calculates the probability 'p' of matching a single theoretical peak
    to any experimental peak by chance.

    Args:
        exp_spectrum: The experimental spectrum (should be peak-picked).
        tolerance_da: The fragment tolerance in Daltons.

    Returns:
        The occurrence probability 'p'. Returns 0 if spectrum is empty or range is invalid.
    """
    peaks = exp_spectrum.get_peaks() # Get (mz, intensity) tuples
    num_peaks = len(peaks)
    if num_peaks == 0:
        return 0.0

    min_mz = exp_spectrum.getMinMZ()
    max_mz = exp_spectrum.getMaxMZ()
    
    # If min/max not stored, calculate from peaks
    if min_mz >= max_mz or min_mz <= 0:
        if num_peaks > 0:
            min_mz = peaks[0][0]
            max_mz = peaks[-1][0] # Assumes peaks are somewhat sorted, or get min/max explicitly
            # Explicit min/max calculation:
            # all_mz = [p[0] for p in peaks]
            # min_mz = min(all_mz)
            # max_mz = max(all_mz)
        else: # Should have been caught by num_peaks == 0
             return 0.0

    mz_range = max_mz - min_mz
    if mz_range <= 0:
        return 0.0 # Avoid division by zero

    # Formula: p = (NumPeaks * 2 * Tolerance) / Range
    # The '2 * tolerance' represents the total window width around a peak.
    occurrence_prob = (num_peaks * 2.0 * tolerance_da) / mz_range

    # Clamp probability between a minimum value and 1.0
    return max(MIN_PROBABILITY, min(occurrence_prob, 1.0 - MIN_PROBABILITY)) # Prevent p=0 or p=1 for logs

# --- Helper: PhosphoRS-like Scoring Function ---
def calculate_phosphors_score(
    theo_spectrum: MSSpectrum,
    exp_spectrum: MSSpectrum,
    occurrence_probability: float,
    tolerance: float,
    is_ppm: bool,
) -> float:
    """
    Calculates a PhosphoRS-like score based on binomial probability.

    Args:
        theo_spectrum: Theoretical spectrum.
        exp_spectrum: Experimental spectrum (should be peak-picked).
        occurrence_probability (p): Probability of random peak match.
        tolerance: Fragment tolerance.
        is_ppm: Whether tolerance is in ppm.

    Returns:
        PhosphoRS-like score (lower is better).
    """
    theo_peaks = sorted(theo_spectrum, key=lambda p: p.getMZ())
    exp_peaks = sorted(exp_spectrum, key=lambda p: p.getMZ()) # Ensure sorted

    if not theo_peaks or not exp_peaks:
        return -LOG10_ZERO_REPLACEMENT * len(theo_peaks) # Penalize heavily if no theo peaks

    p = occurrence_probability
    if p <= MIN_PROBABILITY or p >= 1.0 - MIN_PROBABILITY: # Handle edge cases from get_occurrence_probability
        log_p = LOG10_ZERO_REPLACEMENT if p <= MIN_PROBABILITY else 0.0
        log_1_minus_p = 0.0 if p <= MIN_PROBABILITY else LOG10_ZERO_REPLACEMENT
    else:
        log_p = math.log10(p)
        log_1_minus_p = math.log10(1.0 - p)
        
    score = 0.0
    exp_idx = 0
    matched_exp_indices = set() # Keep track of used experimental peaks

    # Iterate through theoretical peaks and find best match
    for theo_peak in theo_peaks:
        theo_mz = theo_peak.getMZ()

        # Calculate tolerance window in Da
        if is_ppm:
            tol_dalton = theo_mz * tolerance / 1_000_000.0
        else:
            tol_dalton = tolerance
        min_mz = theo_mz - tol_dalton
        max_mz = theo_mz + tol_dalton

        best_match_exp_idx = -1
        min_mass_diff = float('inf')

        # Advance experimental peak index pointer
        while exp_idx < len(exp_peaks) and exp_peaks[exp_idx].getMZ() < min_mz:
            exp_idx += 1
            
        current_exp_idx = exp_idx
        # Search within the tolerance window
        while current_exp_idx < len(exp_peaks) and exp_peaks[current_exp_idx].getMZ() <= max_mz:
            # Check if this exp peak is available and is a better match
            if current_exp_idx not in matched_exp_indices:
                mass_diff = abs(exp_peaks[current_exp_idx].getMZ() - theo_mz)
                if mass_diff < min_mass_diff:
                    min_mass_diff = mass_diff
                    best_match_exp_idx = current_exp_idx
            current_exp_idx += 1

        # Update score based on whether a match was found
        if best_match_exp_idx != -1:
            score += -10.0 * log_p  # Matched
            matched_exp_indices.add(best_match_exp_idx) # Mark experimental peak as used
        else:
            score += -10.0 * log_1_minus_p # Unmatched
            
        # Reset exp_idx pointer slightly for next theoretical peak if needed
        # This optimization depends on the assumption that theoretical peaks are mostly increasing
        # It might be safer to reset exp_idx = 0 for each theo_peak if theo peaks aren't strictly ordered
        # or have large mass gaps, but sorting theo_peaks helps.
        # Let's stick with the advancing pointer for efficiency.

    return score


# --- Main PhosphoRS-like Localization Function ---
def calculate_phospho_localization_compomics_style(
    peptide_hit: PeptideHit,
    spectrum: MSSpectrum,
    modification_name: str = MODIFICATION_NAME,
    potential_sites: set = POTENTIAL_SITES,
    fragment_tolerance: float = FRAGMENT_TOLERANCE,
    fragment_method_ppm: bool = FRAGMENT_METHOD_PPM,
    add_precursor_peak: bool = ADD_PRECURSOR_PEAK,
    add_ion_types: tuple = ADD_ION_TYPES,
    max_ion_charge: int = MAX_ION_CHARGE
):
    """
    Calculates phosphorylation site probabilities using a scoring method
    inspired by the Compomics PhosphoRS implementation.

    Args:
        peptide_hit: The peptide identification result.
        spectrum: The experimental MS/MS spectrum (MUST be peak-picked).
        modification_name: Name of the modification (must be in ModificationsDB).
        potential_sites: Set of single-letter amino acid codes for potential sites.
        fragment_tolerance: Fragment ion tolerance.
        fragment_method_ppm: True if tolerance is ppm, False if Da.
        add_precursor_peak: Whether to add precursor ion peak in theoretical spectra.
        add_ion_types: Ion types to generate (e.g., ('b', 'y')).
        max_ion_charge: Maximum charge for fragment ions.

    Returns:
        dict: A dictionary mapping 0-based residue index to its phosphorylation probability.
              Returns None if localization cannot be performed.
        list: A list of tuples, each containing (isomer_sequence_str, score) for debugging/info.
              Returns None if localization cannot be performed.

    Raises:
        ValueError: If modification name not found or spectrum is unsuitable.
        Exception: If precursor info is missing.
    """
    if fragment_method_ppm and fragment_tolerance <= 0:
         raise ValueError("PPM tolerance must be positive.")
    if not fragment_method_ppm and fragment_tolerance <= 0:
        raise ValueError("Dalton tolerance must be positive.")
        
    mod_db = ModificationsDB()
    mod = mod_db.getModification(modification_name)
    # Try finding by full name if UniMod style name was given
    if mod is None or mod.getName() == "unknown modification":
        mod = mod_db.findModification(modification_name, "", "")
    if mod is None or mod.getName() == "unknown modification":
         raise ValueError(f"Modification '{modification_name}' not found in ModificationsDB.")

    mod_mass = mod.getDiffMonoMass()
    target_residues = set(mod.getOrigin()) if mod.getOrigin() != 'X' else potential_sites
    if potential_sites:
        target_residues = target_residues.intersection(potential_sites)
    if not target_residues:
        print(f"Warning: Modification '{modification_name}' has no target residues specified or compatible with potential_sites.")
        return None, None

    original_sequence = peptide_hit.getSequence()
    if not isinstance(original_sequence, AASequence):
         original_sequence = AASequence.fromString(str(original_sequence))

    # --- Determine number of modifications to place ---
    num_mods_present = 0
    try:
        expected_mass_no_mods = AASequence.fromString(original_sequence.toUnmodifiedString()).getMonoWeight()
        expected_mass_with_mods = original_sequence.getMonoWeight()
        mass_diff = expected_mass_with_mods - expected_mass_no_mods
        
        if abs(mod_mass) > 1e-6:
            num_mods_present = round(mass_diff / mod_mass)
            if not math.isclose(num_mods_present * mod_mass, mass_diff, abs_tol=0.1):
                 print(f"Warning: Mass difference ({mass_diff:.4f}) does not cleanly match multiples of {modification_name} mass ({mod_mass:.4f}). Inferred {num_mods_present} mods. Checking string...")
                 # Fallback check based on string representation (less reliable)
                 mod_string_for_check = f"({modification_name})" 
                 num_mods_present_str = str(original_sequence).count(mod_string_for_check)
                 if num_mods_present_str != num_mods_present:
                      print(f"Warning: String count ({num_mods_present_str}) differs from mass-based count ({num_mods_present}). Using mass-based count.")
                      # Decide which one to trust - mass is usually better unless multiple mods interfere
                      # For now, stick with mass-based.
        else: # Modification mass is (near) zero - cannot use mass difference
             mod_string_for_check = f"({modification_name})" 
             num_mods_present = str(original_sequence).count(mod_string_for_check)
             print(f"Warning: Modification mass is zero. Counting mods based on string representation: {num_mods_present}")

    except Exception as e:
         print(f"Error determining number of modifications: {e}. Attempting string count.")
         try:
              mod_string_for_check = f"({modification_name})" 
              num_mods_present = str(original_sequence).count(mod_string_for_check)
         except Exception as e2:
              raise ValueError(f"Could not determine number of modifications in sequence '{original_sequence}' by mass or string: {e2}") from e
              
    if num_mods_present <= 0:
        print(f"Info: No modification '{modification_name}' found or inferred in sequence '{original_sequence}'. Cannot localize.")
        # Return 0 probability for all potential sites? Or None? Let's return None.
        return None, None

    unmodified_sequence_str = original_sequence.toUnmodifiedString()
    potential_site_indices = [
        i for i, aa in enumerate(unmodified_sequence_str) if aa in target_residues
    ]

    if len(potential_site_indices) < num_mods_present:
        print(f"Warning: Potential sites ({len(potential_site_indices)}) < mods to place ({num_mods_present}). Cannot localize.")
        return None, None
    if len(potential_site_indices) == num_mods_present:
         print(f"Info: Potential sites ({len(potential_site_indices)}) == mods to place ({num_mods_present}). Localization is trivial.")
         trivial_probs = {idx: 1.0 for idx in potential_site_indices}
         isomer_seq = AASequence.fromString(unmodified_sequence_str)
         for idx in potential_site_indices:
             isomer_seq.setModification(idx, modification_name)
         # Assign dummy score (e.g., 0, as lower is better)
         isomer_list = [(str(isomer_seq), 0.0)] 
         return trivial_probs, isomer_list

    # --- Pre-calculate for Scoring ---
    if spectrum.getMSLevel() != 2:
        print("Warning: Spectrum MSLevel is not 2. Results may be incorrect.")
    if len(spectrum) == 0:
         raise ValueError("Experimental spectrum contains no peaks.")
         
    # Ensure tolerance is in Da for occurrence probability calculation
    if fragment_method_ppm:
        # Use a representative M/Z for conversion, e.g., the precursor M/Z or midpoint
        # Using precursor is better if available
        prec_mz = spectrum.getPrecursors()[0].getMZ() if spectrum.getPrecursors() else np.mean(spectrum.get_peaks()[0]) if len(spectrum)>0 else 1000.0 # Guess if no precursor
        tolerance_da_for_p = prec_mz * fragment_tolerance / 1_000_000.0
    else:
        tolerance_da_for_p = fragment_tolerance
        
    occurrence_prob = get_occurrence_probability(spectrum, tolerance_da_for_p)
    if occurrence_prob <= MIN_PROBABILITY:
         print("Warning: Occurrence probability 'p' is extremely low. Scoring might be unstable.")
         # Maybe raise an error or return None? For now, proceed with warning.

    # --- Generate Isomers and Theoretical Spectra ---
    spec_gen = TheoreticalSpectrumGenerator()
    params = spec_gen.getParameters()
    params.setValue("add_metainfo", "true")
    params.setValue("add_precursor_peaks", "true" if add_precursor_peak else "false")
    params.setValue("max_charge", max_ion_charge)
    all_possible_ions = ["a", "b", "c", "x", "y", "z"]
    for ion_type in all_possible_ions:
        params.setValue(f"add_{ion_type}_ions", "true" if ion_type in add_ion_types else "false")
    # Consider adding neutral loss options if relevant for Phospho
    # params.setValue("add_losses", "true")
    # params.setValue("add_phospho_loss", "true") # Check actual param names in pyopenms doc
    spec_gen.setParameters(params)

    isomer_scores = [] # List to store {"isomer_seq": AASeq, "score": float, "sites": set}

    if not spectrum.getPrecursors():
        raise Exception("Experimental spectrum must have precursor information (charge).")
    precursor_charge = spectrum.getPrecursors()[0].getCharge()

    for site_indices_tuple in itertools.combinations(potential_site_indices, num_mods_present):
        isomer_seq = AASequence.fromString(unmodified_sequence_str)
        for site_index in site_indices_tuple:
            isomer_seq.setModification(site_index, modification_name)

        theo_spectrum = MSSpectrum()
        # Ensure spectrum generation uses correct charge state(s) based on precursor
        spec_gen.getSpectrum(theo_spectrum, isomer_seq, 1, precursor_charge) # Generates ions up to precursor_charge

        # Score against experimental spectrum using PhosphoRS-like method
        score = calculate_phosphors_score(
            theo_spectrum, spectrum, occurrence_prob, fragment_tolerance, fragment_method_ppm
        )

        isomer_scores.append({"isomer_seq": isomer_seq, "score": score, "sites": set(site_indices_tuple)})

    # --- Calculate Probabilities using PhosphoRS normalization ---
    if not isomer_scores:
        print("Warning: No isomers were generated or scored.")
        return None, None

    min_score = min(item["score"] for item in isomer_scores) if isomer_scores else 0.0

    total_prob_contribution = 0.0
    for item in isomer_scores:
        # score_diff is non-positive. If score == min_score, score_diff = 0.
        score_diff = item["score"] - min_score 
        # Calculate 10^(-score_diff / 10). Handles potential numerical issues better.
        prob_contrib = 10.0 ** (-score_diff / 10.0) 
        item["prob_contrib"] = prob_contrib
        total_prob_contribution += prob_contrib

    isomer_probabilities = []
    if total_prob_contribution > 1e-9: # Avoid division by zero
        for item in isomer_scores:
            prob = item["prob_contrib"] / total_prob_contribution
            item["probability"] = prob
            isomer_probabilities.append(prob)
    else: # All scores might be extremely high (bad) or identical
        num_isomers = len(isomer_scores)
        prob = 1.0 / num_isomers if num_isomers > 0 else 0.0
        for item in isomer_scores:
            item["probability"] = prob
        isomer_probabilities = [prob] * num_isomers


    # Calculate site probabilities
    site_probabilities = {site_idx: 0.0 for site_idx in potential_site_indices}
    # Also track probabilities for non-potential sites (should remain 0)
    # for i in range(len(unmodified_sequence_str)):
    #      if i not in site_probabilities:
    #           site_probabilities[i] = 0.0

    for i, item in enumerate(isomer_scores):
        prob = item.get("probability", 0.0)
        for site_idx in item["sites"]:
            if site_idx in site_probabilities:
                site_probabilities[site_idx] += prob

    # Format isomer list for return
    isomer_list_out = [(str(item["isomer_seq"]), item["score"]) for item in isomer_scores]

    return site_probabilities, isomer_list_out


# --- Example Usage (Similar to before, but using new function) ---
if __name__ == "__main__":
    # --- 1. Create Dummy Data (Replace with loading real data) ---
    # Using the same peptide as before: TESTPEPTIDESEK
    peptide_sequence_str_unmod = "TESTPEPTIDESEK"
    target_phospho_sites_indices = [1, 3, 9] # Indices for S(1), T(3), S(9)
    num_phospho_to_place = 1

    # Create AASequence with one phospho placed initially (e.g., on S(1))
    initial_seq = AASequence.fromString(peptide_sequence_str_unmod)
    initial_seq.setModification(target_phospho_sites_indices[0], MODIFICATION_NAME)
    
    ph = PeptideHit()
    ph.setSequence(initial_seq)
    ph.setCharge(2)
    ph.setScore(50.0) # Dummy search engine score

    # Example Experimental Spectrum (MUST be peak-picked!)
    exp_spec = MSSpectrum()
    exp_spec.setMSLevel(2)
    exp_spec.setPrecursors([Precursor()])
    # Calculate precursor m/z correctly
    prec_mz = (initial_seq.getMonoWeight() + (ph.getCharge() * Constants.PROTON_MASS_U)) / ph.getCharge()
    exp_spec.getPrecursors()[0].setMZ(prec_mz)
    exp_spec.getPrecursors()[0].setCharge(ph.getCharge())

    # Add peaks favoring T(3) phosphorylation (same peaks as before)
    # m/z values might need slight adjustment depending on exact mass calculation used by pyopenms
    exp_spec.addPeak(MSSpectrum.Peak1D(mz=147.11, intensity=50.0))  # y1 (K)
    exp_spec.addPeak(MSSpectrum.Peak1D(mz=276.15, intensity=40.0))  # y2 (EK)
    exp_spec.addPeak(MSSpectrum.Peak1D(mz=363.18, intensity=30.0))  # y3 (SEK)
    exp_spec.addPeak(MSSpectrum.Peak1D(mz=499.13, intensity=90.0))  # b4 with phospho on T(3)
    exp_spec.addPeak(MSSpectrum.Peak1D(mz=591.30, intensity=70.0))  # y5 (IDESEK)
    exp_spec.addPeak(MSSpectrum.Peak1D(mz=772.31, intensity=100.0)) # y8 with phospho on T(3)
    exp_spec.addPeak(MSSpectrum.Peak1D(mz=990.33, intensity=45.0))  # b9 with phospho on S(9) (charge 1)
    # Add some noise peaks
    exp_spec.addPeak(MSSpectrum.Peak1D(mz=300.1, intensity=10.0))
    exp_spec.addPeak(MSSpectrum.Peak1D(mz=550.2, intensity=15.0))
    exp_spec.addPeak(MSSpectrum.Peak1D(mz=800.3, intensity=12.0))

    # IMPORTANT: Ensure spectrum is sorted by m/z for efficient matching
    exp_spec.sortByPosition()


    # --- 2. Run Localization ---
    print(f"Running Compomics-style localization for: {ph.getSequence()}")
    print(f"Number of mods to place: {num_phospho_to_place}")
    print(f"Potential site indices (0-based): {target_phospho_sites_indices}")
    print(f"Fragment Tolerance: {FRAGMENT_TOLERANCE} {'ppm' if FRAGMENT_METHOD_PPM else 'Da'}")
    print("-" * 20)
    
    try:
        site_probs, isomer_details = calculate_phospho_localization_compomics_style(
            ph,
            exp_spec,
            fragment_tolerance=FRAGMENT_TOLERANCE, # Use global config
            fragment_method_ppm=FRAGMENT_METHOD_PPM, # Use global config
            max_ion_charge=ph.getCharge()
        )

        # --- 3. Print Results ---
        if site_probs is not None:
            print("Isomer Scores (Lower is better):")
            # Sort by score (ascending)
            for seq_str, score in sorted(isomer_details, key=lambda x: x[1]):
                 print(f"  - {seq_str}: {score:.4f}")

            print("\nSite Probabilities:")
            # Sort by index for clarity
            sorted_sites = sorted(site_probs.items())
            total_prob = 0.0
            for site_index, probability in sorted_sites:
                 aa = peptide_sequence_str_unmod[site_index]
                 print(f"  - Site {site_index} ({aa}): {probability:.4f}")
                 total_prob += probability
            print(f"Total Probability Sum: {total_prob:.4f}") # Should be close to 1.0

            # Construct sequence string with probabilities (example format)
            result_seq = ""
            max_prob_site = -1
            max_prob = -1.0
            sites_with_prob = set(site_probs.keys())
            for i, aa in enumerate(peptide_sequence_str_unmod):
                 result_seq += aa
                 if i in sites_with_prob and site_probs[i] > 0.0001: # Only annotate non-zero probs
                     prob = site_probs[i]
                     result_seq += f"({MODIFICATION_NAME};{prob:.2f})" # Use semicolon for clarity
                     if prob > max_prob:
                          max_prob = prob
                          max_prob_site = i
            print(f"\nSequence with probabilities: {result_seq}")
            print(f"(Highest probability site: {max_prob_site} ({peptide_sequence_str_unmod[max_prob_site]}) with P={max_prob:.4f})")

        else:
            print("Localization could not be performed.")
            
    except ValueError as e:
         print(f"ValueError during localization: {e}")
    except Exception as e:
         print(f"An unexpected error occurred: {e}")