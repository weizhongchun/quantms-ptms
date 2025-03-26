from pyopenms import *
import math
import numpy as np
from collections import defaultdict

class AScore(DefaultParamHandler):
    """
    Implementation of the Ascore algorithm for phosphorylation site localization.
    
    For a given peptide sequence and its MS/MS spectrum, it identifies the most 
    probable phosphorylation site(s). For each phosphorylation site, a probability 
    score is calculated.
    
    The algorithm is implemented according to Beausoleil et al. (Nat. Biotechnol. 2006):
    "A probability-based approach for high-throughput protein phosphorylation 
    analysis and site localization"
    """
    
    def __init__(self):
        """Initialize the AScore object with default parameters."""
        DefaultParamHandler.__init__(self, "AScore")
        
        # Set default parameters
        defaults = Param()
        defaults.setValue("fragment_mass_tolerance", 0.05, "Fragment mass tolerance for spectrum comparisons")
        defaults.setMinFloat("fragment_mass_tolerance", 0.0)
        
        defaults.setValue("fragment_mass_unit", "Da", "Unit of fragment mass tolerance")
        defaults.setValidStrings("fragment_mass_unit", ["Da", "ppm"])
        
        advanced = ["advanced"]  # tag for advanced parameters
        
        defaults.setValue("max_peptide_length", 40, 
                         "Restrict scoring to peptides with a length no greater than this value ('0' for 'no restriction')", 
                         advanced)
        defaults.setMinInt("max_peptide_length", 0)
        
        defaults.setValue("max_num_perm", 16384, 
                         "Maximum number of permutations a sequence can have to be processed ('0' for 'no restriction')", 
                         advanced)
        defaults.setMinInt("max_num_perm", 0)
        
        defaults.setValue("add_decoy", "false", 
                         "Include PhosphoDecoy site (A) in phosphorylation site analysis for FLR calculation", 
                         advanced)
        defaults.setValidStrings("add_decoy", ["true", "false"])
        
        self.setDefaults(defaults)
        
        # Initialize member variables
        self.fragment_mass_tolerance_ = 0.05
        self.fragment_tolerance_ppm_ = False
        self.max_peptide_length_ = 40
        self.max_permutations_ = 16384
        self.add_decoys_ = False
        self.base_match_probability_ = 0.0
        self.unambiguous_score_ = 1e6
        
        # Initialize TheoreticalSpectrumGenerator
        self.spectrum_generator_ = TheoreticalSpectrumGenerator()
        p = self.spectrum_generator_.getParameters()
        p.setValue("isotope_model", "none")
        p.setValue("add_first_prefix_ion", "true")
        self.spectrum_generator_.setParameters(p)
        
        self.updateMembers_()
    
    def updateMembers_(self):
        """Update member variables from the parameter object."""
        self.fragment_mass_tolerance_ = self.getParameter("fragment_mass_tolerance").toDouble()
        self.fragment_tolerance_ppm_ = (self.getParameter("fragment_mass_unit").toString() == "ppm")
        self.max_peptide_length_ = self.getParameter("max_peptide_length").toInt()
        self.max_permutations_ = self.getParameter("max_num_perm").toInt()
        self.add_decoys_ = self.getParameter("add_decoy").toBool()
    
    @staticmethod
    def isPhosphoDecoySite(residue):
        """Check if a residue is a PhosphoDecoy site (A)."""
        return residue == 'A'
    
    @staticmethod
    def isPhosphoSite(residue):
        """Check if a residue is a standard phosphorylation site (S, T, Y)."""
        return residue in ['S', 'T', 'Y']
    
    def compute(self, hit, real_spectrum):
        """
        Compute the AScore and return a modified PeptideHit with phosphorylation site scores.
        
        Parameters:
        -----------
        hit : PeptideHit
            The peptide hit containing the sequence with potential phosphorylation sites
        real_spectrum : MSSpectrum
            The experimental MS/MS spectrum
            
        Returns:
        --------
        PeptideHit
            A modified PeptideHit with AScore values for each phosphorylation site
            The score field contains the best AScore value
        """
        phospho = PeptideHit(hit)
        
        # reset phospho
        phospho.setScore(-1)
        phospho.setMetaValue("search_engine_sequence", hit.getSequence().toString())
        
        # Early termination for empty spectra
        if real_spectrum.size() == 0:
            print("Warning: Empty spectrum provided to AScore::compute. Returning original hit.")
            return phospho
        
        sequence_str = phospho.getSequence().toString()
        unmodified_sequence_str = phospho.getSequence().toUnmodifiedString()
        
        # Count phosphorylation events and remove phosphorylations to get the base peptide sequence
        number_of_phosphorylation_events = self.numberOfPhosphoEvents_(sequence_str)
        seq_without_phospho = self.removePhosphositesFromSequence_(sequence_str)
        
        # Initialize counters for regular and decoy phosphorylation events
        regular_phospho_count = sequence_str.count("(Phospho)")
        decoy_phospho_count = sequence_str.count("(PhosphoDecoy)")
        
        if self.add_decoys_:
            print(f"Debug: Found {regular_phospho_count} regular phosphorylation events and "
                 f"{decoy_phospho_count} PhosphoDecoy events in sequence: {sequence_str}")
        elif decoy_phospho_count > 0:
            print("Warning: PhosphoDecoy sites found in sequence, but add_decoys is set to false. "
                 "Please enable add_decoys to include PhosphoDecoy sites in the analysis. "
                 "Returning original hit.")
            return phospho
        
        # Check if peptide is too long for analysis
        if (self.max_peptide_length_ > 0) and (len(unmodified_sequence_str) > self.max_peptide_length_):
            print(f"Debug: Calculation aborted: peptide too long: {seq_without_phospho.toString()}")
            return phospho
        
        # Get all potential phosphorylation sites (S,T,Y and optionally A for decoys)
        sites = self.getSites_(unmodified_sequence_str)
        number_of_sites = len(sites)
        
        if number_of_phosphorylation_events == 0 or number_of_sites == 0:
            return phospho
        
        # Add metadata about phosphorylation types
        if self.add_decoys_:
            phospho.setMetaValue("regular_phospho_count", int(regular_phospho_count))
            phospho.setMetaValue("phospho_decoy_count", int(decoy_phospho_count))
        phospho.setMetaValue("phospho_sites", int(number_of_phosphorylation_events))
        
        # Special case: If all possible sites are phosphorylated, the localization is unambiguous
        if number_of_sites == number_of_phosphorylation_events:
            phospho.setScore(self.unambiguous_score_)
            
            # Create a map with scores of 1.0 for all phosphorylation sites
            site_scores = {}
            for site in sites:
                site_scores[site] = self.unambiguous_score_
            
            # Create ProForma string with scores for all sites
            proforma = self.generateProFormaString_(phospho.getSequence(), site_scores)
            phospho.setMetaValue("ProForma", proforma)
            phospho.setMetaValue("AScore_pep_score", self.unambiguous_score_)
            return phospho
        
        # Generate all possible permutations of phosphorylation sites
        permutations = self.computePermutations_(sites, number_of_phosphorylation_events)
        print(f"Debug: Number of permutations: {len(permutations)}")
        
        # Check if permutations is empty or exceeds the maximum allowed
        if not permutations or ((self.max_permutations_ > 0) and (len(permutations) > self.max_permutations_)):
            print("Debug: Calculation aborted: number of permutations exceeded or early termination")
            return phospho
        
        # Create theoretical spectra for all permutations
        th_spectra = self.createTheoreticalSpectra_(permutations, seq_without_phospho)
        
        # Prepare real spectrum windows
        if not real_spectrum.isSorted():
            real_spectrum.sortByPosition()
        windows_top10 = self.peakPickingPerWindowsInSpectrum_(real_spectrum)
        
        # Compute match probability for a peak depth of 1
        self.base_match_probability_ = self.computeBaseProbability_(real_spectrum[real_spectrum.size()-1].getMZ())
        
        # Calculate peptide score for each possible phospho site permutation
        peptide_site_scores = self.calculatePermutationPeptideScores_(th_spectra, windows_top10)
        
        # Sort peptide permutations by (weighted) peptide score
        ranking = self.rankWeightedPermutationPeptideScores_(peptide_site_scores)
        
        # Get the best scoring permutation and set it as the peptide sequence
        sorted_ranking = sorted(ranking.items(), key=lambda x: x[0], reverse=True)
        
        # Note: there might be multiple permutations with the highest score. If the original
        # sequence is one of them, we take this one. Otherwise, we take the first one.
        best_score = sorted_ranking[0][0]
        best_permutations = [item for item in sorted_ranking if item[0] == best_score]
        
        # Check if original sequence is among the best scoring permutations
        rev = 0  # Index of the best permutation
        for i, (score, idx) in enumerate(best_permutations):
            # Get sequence of this best-scoring permutation
            best_sequence = th_spectra[idx].getName()
            # If the sequence is the same as the original sequence choose this one
            if best_sequence == sequence_str:
                rev = i
                break
        
        # The name of the spectrum contains the peptide sequence
        seq1 = th_spectra[sorted_ranking[rev][1]].getName()
        phospho.setSequence(AASequence.fromString(seq1))
        
        peptide1_score = sorted_ranking[rev][0]
        phospho.setMetaValue("AScore_pep_score", peptide1_score)
        
        seq2 = th_spectra[sorted_ranking[rev+1][1]].getName()
        peptide2_score = sorted_ranking[rev+1][0]
        
        # Determine the highest scoring permutations for each phosphorylation site
        phospho_sites = []
        self.determineHighestScoringPermutations_(peptide_site_scores, phospho_sites, permutations, sorted_ranking)
        
        # Calculate AScore for each phosphorylation site
        rank = 1
        best_Ascore = float('inf')  # The lower the better
        site2score = {}  # Map to store the scores for each phospho site
        
        for phospho_site in phospho_sites:
            Ascore = 0.0
            if peptide1_score == peptide2_score:  # Set Ascore = 0 for each phosphorylation site
                print(f"Debug: Score of best ({seq1}) and second best peptide ({seq2}) are equal ({peptide1_score})")
            else:
                # Calculate site-determining ions
                site_determining_ions = []
                self.computeSiteDeterminingIons_(th_spectra, phospho_site, site_determining_ions)
                N = site_determining_ions[0].size()  # All possibilities have same number, take first one
                p = float(phospho_site['peak_depth']) * self.base_match_probability_
                
                n_first = 0  # Number of matching peaks for first peptide
                # Count matched ions across all windows for the first permutation
                for window_idx in range(len(windows_top10)):
                    n_first += self.numberOfMatchedIons_(site_determining_ions[0], windows_top10[window_idx], phospho_site['peak_depth'])
                
                P_first = self.computeCumulativeScore_(N, n_first, p)
                
                n_second = 0  # Number of matching peaks for second peptide
                # Count matched ions across all windows for the second permutation
                for window_idx in range(len(windows_top10)):
                    n_second += self.numberOfMatchedIons_(site_determining_ions[1], windows_top10[window_idx], phospho_site['peak_depth'])
                
                N2 = site_determining_ions[1].size()
                P_second = self.computeCumulativeScore_(N2, n_second, p)
                
                # Convert probabilities to scores: -10 * log10(P)
                # abs is used to avoid -0 score values
                score_first = abs(-10.0 * math.log10(P_first))
                score_second = abs(-10.0 * math.log10(P_second))
                
                print(f"Debug: first - N: {N}, p: {p}, n: {n_first}, score: {score_first}")
                print(f"Debug: second - N: {N2}, p: {p}, n: {n_second}, score: {score_second}")
                
                # AScore is the difference between the two scores
                Ascore = score_first - score_second
                print(f"Debug: Ascore_{rank}: {Ascore}")
            
            if Ascore < best_Ascore:
                best_Ascore = Ascore
                
            phospho.setMetaValue(f"AScore_{rank}", Ascore)
            site2score[phospho_site['first']] = Ascore
            rank += 1
        
        # Generate a ProForma-like string with phosphorylation site scores
        proforma = self.generateProFormaString_(phospho.getSequence(), site2score)
        phospho.setMetaValue("ProForma", proforma)
        phospho.setScore(best_Ascore)
        
        return phospho
    
    def computeBaseProbability_(self, ppm_reference_mz):
        """Compute the base probability of a random peak match."""
        base_match_probability = 2.0 * self.fragment_mass_tolerance_ / 100.0
        if self.fragment_tolerance_ppm_:
            base_match_probability *= ppm_reference_mz * 1e-6  # Convert to ppm
        return base_match_probability
    
    def computeCumulativeScore_(self, N, n, p):
        """
        Compute the cumulative binomial probability P(X ≥ n).
        
        Parameters:
        -----------
        N : int
            Total number of trials (theoretical peaks)
        n : int
            Number of successes (matched peaks)
        p : float
            Probability of success for a single trial
            
        Returns:
        --------
        float
            Cumulative binomial probability P(X ≥ n)
        """
        assert n <= N, "The number of matched ions (n) can be at most as large as the number of trials (N)."
        assert 0 <= p <= 1.0, "p must be a probability [0,1]."
        
        # Use scipy's implementation for binomial CDF complement if available
        try:
            from scipy.stats import binom
            return 1.0 - binom.cdf(n-1, N, p)
        except ImportError:
            # Fallback implementation: Calculate P(X ≥ n) directly
            # This is less numerically stable but works without scipy
            total_prob = 0.0
            for k in range(n, N+1):
                # Calculate binomial coefficient * p^k * (1-p)^(N-k)
                coef = 1.0
                for i in range(1, k+1):
                    coef *= (N - k + i) / i
                prob = coef * (p**k) * ((1-p)**(N-k))
                total_prob += prob
            return total_prob
    
    def determineHighestScoringPermutations_(self, peptide_site_scores, sites, permutations, ranking):
        """
        Determine the highest scoring permutations for each phosphorylation site.
        
        For each phosphorylation site in the highest-scoring permutation, this function:
        1. Finds the next best permutation where this site is not phosphorylated
        2. Determines the peak depth that maximizes the score difference between these permutations
        3. Stores this information for AScore calculation
        """
        sites.clear()
        
        # Take first set of phospho site assignments
        best_peptide_sites = permutations[ranking[0][1]]  # Sites of the assignment that achieved the highest weighted score
        
        for i in range(len(best_peptide_sites)):  # For each phosphorylated site
            site_info = {'first': best_peptide_sites[i], 'seq_1': ranking[0][1]}
            continue_search = True
            
            rank_index = 1  # Start with the next-best scoring peptide
            # Iterate from best scoring peptide to the first peptide that doesn't contain the current phospho site
            while continue_search and rank_index < len(ranking):
                current_permutation = permutations[ranking[rank_index][1]]
                
                # Check if this permutation has each of the phosphorylation sites as the best one,
                # except for the current site (i) which should be different
                for j in range(len(best_peptide_sites)):
                    if j == i:  # The site we are interested in
                        if best_peptide_sites[j] in current_permutation:
                            continue_search = True  # Permutation is also phosphorylated -> skip
                            break
                        else:
                            continue_search = False
                    else:  # The other sites
                        if best_peptide_sites[j] not in current_permutation:
                            continue_search = True  # Phosphorylation state should be the same but isn't
                            break
                        else:
                            continue_search = False
                
                if continue_search:
                    rank_index += 1
            
            # Store permutation of peptide without the phospho site i (seq_2)
            if rank_index < len(ranking):
                site_info['seq_2'] = ranking[rank_index][1]
                
                # Store phospho site location that is not contained in seq_1 but in seq_2
                for position in permutations[site_info['seq_2']]:
                    if position not in best_peptide_sites:
                        site_info['second'] = position
                        break
                
                # Find the peak depth that maximizes the score difference
                maximum_score_difference = 0.0
                site_info['peak_depth'] = 1
                
                for depth in range(1, len(peptide_site_scores[site_info['seq_1']])+1):
                    phospho_at_site_score = peptide_site_scores[site_info['seq_1']][depth-1]
                    no_phospho_at_site_score = peptide_site_scores[site_info['seq_2']][depth-1]
                    score_difference = phospho_at_site_score - no_phospho_at_site_score
                    
                    if score_difference > maximum_score_difference:
                        maximum_score_difference = score_difference
                        site_info['peak_depth'] = depth
                        site_info['ascore'] = maximum_score_difference
                
                sites.append(site_info)
    
    def computeSiteDeterminingIons_(self, th_spectra, candidates, site_determining_ions):
        """
        Compute the site-determining ions for phosphorylation site candidates.
        
        Site-determining ions are fragment ions that can distinguish between different
        phosphorylation site localizations. This method identifies ions that are unique
        to each of the two best-scoring permutations for a given phosphorylation site.
        """
        site_determining_ions.clear()
        site_determining_ions.extend([MSSpectrum(), MSSpectrum()])
        
        spectrum_first = th_spectra[candidates['seq_1']]
        spectrum_second = th_spectra[candidates['seq_2']]
        
        # Find peaks that are unique to the first spectrum (not in the second)
        spectrum_first_diff = MSSpectrum()
        self.getSpectrumDifference_(spectrum_first, spectrum_second, spectrum_first_diff)
        
        # Find peaks that are unique to the second spectrum (not in the first)
        spectrum_second_diff = MSSpectrum()
        self.getSpectrumDifference_(spectrum_second, spectrum_first, spectrum_second_diff)
        
        site_determining_ions[0] = spectrum_first_diff
        site_determining_ions[1] = spectrum_second_diff
        
        site_determining_ions[0].sortByPosition()
        site_determining_ions[1].sortByPosition()
    
    def getSpectrumDifference_(self, spectrum1, spectrum2, result):
        """
        Find the difference between two spectra based on m/z values.
        
        This method finds peaks that are in spectrum1 but not in spectrum2,
        within the fragment mass tolerance.
        """
        result.clear()
        
        i, j = 0, 0
        while i < spectrum1.size() and j < spectrum2.size():
            mz1 = spectrum1[i].getMZ()
            mz2 = spectrum2[j].getMZ()
            val = self.compareMZ_(mz1, mz2)
            
            if val == -1:  # mz1 < mz2 (outside tolerance)
                result.push_back(spectrum1[i])
                i += 1
            elif val == 1:  # mz1 > mz2 (outside tolerance)
                j += 1
            else:  # within tolerance, check if more ions are within tolerance
                # Check mz2 until no match
                j += 1
                while j < spectrum2.size() and self.compareMZ_(mz1, spectrum2[j].getMZ()) == 0:
                    j += 1
                
                # Check mz1 until no match
                i += 1
                while i < spectrum1.size() and self.compareMZ_(spectrum1[i].getMZ(), mz2) == 0:
                    i += 1
        
        # Add remaining peaks from spectrum1
        while i < spectrum1.size():
            result.push_back(spectrum1[i])
            i += 1
    
    def numberOfMatchedIons_(self, th, window, depth):
        """
        Count the number of matched ions between a theoretical spectrum and experimental spectrum window.
        
        Parameters:
        -----------
        th : MSSpectrum
            Theoretical spectrum (must be sorted by m/z)
        window : MSSpectrum
            Experimental spectrum window (must be sorted by m/z)
        depth : int
            Maximum number of peaks to consider from the window (peak depth)
            
        Returns:
        --------
        int
            Number of matched ions
        """
        window_reduced = MSSpectrum(window)
        if window_reduced.size() > depth:
            # Resize by creating a new spectrum with only the first 'depth' peaks
            temp = MSSpectrum()
            for i in range(depth):
                if i < window_reduced.size():
                    temp.push_back(window_reduced[i])
            window_reduced = temp
        
        window_reduced.sortByPosition()
        matched_peaks = 0
        
        # Manually implement the matching logic
        for i in range(th.size()):
            th_mz = th[i].getMZ()
            for j in range(window_reduced.size()):
                window_mz = window_reduced[j].getMZ()
                # Check if within tolerance
                if self.compareMZ_(th_mz, window_mz) == 0:
                    matched_peaks += 1
                    break
        
        return matched_peaks
    
    def compareMZ_(self, mz1, mz2):
        """
        Compare two m/z values using the fragment mass tolerance.
        
        Returns:
        --------
        int
            -1 if mz1 < mz2 (outside tolerance)
            1 if mz1 > mz2 (outside tolerance)
            0 if they are within tolerance
        """
        tolerance = self.fragment_mass_tolerance_
        error = mz1 - mz2
        
        if self.fragment_tolerance_ppm_:
            avg_mass = (mz1 + mz2) / 2.0
            tolerance = tolerance * avg_mass / 1.0e6
        
        if error < -tolerance:
            return -1
        elif error > tolerance:
            return 1
        else:
            return 0
    
    def peptideScore_(self, scores):
        """
        Compute the weighted peptide score according to Beausoleil et al.
        
        Calculates a weighted average of scores at different peak depths according to
        the formula described in Beausoleil et al. (page 1291).
        
        Parameters:
        -----------
        scores : list
            List of 10 scores at different peak depths
            
        Returns:
        --------
        float
            Weighted peptide score
        """
        assert len(scores) == 10, "Scores vector must contain a score for every peak level."
        
        return (scores[0] * 0.5 +
                scores[1] * 0.75 +
                scores[2] +
                scores[3] +
                scores[4] +
                scores[5] +
                scores[6] * 0.75 +
                scores[7] * 0.5 +
                scores[8] * 0.25 +
                scores[9] * 0.25) / 7.0
    
    def getSites_(self, unmodified):
        """
        Get all potential phosphorylation sites in a peptide sequence.
        
        Parameters:
        -----------
        unmodified : str
            The unmodified peptide sequence
            
        Returns:
        --------
        list
            Positions (0-based) of potential phosphorylation sites
        """
        # Get both standard phosphorylation sites and PhosphoDecoy sites
        phospho_sites = self.getPhosphoSites_(unmodified)
        decoy_sites = []
        
        if self.add_decoys_:
            decoy_sites = self.getPhosphoDecoySites_(unmodified)
        
        # Combine both types of sites
        all_sites = phospho_sites + decoy_sites
        all_sites.sort()
        
        return all_sites
    
    def getPhosphoSites_(self, unmodified):
        """
        Identify standard phosphorylation sites (S, T, Y) in a peptide sequence.
        
        Parameters:
        -----------
        unmodified : str
            The unmodified peptide sequence
            
        Returns:
        --------
        list
            Positions (0-based) of standard phosphorylation sites
        """
        ret = []
        for i in range(len(unmodified)):
            if self.isPhosphoSite(unmodified[i]):
                ret.append(i)
        return ret
    
    def getPhosphoDecoySites_(self, unmodified):
        """
        Identify PhosphoDecoy sites (A) in a peptide sequence.
        
        Parameters:
        -----------
        unmodified : str
            The unmodified peptide sequence
            
        Returns:
        --------
        list
            Positions (0-based) of PhosphoDecoy sites
        """
        ret = []
        for i in range(len(unmodified)):
            if self.isPhosphoDecoySite(unmodified[i]):
                ret.append(i)
        return ret
    
    def numberOfPhosphoEvents_(self, sequence):
        """
        Count the number of phosphorylation events in a peptide sequence.
        
        Parameters:
        -----------
        sequence : str
            The peptide sequence as a string
            
        Returns:
        --------
        int
            Number of phosphorylation events
        """
        cnt_phospho_events = 0
        
        # Count regular phosphorylation events
        pos = 0
        while True:
            pos = sequence.find("(Phospho)", pos)
            if pos == -1:
                break
            cnt_phospho_events += 1
            pos += 9  # Move past "(Phospho)"
        
        # Count PhosphoDecoy events
        pos = 0
        while True:
            pos = sequence.find("(PhosphoDecoy)", pos)
            if pos == -1:
                break
            cnt_phospho_events += 1
            pos += 14  # Move past "(PhosphoDecoy)"
        
        print(f"Debug: Found {cnt_phospho_events} phosphorylation events in sequence: {sequence}")
        
        return cnt_phospho_events
    
    def computePermutations_(self, sites, n_phosphorylation_events):
        """
        Generate all possible combinations of phosphorylation sites.
        
        Parameters:
        -----------
        sites : list
            Positions of potential phosphorylation sites
        n_phosphorylation_events : int
            Number of phosphorylation events to place
            
        Returns:
        --------
        list
            All possible permutations of phosphorylation sites
        """
        permutations = []
        
        # Early termination if we can estimate that the number of permutations will exceed the maximum
        if self.max_permutations_ > 0 and n_phosphorylation_events >= 1:
            # Check if k > n (more phosphorylation events than sites)
            if n_phosphorylation_events > len(sites):
                print(f"Debug: Early termination: more phosphorylation events ({n_phosphorylation_events}) "
                     f"than available sites ({len(sites)})")
                return permutations
            else:
                # Estimate the number of permutations using the binomial coefficient (n choose k)
                try:
                    from math import comb
                    estimated_permutations = comb(len(sites), n_phosphorylation_events)
                except (ImportError, OverflowError):
                    # Fallback calculation
                    try:
                        estimated_permutations = 1
                        for i in range(1, n_phosphorylation_events + 1):
                            estimated_permutations *= (len(sites) - i + 1) / i
                    except OverflowError:
                        estimated_permutations = float('inf')
                
                if estimated_permutations > self.max_permutations_:
                    print(f"Debug: Early termination: estimated permutations ({estimated_permutations}) "
                         f"exceeds maximum ({self.max_permutations_})")
                    return permutations
        
        if n_phosphorylation_events == 0:
            return permutations
        elif n_phosphorylation_events == 1:
            for i in range(len(sites)):
                permutations.append([sites[i]])
            return permutations
        # All sites are phosphorylated? Return one permutation containing all sites at once.
        elif len(sites) == n_phosphorylation_events:
            permutations.append(list(sites))
            return permutations
        else:
            # Generate all n_phosphorylation_events sized combinations from sites
            # Use itertools combinations if available
            try:
                from itertools import combinations
                for combo in combinations(sites, n_phosphorylation_events):
                    permutations.append(list(combo))
                    # Early termination check
                    if self.max_permutations_ > 0 and len(permutations) > self.max_permutations_:
                        print(f"Debug: Early termination during iteration: current permutations ({len(permutations)}) "
                             f"exceeds maximum ({self.max_permutations_})")
                        return []
                return permutations
            except ImportError:
                # Fallback implementation using a recursive algorithm
                def backtrack(start, current_combo):
                    if len(current_combo) == n_phosphorylation_events:
                        permutations.append(list(current_combo))
                        return
                    
                    for i in range(start, len(sites)):
                        current_combo.append(sites[i])
                        backtrack(i + 1, current_combo)
                        current_combo.pop()
                        
                        # Early termination check
                        if self.max_permutations_ > 0 and len(permutations) > self.max_permutations_:
                            return
                
                backtrack(0, [])
                
                # Check if we terminated early
                if self.max_permutations_ > 0 and len(permutations) > self.max_permutations_:
                    return []
                    
                return permutations
    
    def removePhosphositesFromSequence_(self, sequence):
        """
        Create a variant of the peptide with all phosphorylations removed.
        
        Parameters:
        -----------
        sequence : str
            The peptide sequence as a string
            
        Returns:
        --------
        AASequence
            Base peptide sequence without phosphorylations
        """
        seq = sequence
        # Remove both regular and decoy phosphorylation markers
        seq = seq.replace("(Phospho)", "")
        seq = seq.replace("(PhosphoDecoy)", "")
        without_phospho = AASequence.fromString(seq)
        
        return without_phospho
    
    def createTheoreticalSpectra_(self, permutations, seq_without_phospho):
        """
        Create theoretical spectra for all phosphorylation site permutations.
        
        For each permutation of phosphorylation sites, this method:
        1. Creates a peptide sequence with phosphorylations at the specified positions
        2. Generates a theoretical spectrum with b- and y-ions
        3. Sets the spectrum name to the peptide sequence
        
        Parameters:
        -----------
        permutations : list
            List of all phosphorylation site permutations
        seq_without_phospho : AASequence
            Base peptide sequence without phosphorylations
            
        Returns:
        --------
        list
            List of theoretical spectra for all permutations
        """
        th_spectra = []
        for _ in range(len(permutations)):
            th_spectra.append(MSSpectrum())
        
        # Get the unmodified sequence as a string to check residue types
        seq_string = seq_without_phospho.toUnmodifiedString()
        
        for i in range(len(permutations)):
            seq = AASequence(seq_without_phospho)
            permu = 0
            
            for as_pos in range(seq.size()):
                if permu < len(permutations[i]) and as_pos == permutations[i][permu]:
                    # Apply the appropriate modification based on residue type
                    residue = seq_string[as_pos]
                    
                    # Determine which modification to apply based on residue type
                    if self.add_decoys_ and self.isPhosphoDecoySite(residue):
                        # This is a PhosphoDecoy site (A)
                        seq.setModification(as_pos, "PhosphoDecoy")
                        print(f"Debug: Set PhosphoDecoy modification at position {as_pos} (residue {residue})")
                    elif residue in ['S', 'T', 'Y']:
                        # This is a standard phosphorylation site (S, T, Y)
                        seq.setModification(as_pos, "Phospho")
                        print(f"Debug: Set Phospho modification at position {as_pos} (residue {residue})")
                    else:
                        # This should not happen
                        print(f"Warning: Attempted to set phosphorylation on non-phosphorylatable residue "
                             f"{residue} at position {as_pos}")
                    
                    permu += 1
                
                if permu == len(permutations[i]):
                    break
            
            # Generate b- and y-ions for charge 1
            self.spectrum_generator_.getSpectrum(th_spectra[i], seq, 1, 1)
            th_spectra[i].setName(seq.toString())
        
        return th_spectra
    
    def peakPickingPerWindowsInSpectrum_(self, real_spectrum):
        """
        Pick the top 10 intensity peaks for each 100 Da window in the spectrum.
        
        Parameters:
        -----------
        real_spectrum : MSSpectrum
            The experimental MS/MS spectrum
            
        Returns:
        --------
        list
            List of spectra containing the top 10 peaks for each window
        """
        windows_top10 = []
        
        spect_lower_bound = math.floor(real_spectrum[0].getMZ() / 100) * 100
        spect_upper_bound = math.ceil(real_spectrum[real_spectrum.size()-1].getMZ() / 100) * 100
        
        number_of_windows = int(math.ceil((spect_upper_bound - spect_lower_bound) / 100))
        for _ in range(number_of_windows):
            windows_top10.append(MSSpectrum())
        
        it_current_peak = 0
        window_upper_bound = spect_lower_bound + 100
        
        for current_window in range(number_of_windows):
            real_window = MSSpectrum()
            
            while it_current_peak < real_spectrum.size() and real_spectrum[it_current_peak].getMZ() <= window_upper_bound:
                real_window.push_back(real_spectrum[it_current_peak])
                it_current_peak += 1
            
            # Sort window by intensity (descending)
            real_window.sortByIntensity(True)
            
            # Take top 10 peaks (or fewer if window has less than 10 peaks)
            for i in range(min(10, real_window.size())):
                windows_top10[current_window].push_back(real_window[i])
            
            window_upper_bound += 100
        
        return windows_top10
    
    def calculatePermutationPeptideScores_(self, th_spectra, windows_top10):
        """
        Calculate scores for each permutation at 10 different peak depths.
        
        Parameters:
        -----------
        th_spectra : list
            List of theoretical spectra for all permutations
        windows_top10 : list
            List of experimental spectra windows with top 10 peaks
            
        Returns:
        --------
        list
            List of scores for each permutation at each peak depth
        """
        permutation_peptide_scores = []
        for _ in range(len(th_spectra)):
            permutation_peptide_scores.append([0.0] * 10)
        
        # For each phospho site assignment
        for idx, spectrum in enumerate(th_spectra):
            # The number of theoretical peaks (all b- and y-ions) correspond to the number of trials N
            N = spectrum.size()
            
            for i in range(1, 11):  # Peak depths 1-10
                n = 0  # Matched ions counter
                for window in windows_top10:  # Count matched ions over all 100 Da windows
                    n += self.numberOfMatchedIons_(spectrum, window, i)
                
                p = float(i) * self.base_match_probability_
                cumulative_score = self.computeCumulativeScore_(N, n, p)
                
                # abs is used to avoid -0 score values
                permutation_peptide_scores[idx][i-1] = abs((-10.0 * math.log10(cumulative_score)))
        
        return permutation_peptide_scores
    
    def rankWeightedPermutationPeptideScores_(self, peptide_site_scores):
        """
        Rank permutations by their weighted peptide scores.
        
        Parameters:
        -----------
        peptide_site_scores : list
            Scores for each permutation at each peak depth
            
        Returns:
        --------
        dict
            Dictionary mapping peptide scores to permutation indices
        """
        ranking = {}
        
        for i in range(len(peptide_site_scores)):
            weighted_score = self.peptideScore_(peptide_site_scores[i])
            ranking[weighted_score] = i
        
        return ranking
    
    def generateProFormaString_(self, peptide, ascores):
        """
        Generate a ProForma-like string with phosphorylation site localization scores.
        
        Parameters:
        -----------
        peptide : AASequence
            The peptide sequence
        ascores : dict
            Map of positions to AScore values
            
        Returns:
        --------
        str
            ProForma-like string with phosphorylation site scores
        """
        # Get the unmodified sequence as a string
        unmodified_str = peptide.toUnmodifiedString()
        
        # Create a map to store scores for each position
        position_scores = {}
        
        # Convert AScores to probabilities (0-1 range) for each phosphorylation site
        for position, ascore in ascores.items():
            # AScore is a -log10 scale, convert to probability
            # Use a simplified conversion that spreads values better
            probability = 1.0 - (1.0 / (1.0 + ascore))
            
            # Cap probability between 0 and 1
            probability = max(0.0, min(1.0, probability))
            
            # Store the probability for this position
            position_scores[position] = probability
        
        # Build the ProForma string
        result = ""
        for i in range(len(unmodified_str)):
            result += unmodified_str[i]
            
            # Check if this position has a phosphorylation
            if i in position_scores:
                # Determine the modification type based on the residue
                residue = unmodified_str[i]
                mod_name = ""
                
                if self.isPhosphoSite(residue):
                    # Standard phosphorylation site (S, T, Y)
                    mod_name = "Phospho"
                elif self.add_decoys_ and self.isPhosphoDecoySite(residue):
                    # PhosphoDecoy site (A)
                    mod_name = "PhosphoDecoy"
                else:
                    continue  # Skip if not a valid phosphorylation site
                
                # Format the score with 4 decimal places
                score_str = "{:.4f}".format(position_scores[i])
                # Add the ProForma-like annotation
                result += "[" + mod_name + "|score=" + score_str + "]"
        
        return result
