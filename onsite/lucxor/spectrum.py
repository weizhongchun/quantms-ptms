"""
Spectrum module.

This module contains the Spectrum class for handling mass spectrometry data.
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Optional
import pyopenms

logger = logging.getLogger(__name__)

class Spectrum:
    """Class representing a mass spectrum."""
    
    def __init__(self, mz_array=None, intensity_array=None):
        """
        Initialize a new Spectrum instance.
        
        Args:
            mz_array: Array of m/z values
            intensity_array: Array of intensity values
        """
        self.N = 0  # number of peaks
        self.max_i_index = 0  # index of max intensity peak
        self.max_i = 0.0  # max intensity value
        
        # Peak arrays
        self.mz = []  # m/z values
        self.raw_intensity = []  # raw intensity values
        self.rel_intensity = []  # relative intensity values (0-100)
        self.norm_intensity = []  # normalized intensity values (log scale)
        
        if mz_array is not None and intensity_array is not None:
            self._init_peaks(mz_array, intensity_array)
    
    def _init_peaks(self, mz_array, intensity_array):
        """Initialize peak arrays from input data."""
        # Ensure arrays have same length
        N = min(len(mz_array), len(intensity_array))
        
        # Filter out zero intensity peaks
        cand_peaks = []
        for i in range(N):
            if intensity_array[i] > 0:
                cand_peaks.append(i)
        
        self.N = len(cand_peaks)
        
        # Initialize arrays
        self.mz = np.zeros(self.N)
        self.raw_intensity = np.zeros(self.N)
        self.rel_intensity = np.zeros(self.N)
        self.norm_intensity = np.zeros(self.N)
        
        # Record peak data
        for i, idx in enumerate(cand_peaks):
            intensity = intensity_array[idx]
            if intensity > self.max_i:
                self.max_i = intensity
                self.max_i_index = idx
            
            self.mz[i] = mz_array[idx]
            self.raw_intensity[i] = intensity
        
        # Calculate relative intensities
        self.calc_relative_intensity()
    
    def calc_relative_intensity(self):
        """Calculate relative intensities (0-100) based on max intensity."""
        if self.max_i > 0:
            self.rel_intensity = (self.raw_intensity / self.max_i) * 100.0
    
    def median_normalize_spectra(self):
        """
        Normalize intensities using median.
        Computes log(rel_intensity/median_intensity).
        """
        if self.N == 0:
            return
            
        # Sort relative intensities
        sorted_intensities = sorted(self.rel_intensity)
        
        # Calculate median
        mid = self.N // 2
        if self.N % 2 == 0:  # even number of peaks
            a = sorted_intensities[mid - 1]
            b = sorted_intensities[mid]
            median_i = (a + b) / 2.0
        else:  # odd number of peaks
            median_i = sorted_intensities[mid]
        
        # Calculate normalized intensities
        for i in range(self.N):
            d = self.rel_intensity[i] / median_i
            if d > 0:  # avoid log(0)
                self.norm_intensity[i] = np.log(d)
            else:
                self.norm_intensity[i] = float('-inf')
    
    def get_peaks(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get peak arrays.
        
        Returns:
            Tuple of (m/z array, intensity array)
        """
        return self.mz, self.raw_intensity
    
    def is_empty(self) -> bool:
        """
        Check if spectrum is empty.
        
        Returns:
            True if spectrum has no peaks
        """
        return self.N == 0
    
    def find_index_by_mz(self, mz: float) -> int:
        """
        Find peak index by m/z value.
        
        Args:
            mz: m/z value to find
            
        Returns:
            Index of peak with matching m/z, or -1 if not found
        """
        for i in range(self.N):
            if abs(self.mz[i] - mz) < 1e-6:  # floating point comparison
                return i
        return -1 