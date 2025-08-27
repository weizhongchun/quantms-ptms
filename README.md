# OnSite

Mass spectrometry post-translational modification localization tool.

This package provides tools for phosphorylation site localization and scoring using various algorithms including AScore and PhosphoRS.

## Features

- **AScore Algorithm**: Implementation of the AScore algorithm for phosphorylation site localization
- **PhosphoRS Algorithm**: Implementation of the PhosphoRS algorithm for phosphorylation site localization  
- **PhosphoScoring Pipeline**: Complete workflow tool for processing MS/MS data files
- **LucXor**: Python port of LuciPHOr2 for accurate PTM localization with probabilistic scoring

## Installation

### Using Poetry (Recommended)

```bash
# Clone the repository
git clone https://github.com/bigbio/onsite.git
cd onsite

# Install with Poetry
poetry install

# Activate the virtual environment
poetry shell
```

### Using pip

```bash
pip install onsite
```

## Usage

### Command Line Interface

#### PhosphoScoring Pipeline

Process MS/MS data files using the PhosphoScoring pipeline:

```bash
# Using Poetry
poetry run phospho-scoring -in spectra.mzML -id identifications.idXML -out results.idXML

# Or after installation
phospho-scoring -in spectra.mzML -id identifications.idXML -out results.idXML
```

#### LucXor Tool

Use the LucXor tool for PTM localization:

```bash
# Using Poetry
poetry run lucxor -in spectra.mzML -id identifications.idXML -out results.idXML

# Or after installation
lucxor -in spectra.mzML -id identifications.idXML -out results.idXML

# With custom parameters
python -m onsite.lucxor.cli \
    -in spectra.mzML \
    -id identifications.idXML \
    -out results.idXML \
    --fragment_tolerance 0.05 \
    --fragment_tolerance_unit Da \
    --min_peptide_length 6 \
    --max_peptide_length 50
```

### Python API

#### AScore Algorithm

```python
from pyopenms import *
from onsite import AScore

# Initialize AScore
ascore = AScore()

# Set parameters
ascore.setParameter("fragment_mass_tolerance", 0.05)
ascore.setParameter("fragment_mass_unit", "Da")

# Process a peptide hit
result = ascore.compute(peptide_hit, spectrum)
```

#### PhosphoRS Algorithm

```python
from onsite import calculate_phospho_localization_compomics_style

# Calculate phosphorylation site probabilities
site_probs, isomer_details = calculate_phospho_localization_compomics_style(
    peptide_hit, 
    spectrum,
    fragment_tolerance=0.05
)
```

#### LucXor API

```python
from onsite.lucxor import PyLuciPHOr2, LucXor, Peptide, PSM
from onsite.lucxor.models import CIDModel, HCDModel
from onsite.lucxor.spectrum import Spectrum
from onsite.lucxor.flr import FLRCalculator

# Initialize the processor
processor = PyLuciPHOr2()

# Configure parameters
config = LucXorConfig()
config.fragment_tolerance = 0.05
config.fragment_tolerance_unit = "Da"

# Process data
results = processor.process(spectra, identifications)
```

#### Advanced LucXor Usage

```python
from onsite.lucxor.core import CoreProcessor
from onsite.lucxor.config import LucXorConfig

# Custom configuration
config = LucXorConfig()
config.fragment_tolerance = 0.05
config.fragment_tolerance_unit = "Da"
config.min_peptide_length = 6
config.max_peptide_length = 50
config.threads = 4

# Initialize processor
processor = CoreProcessor(config)

# Process data
results = processor.process_psms(psm_list, spectrum_map)
```

## LucXor Algorithm

LucXor implements the complete LuciPHOr2 algorithm:

1. **Spectrum Processing**: Peak picking and noise filtering
2. **Theoretical Spectrum Generation**: Generate all possible modification site combinations
3. **Fragment Matching**: Match experimental and theoretical fragments
4. **Probabilistic Scoring**: Calculate localization probabilities using binomial statistics
5. **False Localization Rate**: Estimate confidence in localization results
6. **Result Ranking**: Rank results by probability and confidence

## LucXor Configuration

Key parameters:

- `fragment_tolerance`: Mass tolerance for fragment matching (default: 0.05 Da)
- `fragment_tolerance_unit`: Unit for tolerance (Da or ppm)
- `min_peptide_length`: Minimum peptide length to process
- `max_peptide_length`: Maximum peptide length to process
- `threads`: Number of parallel processing threads

## LucXor Performance

- **Processing Speed**: ~100-1000 PSMs/second depending on data complexity
- **Memory Usage**: ~1-2 GB for typical datasets
- **Scalability**: Linear scaling with dataset size and number of threads

## LucXor Troubleshooting

### Common Issues

1. **Memory Errors**: Reduce `max_peptide_length` or use fewer threads
2. **Slow Processing**: Increase number of threads or reduce fragment tolerance
3. **Poor Results**: Check fragment tolerance settings and data quality

## Dependencies

- Python 3.11+
- pyOpenMS 3.4.0+
- NumPy 2.3.2+
- SciPy 1.16.1+

## Development

### Setup Development Environment

```bash
# Install development dependencies
poetry install --with dev

# Run tests
poetry run pytest

# Format code
poetry run black .

# Lint code
poetry run flake8
```

### Project Structure

```
onsite/
├── onsite/                 # Main package directory
│   ├── __init__.py        # Package initialization
│   ├── ascore.py          # AScore algorithm implementation
│   ├── phosphors.py       # PhosphoRS algorithm implementation
│   └── lucxor/            # LucXor utilities
│       ├── cli.py         # Command-line interface
│       ├── core.py        # Core processing logic
│       ├── models.py      # Fragmentation models
│       └── ...            # Other modules
├── PhosphoScoring.py      # Processing pipeline tool
├── pyproject.toml         # Poetry configuration
├── README.md              # This file
└── LICENSE                # MIT License
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## Citation

If you use this software in your research, please cite:

```
BigBio Stack. (2025). OnSite: Mass spectrometry post-translational 
localization tool. https://github.com/bigbio/onsite
```

If you use LucXor specifically, please also cite the original LuciPHOr2 paper:

```
Fermin, D., et al. (2013). LuciPHOr2: site identification of generic 
post-translational modifications from tandem mass spectrometry data. 
Nature Methods, 10(8), 744-746.
```
