"""
Test configuration and fixtures for OnSite tests.
"""

import pytest
import sys
import os
from pathlib import Path

# Add the parent directory to the path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@pytest.fixture(scope="session")
def data_dir():
    """Get the data directory path."""
    return Path(__file__).parent.parent / "data"


@pytest.fixture(scope="session")
def idxml_file(data_dir):
    """Get the idXML file path."""
    return data_dir / "1_consensus_fdr_filter_pep.idXML"


@pytest.fixture(scope="session")
def mzml_file(data_dir):
    """Get the mzML file path."""
    return data_dir / "1.mzML"


@pytest.fixture(scope="session")
def data_files_exist(data_dir, idxml_file, mzml_file):
    """Check that required data files exist."""
    if not data_dir.exists():
        pytest.skip("Data directory does not exist")
    if not idxml_file.exists():
        pytest.skip("idXML file does not exist")
    if not mzml_file.exists():
        pytest.skip("mzML file does not exist")
    return True


@pytest.fixture
def temp_output_dir():
    """Create a temporary directory for output files."""
    import tempfile
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir


@pytest.fixture
def sample_output_file(temp_output_dir):
    """Get a sample output file path."""
    return os.path.join(temp_output_dir, "test_output.idXML")


# Markers for different test categories
def pytest_configure(config):
    """Configure pytest markers."""
    config.addinivalue_line("markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')")
    config.addinivalue_line("markers", "integration: marks tests as integration tests")
    config.addinivalue_line("markers", "data: marks tests that require data files")
    config.addinivalue_line("markers", "cli: marks tests that test CLI functionality")
    config.addinivalue_line("markers", "algorithm: marks tests that test algorithm functionality")


# Skip tests if required data files are missing
def pytest_collection_modifyitems(config, items):
    """Modify test collection to skip tests if data files are missing."""
    data_dir = Path(__file__).parent.parent / "data"
    idxml_file = data_dir / "1_consensus_fdr_filter_pep.idXML"
    mzml_file = data_dir / "1.mzML"
    
    if not data_dir.exists() or not idxml_file.exists() or not mzml_file.exists():
        for item in items:
            if "data" in item.keywords:
                item.add_marker(pytest.mark.skip(reason="Required data files not found"))
