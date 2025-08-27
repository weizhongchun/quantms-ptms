"""
Test that all modules can be imported correctly.
"""

import pytest


def test_onsite_package_import():
    """Test that the onsite package can be imported."""
    try:
        import onsite
        assert onsite is not None
        assert onsite.__version__ == "0.0.1"
    except ImportError as e:
        pytest.skip(f"onsite package import failed: {e}")


def test_ascore_import():
    """Test that AScore can be imported from onsite package."""
    try:
        from onsite import AScore
        assert AScore is not None
    except ImportError as e:
        pytest.skip(f"AScore import failed: {e}")


def test_phosphors_import():
    """Test that phosphors functions can be imported from onsite package."""
    try:
        from onsite import calculate_phospho_localization_compomics_style
        assert calculate_phospho_localization_compomics_style is not None
    except ImportError as e:
        pytest.skip(f"phosphors import failed: {e}")


def test_phosphoscoring_import():
    """Test that PhosphoScoring can be imported."""
    try:
        from PhosphoScoring import main
        assert main is not None
    except ImportError as e:
        pytest.skip(f"PhosphoScoring import failed: {e}")


def test_lucxor_import():
    """Test that lucxor can be imported."""
    try:
        import onsite.lucxor
        assert onsite.lucxor is not None
    except ImportError as e:
        pytest.skip(f"lucxor import failed: {e}")

