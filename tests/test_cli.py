"""
Test CLI functionality.
"""

import pytest
import sys
import os
from unittest.mock import patch, MagicMock

# Add the parent directory to the path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def test_cli_help():
    """Test CLI help output."""
    from onsite.cli import create_parser
    
    parser = create_parser()
    
    # Test that parser has subcommands
    subparsers_action = None
    for action in parser._actions:
        if hasattr(action, 'choices') and action.choices:
            subparsers_action = action
            break
    
    assert subparsers_action is not None
    assert 'ascore' in subparsers_action.choices
    assert 'phosphors' in subparsers_action.choices
    assert 'lucxor' in subparsers_action.choices

def test_cli_ascore_parser():
    """Test AScore CLI parser."""
    from onsite.cli import create_parser
    
    parser = create_parser()
    # Test that ascore subcommand exists and can be parsed
    args = parser.parse_args(['ascore', '-in', 'test.mzML', '-id', 'test.idXML', '-out', 'result.idXML'])
    assert args.algorithm == 'ascore'
    assert args.in_file == 'test.mzML'

def test_cli_phosphors_parser():
    """Test PhosphoRS CLI parser."""
    from onsite.cli import create_parser
    
    parser = create_parser()
    # Test that phosphors subcommand exists and can be parsed
    args = parser.parse_args(['phosphors', '-in', 'test.mzML', '-id', 'test.idXML', '-out', 'result.idXML'])
    assert args.algorithm == 'phosphors'
    assert args.in_file == 'test.mzML'

def test_cli_lucxor_parser():
    """Test LucXor CLI parser."""
    from onsite.cli import create_parser
    
    parser = create_parser()
    # Test that lucxor subcommand exists and can be parsed
    args = parser.parse_args(['lucxor', '-in', 'test.mzML', '-id', 'test.idXML', '-out', 'result.idXML'])
    assert args.algorithm == 'lucxor'
    assert args.input_spectrum == 'test.mzML'

def test_cli_no_algorithm():
    """Test CLI with no algorithm specified."""
    from onsite.cli import create_parser
    
    parser = create_parser()
    args = parser.parse_args([])
    assert args.algorithm is None

def test_cli_unknown_algorithm():
    """Test CLI with unknown algorithm."""
    from onsite.cli import create_parser
    
    parser = create_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(['unknown_algorithm'])

@patch('onsite.cli.run_ascore')
def test_cli_ascore_execution(mock_run_ascore):
    """Test AScore execution through CLI."""
    from onsite.cli import main
    
    # Mock sys.argv
    with patch('sys.argv', ['onsite', 'ascore', '-in', 'test.mzML', '-id', 'test.idXML', '-out', 'result.idXML']):
        main()
        mock_run_ascore.assert_called_once()

@patch('onsite.cli.run_phosphors')
def test_cli_phosphors_execution(mock_run_phosphors):
    """Test PhosphoRS execution through CLI."""
    from onsite.cli import main
    
    # Mock sys.argv
    with patch('sys.argv', ['onsite', 'phosphors', '-in', 'test.mzML', '-id', 'test.idXML', '-out', 'result.idXML']):
        main()
        mock_run_phosphors.assert_called_once()

@patch('onsite.cli.run_lucxor')
def test_cli_lucxor_execution(mock_run_lucxor):
    """Test LucXor execution through CLI."""
    from onsite.cli import main
    
    # Mock sys.argv
    with patch('sys.argv', ['onsite', 'lucxor', '-in', 'test.mzML', '-id', 'test.idXML', '-out', 'result.idXML']):
        main()
        mock_run_lucxor.assert_called_once()
