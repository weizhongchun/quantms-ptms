#!/usr/bin/env python3
"""
OnSite CLI - Integrated command line interface for phosphorylation site localization tools.

This module provides a unified interface for accessing AScore, PhosphoRS, and LucXor algorithms.
"""

import click
import sys
from onsite.lucxor.cli import lucxor
from onsite.phosphors.cli import phosphors
from onsite.ascore.cli import ascore


@click.group()
@click.version_option(version="0.0.1")
def cli():
    """
    OnSite: Mass spectrometry post-translational modification localization tool

    Available algorithms:
      ascore      AScore algorithm for phosphorylation site localization
      phosphors   PhosphoRS algorithm for phosphorylation site localization
      lucxor      LucXor (LuciPHOr2) algorithm for PTM localization

    Examples:
      onsite ascore -in spectra.mzML -id identifications.idXML -out results.idXML
      onsite phosphors -in spectra.mzML -id identifications.idXML -out results.idXML
      onsite lucxor -in spectra.mzML -id identifications.idXML -out results.idXML
    """
    pass


# Add the individual CLI commands to the main CLI group
cli.add_command(ascore)
cli.add_command(phosphors)
cli.add_command(lucxor)


def main():
    """Main entry point for OnSite CLI."""
    try:
        cli()
    except KeyboardInterrupt:
        click.echo("\nOperation cancelled by user")
        sys.exit(1)
    except Exception as e:
        click.echo(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
