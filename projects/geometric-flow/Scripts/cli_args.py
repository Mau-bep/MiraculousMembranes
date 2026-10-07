"""Shared argparse helpers for the Create_subjob_*.py scripts."""
import argparse


def number_str(text):
    """Accept a number but keep the text as typed.

    The subjob and config file names are built from the argument text, so "5"
    must stay "5" (not 5.0): the Mem_*.sh scripts sbatch those names.
    """
    try:
        float(text)
    except ValueError:
        raise argparse.ArgumentTypeError("{!r} is not a number".format(text))
    return text
