"""
Filip Mazurek - 9/1/2019

Utility for suppressing output when calling functions
Based on answer in https://stackoverflow.com/questions/8391411/suppress-calls-to-print-python

Usage:

with hidden_prints():
    call_function()
"""

import os
import sys


class hidden_prints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
