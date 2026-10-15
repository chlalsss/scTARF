"""Import wrappers without running scDRS and validate PATH command resolution."""
import os
from pathlib import Path
import runpy
import shutil
import unittest


class CommandPaths(unittest.TestCase):
    def test_python_command_is_resolved_from_path(self):
        prior = os.environ.get('SCTARF_PYTHON')
        os.environ['SCTARF_PYTHON'] = 'python'
        try:
            self.assertIsNotNone(shutil.which('python'))
            wrapper = runpy.run_path('scripts/03_nonpaired_analysis/22_run_scdrs_cite.py', run_name='path_test')
            self.assertTrue(Path(wrapper['PYTHON']).is_file(), 'Bare command must resolve through PATH')
        finally:
            if prior is None:
                os.environ.pop('SCTARF_PYTHON', None)
            else:
                os.environ['SCTARF_PYTHON'] = prior


if __name__ == '__main__':
    unittest.main()
