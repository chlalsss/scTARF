"""Exercise the actual LSI CLI with an ordered sparse count fixture."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from scipy import sparse
from scipy.io import mmwrite


class LsiCliTest(unittest.TestCase):
    def test_shape_and_existing_output_guard(self):
        script = Path('scripts/03_nonpaired_analysis/02_build_granja_lsi.py').resolve()
        with tempfile.TemporaryDirectory(prefix='lsi space ') as tmp:
            out = Path(tmp) / 'multiomics_module1_pairing/method_B_scSHEFT/scsheft_input'
            out.mkdir(parents=True)
            counts = sparse.csr_matrix(np.random.default_rng(7).poisson(0.5, size=(55, 60)))
            mmwrite(out / 'atac_counts.mtx', counts)
            (out / 'atac_cells.csv').write_text('cell\n' + '\n'.join(f'cell_{i}' for i in range(60)) + '\n')
            env = dict(os.environ, SCTARF_BLOOD_ROOT=tmp, PYTHONIOENCODING='utf-8')
            subprocess.run([sys.executable, str(script)], env=env, check=True, capture_output=True)
            embedding = np.load(out / 'blood_lsi.npy')
            self.assertEqual(embedding.shape, (60, 50))
            self.assertTrue(np.isfinite(embedding).all())
            np.testing.assert_allclose(np.linalg.norm(embedding, axis=1), 1, atol=1e-6)
            saved = (out / 'blood_lsi.npy').read_bytes()
            rerun = subprocess.run([sys.executable, str(script)], env=env, capture_output=True)
            self.assertNotEqual(rerun.returncode, 0)
            self.assertEqual(saved, (out / 'blood_lsi.npy').read_bytes())


if __name__ == '__main__':
    unittest.main()
