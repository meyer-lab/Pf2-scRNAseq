"""
Test the cross validation accuracy.
"""

import numpy as np
import pytest

from ..imports import import_CCLE, import_citeseq, import_HTAN


@pytest.mark.parametrize(
    "import_func",
    [import_citeseq, import_HTAN, import_CCLE],
)
def test_imports(import_func):
    """Test import functions."""
    X = import_func()
    print(f"Data shape: {X.shape}")
    assert X.X.dtype == np.float32
