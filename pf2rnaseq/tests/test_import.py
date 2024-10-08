"""
Test the cross validation accuracy.
"""

import pytest
import numpy as np
from ..imports import (
    import_citeseq,
    import_HTAN,
    import_CCLE
)


@pytest.mark.parametrize(
    "import_func",
    [import_citeseq, import_HTAN, import_CCLE],
)
def test_imports(import_func):
    """Test import functions."""
    X = import_func()
    print(f"Data shape: {X.shape}")
    assert X.X.dtype == np.float32
