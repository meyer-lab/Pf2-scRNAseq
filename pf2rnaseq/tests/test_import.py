"""
Test the cross validation accuracy.
"""

import numpy as np
import pytest
from scipy.sparse import csr_array

from ..imports import import_cytokine, import_Heiser, import_MouseImmune


@pytest.mark.parametrize(
    "import_func",
    [import_cytokine, import_MouseImmune],
)
def test_imports(import_func):
    """Test import functions."""
    X = import_func()
    print(f"Data shape: {X.shape}")
    assert X.X.dtype == np.float32


@pytest.mark.parametrize(
    "deviance",
    [False, True],
)
def test_import_Heiser(deviance):
    """Test import functions."""
    X = import_Heiser(deviance)
    print(f"Data shape: {X.shape}")

    Xdata: csr_array | np.ndarray = X.X # type: ignore

    if isinstance(X.X, csr_array):
        assert np.all(np.isfinite(Xdata.data))
    else:
        assert np.all(np.isfinite(Xdata)) # type: ignore