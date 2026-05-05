"""
Basic unit tests for the CpG Methylation ML pipeline.
Run with: pytest tests/ -v
"""

import pytest
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier


def make_mock_methylation(n_samples=20, n_cpgs=100, random_state=42):
    """Generate a small mock methylation matrix for testing."""
    rng = np.random.default_rng(random_state)
    X = pd.DataFrame(
        rng.uniform(0, 1, size=(n_samples, n_cpgs)),
        columns=[f"CpG_{i}" for i in range(n_cpgs)]
    )
    y = pd.Series([0] * (n_samples // 2) + [1] * (n_samples // 2))
    return X, y


def test_mock_data_shape():
    X, y = make_mock_methylation()
    assert X.shape == (20, 100)
    assert len(y) == 20


def test_beta_value_range():
    X, _ = make_mock_methylation()
    assert (X >= 0).all().all(), "Beta values should be >= 0"
    assert (X <= 1).all().all(), "Beta values should be <= 1"


def test_class_balance():
    _, y = make_mock_methylation()
    assert y.sum() == len(y) // 2, "Should have equal classes"


def test_random_forest_fits():
    X, y = make_mock_methylation()
    rf = RandomForestClassifier(n_estimators=10, random_state=42)
    rf.fit(X, y)
    preds = rf.predict(X)
    assert len(preds) == len(y)


def test_feature_importance_length():
    X, y = make_mock_methylation()
    rf = RandomForestClassifier(n_estimators=10, random_state=42)
    rf.fit(X, y)
    assert len(rf.feature_importances_) == X.shape[1]


def test_no_nan_in_mock_data():
    X, _ = make_mock_methylation()
    assert not X.isnull().any().any(), "Mock data should have no NaN values"
