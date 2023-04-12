from pathlib import Path

import pytest


@pytest.fixture
def sample():
    return Path(__file__).parent / 'sample.mzML'
