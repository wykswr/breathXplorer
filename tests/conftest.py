import pytest
from pathlib import Path


@pytest.fixture
def sample():
    return Path(__file__).parent / 'sample.mzML'
