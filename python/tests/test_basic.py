"""
Basic tests for GCTB Python interface
"""

import pytest
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import gctb

def test_import():
    """Test that the module imports successfully"""
    assert gctb is not None
    assert hasattr(gctb, '__version__')

def test_snp_info_creation():
    """Test creating a SnpInfo object"""
    snp = gctb.SnpInfo(
        idx=0,
        id="rs12345",
        allele1="A",
        allele2="G",
        chr=1,
        gpos=0.5,
        ppos=1000000
    )
    assert snp.ID == "rs12345"
    assert snp.chrom == 1
    assert snp.physPos == 1000000
    assert snp.a1 == "A"
    assert snp.a2 == "G"

def test_ind_info_creation():
    """Test creating an IndInfo object"""
    ind = gctb.IndInfo(
        idx=0,
        fid="FAM001",
        pid="IND001",
        dad="0",
        mom="0",
        sex=1
    )
    assert ind.famID == "FAM001"
    assert ind.indID == "IND001"
    assert ind.sex == 1

def test_data_creation():
    """Test creating a Data object"""
    data = gctb.Data()
    assert data is not None
    assert data.num_snps == 0
    assert data.num_inds == 0

def test_timer():
    """Test Timer utility"""
    timer = gctb.Timer()
    timer.set_time()
    date = timer.get_date()
    assert isinstance(date, str)
    assert len(date) > 0

if __name__ == "__main__":
    pytest.main([__file__, "-v"])

