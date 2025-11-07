"""
Test data I/O functionality with real test data
"""

import pytest
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import gctb

# Path to test data
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test" / "data"
TEST_BFILE = TEST_DATA_DIR / "uk10k_chr1_1mb"

@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason="Test data not found")
class TestDataIO:
    
    def test_read_fam_file(self):
        """Test reading .fam file"""
        data = gctb.Data()
        fam_file = str(TEST_BFILE) + ".fam"
        data.read_fam_file(fam_file)
        
        # Check that individuals were loaded
        assert data.num_inds > 0
        print(f"Loaded {data.num_inds} individuals")
        
        # Check individual info
        ind_vec = data.get_ind_info_vec()
        assert len(ind_vec) > 0
        
        # Check first individual
        first_ind = ind_vec[0]
        assert first_ind.famID is not None
        assert first_ind.indID is not None
    
    def test_read_bim_file(self):
        """Test reading .bim file"""
        data = gctb.Data()
        bim_file = str(TEST_BFILE) + ".bim"
        data.read_bim_file(bim_file)
        
        # Check that SNPs were loaded
        assert data.num_snps > 0
        print(f"Loaded {data.num_snps} SNPs")
        
        # Check SNP info
        snp_vec = data.get_snp_info_vec()
        assert len(snp_vec) > 0
        
        # Check first SNP
        first_snp = snp_vec[0]
        assert first_snp.ID is not None
        assert first_snp.chrom > 0
        print(f"First SNP: {first_snp.ID}, chr{first_snp.chrom}:{first_snp.physPos}")
    
    def test_read_plink_data(self):
        """Test reading complete PLINK dataset"""
        data = gctb.Data()
        
        # Read .fam
        fam_file = str(TEST_BFILE) + ".fam"
        data.read_fam_file(fam_file)
        
        # Read .bim
        bim_file = str(TEST_BFILE) + ".bim"
        data.read_bim_file(bim_file)
        
        # Read .bed (genotypes)
        bed_file = str(TEST_BFILE) + ".bed"
        data.read_bed_file(False, bed_file)
        
        print(f"Dataset: {data.num_snps} SNPs x {data.num_inds} individuals")
        assert data.num_snps > 0
        assert data.num_inds > 0

if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])

