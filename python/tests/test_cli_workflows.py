"""CLI workflow tests for the Python interface."""

import sys
from pathlib import Path

import pytest
from click.testing import CliRunner


# Ensure package imports resolve relative to repository root
sys.path.insert(0, str(Path(__file__).parent.parent))

from gctb import cli  # noqa: E402


TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test" / "data"
TEST_BFILE = TEST_DATA_DIR / "uk10k_chr1_1mb"
TEST_PHENO = TEST_DATA_DIR / "test.phen"
TEST_INCLUDE_SNP = TEST_DATA_DIR / "causal.snplist"


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason="Test data not found")
def test_load_plink_data_with_filters():
    data = cli.load_plink_data(
        bfile=str(TEST_BFILE),
        pheno=str(TEST_PHENO),
        include_chr=1,
        include_snp_file=str(TEST_INCLUDE_SNP),
        verbose=False,
    )

    assert data.num_incd_snps > 0
    assert data.num_kept_inds > 0


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason="Test data not found")
def test_ldmatrix_make_command(tmp_path):
    runner = CliRunner()
    out_prefix = tmp_path / "test_ldm"

    result = runner.invoke(
        cli.main,
        [
            "ldmatrix",
            "make",
            "--bfile",
            str(TEST_BFILE),
            "--pheno",
            str(TEST_PHENO),
            "--out",
            str(out_prefix),
            "--type",
            "sparse",
            "--chisq-threshold",
            "5",
            "--ld-threshold",
            "0.0",
            "--window-width",
            "500000",
            "--verbose",
        ],
    )

    assert result.exit_code == 0, result.output

    bin_file = tmp_path / "test_ldm.ldm.sparse.bin"
    info_file = tmp_path / "test_ldm.ldm.sparse.info"

    assert bin_file.exists()
    assert info_file.exists()

