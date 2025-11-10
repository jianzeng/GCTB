"""CLI workflow tests for the Python interface."""

import sys
from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner


# Ensure package imports resolve relative to repository root
sys.path.insert(0, str(Path(__file__).parent.parent))

from gctb import cli  # noqa: E402
from gctb import workflows  # noqa: E402


TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test" / "data"
TEST_BFILE = TEST_DATA_DIR / "uk10k_chr1_1mb"
TEST_PHENO = TEST_DATA_DIR / "test.phen"
TEST_INCLUDE_SNP = TEST_DATA_DIR / "causal.snplist"
TEST_LDM = TEST_DATA_DIR / "test_ldm.ldm.sparse"
TEST_GWAS = TEST_DATA_DIR / "test_gwas_summary.ma"
TEST_SNP_RES = Path(__file__).parent.parent.parent / "test" / "res" / "test.snpRes"


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


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason="Test data not found")
def test_window_credible_sets():
    data = cli.load_summary_data(
        ldm_prefix=str(TEST_LDM),
        gwas_file=str(TEST_GWAS),
        verbose=False,
    )

    result = workflows.compute_window_credible_sets(
        data,
        snp_results=str(TEST_SNP_RES),
        window_width=500_000,
        pip_threshold=0.8,
    )

    assert "windows" in result and result["windows"]
    first = result["windows"][0]
    assert "credible_set" in first
    assert isinstance(first["credible_set"], list)


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason="Test data not found")
def test_window_pip_with_overlap():
    data = cli.load_summary_data(
        ldm_prefix=str(TEST_LDM),
        gwas_file=str(TEST_GWAS),
        verbose=False,
    )

    windows = workflows.compute_window_pip(
        data,
        snp_results=str(TEST_SNP_RES),
        window_width=500_000,
        step_size=250_000,
    )

    assert windows, "Expected at least one window summary"
    assert "snps" in windows[0] and windows[0]["snps"]
    assert 0.0 <= windows[0]["window_pip"] <= 1.0


def test_pip_to_pvalues_roundtrip():
    pips = [0.3, 0.6, 0.1]
    pvals = workflows.pip_to_pvalues(pips, prop_null=0.9)
    assert len(pvals) == len(pips)
    assert all(p >= 0.0 for p in pvals)
    # Higher PIP should correspond to smaller p-value proxy
    sorted_pairs = sorted(zip(pips, pvals), key=lambda t: t[0], reverse=True)
    assert sorted_pairs[0][1] <= sorted_pairs[-1][1]


def test_mcmc_sparse_extraction_helpers():
    class _StubSamples:
        storage_mode = "sparse"

        def sparse_data(self):
            import numpy as np  # local import for stub
            rows = np.array([0, 1, 1], dtype=np.int32)
            cols = np.array([0, 0, 2], dtype=np.int32)
            data = np.array([1.0, 0.5, -0.25], dtype=np.float32)
            shape = (2, 3)
            return rows, cols, data, shape

    stub = _StubSamples()
    rows, cols, vals, shape = workflows.mcmc_samples_sparse_matrix(stub)

    assert rows.tolist() == [0, 1, 1]
    assert cols.tolist() == [0, 0, 2]
    assert np.allclose(vals, [1.0, 0.5, -0.25])
    assert shape == (2, 3)

    if workflows.sp is not None:
        csr = workflows.mcmc_samples_to_csr(stub)
        assert csr.shape == (2, 3)
        assert csr.nnz == 3


@pytest.mark.skipif(not TEST_DATA_DIR.exists(), reason="Test data not found")
@pytest.mark.skip(reason="Known limitation: pybind11 MCMC bindings abort with Trace/BPT trap when run in-process on macOS (see KNOWN_ISSUES).")
def test_multi_chain_sbayes(tmp_path):
    data = cli.load_summary_data(
        ldm_prefix=str(TEST_LDM),
        gwas_file=str(TEST_GWAS),
        verbose=False,
    )

    chains = workflows.run_multi_chain_sbayes(
        data,
        sbayes_type="R",
        heritability=0.5,
        pi=0.01,
        num_chains=2,
        chain_length=20,
        burnin=5,
        thin=5,
        output_prefix=str(tmp_path / "multi_chain"),
        random_start=True,
        verbose=False,
    )

    assert len(chains) == 2
    assert all(isinstance(chain, list) and len(chain) > 0 for chain in chains)

    gelman = tmp_path / "multi_chain.gelman"
    assert gelman.exists()

