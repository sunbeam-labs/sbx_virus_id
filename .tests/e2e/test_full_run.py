import pytest
import shutil
import subprocess as sp
import sys
from pathlib import Path


@pytest.fixture
def setup(tmpdir):
    reads_fp = Path(".tests/data/reads/").resolve()
    hosts_fp = Path(".tests/data/hosts/").resolve()
    db_fp = Path(".tests/data/db/").resolve()

    project_dir = tmpdir / "project"

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = project_dir / "sunbeam_config.yml"

    config_str = f"sbx_virus_id: {{cenote_taker_db: {db_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield tmpdir, project_dir


@pytest.fixture
def run_sunbeam(setup):
    tmpdir, project_dir = setup

    output_fp = project_dir / "sunbeam_output"

    try:
        # Run the test job
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--conda-frontend",
                "conda",
                "--profile",
                project_dir,
                "all_virus_id",
                "--directory",
                tmpdir,
            ]
        )
    except sp.CalledProcessError as e:
        shutil.copytree(output_fp / "logs", "logs/")
        shutil.copytree(project_dir / "stats", "stats/")
        sys.exit(e)

    shutil.copytree(output_fp / "logs", "logs/")
    shutil.copytree(project_dir / "stats", "stats/")

    benchmarks_fp = project_dir / "stats"

    yield output_fp, benchmarks_fp


# def test_full_run(run_sunbeam):
#    output_fp, benchmarks_fp = run_sunbeam

#    all_align_fp = output_fp / "virus" / "summary" / "all_align_summary.txt"

# Check output
#    assert all_align_fp.exists()

#    with open(all_align_fp) as f:
#        header_line = f.readline()


@pytest.fixture
def dry_run_sunbeam(setup):
    tmpdir, project_dir = setup

    output_fp = project_dir / "sunbeam_output"

    try:
        # Run the test job
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--conda-frontend",
                "conda",
                "--profile",
                project_dir,
                "all_virus_id",
                "--directory",
                tmpdir,
                "-n",
            ]
        )
    except sp.CalledProcessError as e:
        shutil.copytree(output_fp / "logs", "logs/")
        shutil.copytree(project_dir / "stats", "stats/")
        sys.exit(e)

    shutil.copytree(output_fp / "logs", "logs/")
    shutil.copytree(project_dir / "stats", "stats/")

    benchmarks_fp = project_dir / "stats"

    yield output_fp, benchmarks_fp


def test_dry_run(dry_run_sunbeam):
    assert True
