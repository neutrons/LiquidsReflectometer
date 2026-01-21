import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from lr_reduction.output import RunCollection, read_file


class TestRunCollection:
    """Test cases for RunCollection class"""

    def test_add_with_dq(self):
        """Test adding data with explicit dq values"""
        rc = RunCollection()
        q = np.array([0.1, 0.2, 0.3])
        r = np.array([1.0, 0.5, 0.2])
        dr = np.array([0.1, 0.05, 0.02])
        dq = np.array([0.01, 0.02, 0.03])
        meta = {"experiment": "test", "run_number": 1}

        rc.add(q, r, dr, meta, dq=dq)

        assert len(rc.collection) == 1
        assert np.array_equal(rc.collection[0]["q"], q)
        assert np.array_equal(rc.collection[0]["dq"], dq)

    def test_add_without_dq(self):
        """Test adding data without dq, should compute from metadata"""
        rc = RunCollection()
        q = np.array([0.1, 0.2, 0.3])
        r = np.array([1.0, 0.5, 0.2])
        dr = np.array([0.1, 0.05, 0.02])
        meta = {"experiment": "test", "run_number": 1, "dq_over_q": 0.05}

        rc.add(q, r, dr, meta)

        expected_dq = 0.05 * q
        assert np.allclose(rc.collection[0]["dq"], expected_dq)

    def test_merge_single_run(self):
        """Test merging with a single run"""
        rc = RunCollection()
        q = np.array([0.1, 0.2, 0.3])
        r = np.array([1.0, 0.5, 0.2])
        dr = np.array([0.1, 0.05, 0.02])
        dq = np.array([0.01, 0.02, 0.03])
        meta = {"experiment": "test", "run_number": 1}

        rc.add(q, r, dr, meta, dq=dq)
        rc.merge()

        assert np.array_equal(rc.qz_all, q)
        assert np.array_equal(rc.refl_all, r)
        assert np.array_equal(rc.d_refl_all, dr)
        assert np.array_equal(rc.d_qz_all, dq)

    def test_merge_multiple_runs_sorted(self):
        """Test merging multiple runs with sorting"""
        rc = RunCollection()
        q1 = np.array([0.3, 0.1])
        r1 = np.array([0.2, 1.0])
        dr1 = np.array([0.02, 0.1])
        dq1 = np.array([0.03, 0.01])
        meta1 = {"experiment": "test", "run_number": 1}

        q2 = np.array([0.2])
        r2 = np.array([0.5])
        dr2 = np.array([0.05])
        dq2 = np.array([0.02])
        meta2 = {"experiment": "test", "run_number": 2}

        rc.add(q1, r1, dr1, meta1, dq=dq1)
        rc.add(q2, r2, dr2, meta2, dq=dq2)
        rc.merge()

        # Check sorting by q
        expected_q = np.array([0.1, 0.2, 0.3])
        assert np.array_equal(rc.qz_all, expected_q)
        assert np.array_equal(rc.refl_all, np.array([1.0, 0.5, 0.2]))

    def test_merge_with_average_overlap(self):
        """Test merging with averaging overlapping points"""
        rc = RunCollection(average_overlap=True)

        # Add two runs with overlapping q values
        q1 = np.array([0.1, 0.2])
        r1 = np.array([1.0, 2.0])
        dr1 = np.array([0.1, 0.2])
        dq1 = np.array([0.01, 0.01])
        meta1 = {"experiment": "test", "run_number": 1}

        q2 = np.array([0.2])
        r2 = np.array([0.5])
        dr2 = np.array([0.05])
        dq2 = np.array([0.02])
        meta2 = {"experiment": "test", "run_number": 2}

        rc.add(q1, r1, dr1, meta1, dq=dq1)
        rc.add(q2, r2, dr2, meta2, dq=dq2)
        rc.merge()

        assert np.array_equal(rc.qz_all, np.array([0.1, 0.2]))
        assert np.array_equal(rc.refl_all, np.array([1.0, (2.0 + 0.5) / 2]))

    def test_save_ascii(self):
        """Test saving data to ASCII file"""
        rc = RunCollection()
        q = np.array([0.1, 0.2])
        r = np.array([1.0, 0.5])
        dr = np.array([0.1, 0.05])
        dq = np.array([0.01, 0.02])
        meta = {
            "experiment": "test_exp",
            "run_number": 1,
            "run_title": "Test Run",
            "start_time": "2021-01-01 10:00:00",
            "time": "2021-01-01 10:05:00",
            "theta": 0.01,
            "norm_run": 0,
            "wl_min": 1.0,
            "wl_max": 10.0,
            "q_min": 0.01,
            "q_max": 1.0,
        }

        rc.add(q, r, dr, meta, dq=dq)

        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = Path(tmpdir) / "output.txt"
            rc.save_ascii(str(file_path))

            with open(file_path, "r") as f:
                content = f.read()
                assert "# Experiment test_exp Run 1" in content
                assert "# Run title: Test Run" in content
                assert "Q [1/Angstrom]" in content
                assert "0.1" in content

    def test_add_from_file(self):
        """Test reading and adding data from file"""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = Path(tmpdir) / "input.txt"

            # Create a test file
            meta = {"test": "data"}
            with open(file_path, "w") as f:
                f.write(f"# Meta:{json.dumps(meta)}\n")
                f.write("0.1  1.0  0.1  0.01\n")
                f.write("0.2  0.5  0.05  0.02\n")

            rc = RunCollection()
            rc.add_from_file(str(file_path))

            assert len(rc.collection) == 1
            assert rc.collection[0]["info"]["test"] == "data"

    def test_read_file_valid(self):
        """Test reading a valid data file"""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = Path(tmpdir) / "data.txt"
            meta = {"run": 1, "exp": "test"}

            with open(file_path, "w") as f:
                f.write(f"# Meta:{json.dumps(meta)}\n")
                f.write("0.1  1.0  0.1  0.01\n")
                f.write("0.2  0.5  0.05  0.02\n")

            q, r, dr, dq, read_meta = read_file(str(file_path))

            assert len(q) == 2
            assert q[0] == pytest.approx(0.1)
            assert read_meta == meta

    def test_read_file_empty(self):
        """Test reading an empty or invalid data file"""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = Path(tmpdir) / "empty.txt"

            with open(file_path, "w") as f:
                f.write("# No data\n")

            q, r, dr, dq, meta = read_file(str(file_path))

            assert len(q) == 0
            assert meta == {}
