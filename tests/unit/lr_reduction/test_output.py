import json
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pytest

from lr_reduction.output import RunCollection, read_file

# What read_file prints when it finds no points, at the base and after (U4).
_NO_POINTS = "Could not read file. It may have no points\n"


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

    def test_read_file_empty(self, capsys):
        """Test reading an empty or invalid data file: four empties, no meta, the message, and no warning (U4)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            file_path = Path(tmpdir) / "empty.txt"

            with open(file_path, "w") as f:
                f.write("# No data\n")

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                q, r, dr, dq, meta = read_file(str(file_path))

            assert [str(w.message) for w in caught] == []
            assert (q, r, dr, dq) == ([], [], [], [])
            assert meta == {}
            assert capsys.readouterr().out == _NO_POINTS

    @pytest.mark.parametrize(
        "text, meta",
        [
            ('# Meta:{"run": 1}\n', {"run": 1}),
            ('# Meta:{"run": 1}\n\n   \n  # a comment after blanks\n', {"run": 1}),
        ],
        ids=["meta-only", "meta-blanks-and-comments"],
    )
    def test_read_file_without_a_data_row(self, tmp_path, capsys, text, meta):
        """U4: a file with no data row gives four empties and its meta, with the base's message and no warning"""
        file_path = tmp_path / "empty.txt"
        file_path.write_text(text)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = read_file(str(file_path))

        assert [str(w.message) for w in caught] == []
        assert result == ([], [], [], [], meta)
        assert capsys.readouterr().out == _NO_POINTS

    @pytest.mark.parametrize(
        "text, meta",
        [
            ('# Meta:{"run": 2}\n0.1 1.0 0.1\n0.2 0.5 0.05\n', {"run": 2}),
            ('# Meta:{"run": 2}\n0.1 1.0 0.1 0.01 7\n', {"run": 2}),
            ('# Meta:{"run": 2}\n0.1 1.0 0.1\n', {"run": 2}),
            ('# Meta:{"run": 2}\n0.1\n0.2\n0.3\n', {"run": 2}),
            ('# Meta:{"run": 2}\n0.1 1.0 0.1 0.01\n0.2 0.5\n', {"run": 2}),
            ('# Meta:{"run": 2}\nq r dr dq\n', {"run": 2}),
            ("5\n", {}),
            ('# Meta:{"start_time": "x"}\n0.01\n', {"start_time": "x"}),
        ],
        ids=[
            "rows-of-three", "one-row-of-five", "one-row-of-three", "three-rows-of-one", "ragged", "not-numbers",
            "a-single-number", "a-single-number-after-meta",
        ],
    )
    def test_read_file_whose_rows_are_not_four_numbers(self, tmp_path, capsys, text, meta):
        """U4: every shape the base read as "no points" gives four empty lists and the meta, with the base's
        message. That includes the 0-d shape: np.loadtxt reads a single number as a 0-d array, whose unpack
        raises TypeError, not ValueError (v2, B-1: v1 let it out)."""
        file_path = tmp_path / "bad.txt"
        file_path.write_text(text)

        result = read_file(str(file_path))

        assert result == ([], [], [], [], meta)
        assert all(type(column) is list for column in result[:4])
        assert capsys.readouterr().out == _NO_POINTS

    def test_read_file_on_a_missing_file_raises_file_not_found(self, tmp_path):
        """U4 (v2, test advisory A4): no file is not "no points": FileNotFoundError propagates, as at the base"""
        with pytest.raises(FileNotFoundError):
            read_file(str(tmp_path / "absent.txt"))

    def test_read_file_lets_an_interrupt_through(self, tmp_path, monkeypatch):
        """U4: only a ValueError reads as "no points"; a KeyboardInterrupt raised while reading propagates"""
        file_path = tmp_path / "data.txt"
        file_path.write_text("0.1  1.0  0.1  0.01\n0.2  0.5  0.05  0.02\n")

        def interrupted(*_args, **_kwargs):
            raise KeyboardInterrupt

        monkeypatch.setattr(np, "loadtxt", interrupted)
        with pytest.raises(KeyboardInterrupt):
            read_file(str(file_path))

    def test_read_file_reads_data_rows_as_the_base_did(self, tmp_path):
        """U4: with data rows, among comments and blanks, the four columns are np.loadtxt(path).T, as at the base"""
        file_path = tmp_path / "data.txt"
        file_path.write_text(
            '# Meta:{"run": 3}\n# q r dr dq\n0.1 1.0 0.1 0.01\n\n0.2 0.5 0.05 0.02  # a note\n0.3 0.25 0.02 0.03\n'
        )

        q, r, dr, dq, meta = read_file(str(file_path))

        expected = np.loadtxt(file_path).T
        assert all(np.array_equal(got, want) for got, want in zip((q, r, dr, dq), expected, strict=True))
        assert all(isinstance(column, np.ndarray) and column.shape == (3,) for column in (q, r, dr, dq))
        assert meta == {"run": 3}

    @pytest.mark.parametrize("text", ["0.1 1.0 0.1 0.01\n", "0.1\n1.0\n0.1\n0.01\n"], ids=["one-row", "four-rows-of-one"])
    def test_read_file_with_one_row_returns_four_scalars(self, tmp_path, text):
        """Plan A2, kept as at the base: four numbers in one row, or in four rows of one, give four scalars, not
        four arrays of one element (np.loadtxt reads both as one dimension of four)"""
        file_path = tmp_path / "one.txt"
        file_path.write_text(text)

        q, r, dr, dq, _ = read_file(str(file_path))

        assert [np.ndim(v) for v in (q, r, dr, dq)] == [0, 0, 0, 0]
        assert (q, r, dr, dq) == (0.1, 1.0, 0.1, 0.01)
