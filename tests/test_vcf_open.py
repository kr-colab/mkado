"""Resource cleanup when opening VCF inputs."""

import errno
import logging
import os

import cyvcf2
import pytest

from mkado.io.vcf import _open_vcf


@pytest.mark.parametrize("fails", [False, True])
def test_open_vcf_closes_pipe_and_restores_stderr(monkeypatch, caplog, fails):
    pipe_fds = []
    real_pipe = os.pipe
    original_stderr = os.fstat(2)
    opened = object()
    failure = ValueError("invalid VCF")

    def tracked_pipe():
        pair = real_pipe()
        pipe_fds.extend(pair)
        return pair

    def open_input(path):
        assert path == "input.vcf"
        os.write(2, b"test htslib warning\n")
        if fails:
            raise failure
        return opened

    monkeypatch.setattr(os, "pipe", tracked_pipe)
    monkeypatch.setattr(cyvcf2, "VCF", open_input)
    caplog.set_level(logging.DEBUG, logger="mkado.io.vcf")
    try:
        if fails:
            with pytest.raises(ValueError) as caught:
                _open_vcf("input.vcf")
            assert caught.value is failure
        else:
            assert _open_vcf("input.vcf") is opened
            assert "htslib: test htslib warning" in caplog.text
        restored = os.fstat(2)
        assert (restored.st_dev, restored.st_ino) == (
            original_stderr.st_dev,
            original_stderr.st_ino,
        )
        for fd in pipe_fds:
            with pytest.raises(OSError) as closed:
                os.fstat(fd)
            assert closed.value.errno == errno.EBADF
    finally:
        # Keep a failing regression run from leaking descriptors itself.
        for fd in pipe_fds:
            try:
                os.close(fd)
            except OSError:
                pass
