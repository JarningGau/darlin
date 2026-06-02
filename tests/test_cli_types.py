import argparse

import pytest

from darlin.cli_types import positive_float, positive_int


def test_positive_int_accepts_one() -> None:
    assert positive_int("1") == 1
    assert positive_int("42") == 42


@pytest.mark.parametrize("value", ["0", "-3", "abc"])
def test_positive_int_rejects_invalid(value: str) -> None:
    with pytest.raises(argparse.ArgumentTypeError):
        positive_int(value)


def test_positive_float_accepts_one() -> None:
    assert positive_float("1") == 1.0
    assert positive_float("0.01") == 0.01


@pytest.mark.parametrize("value", ["0", "-0.5", "abc"])
def test_positive_float_rejects_invalid(value: str) -> None:
    with pytest.raises(argparse.ArgumentTypeError):
        positive_float(value)
