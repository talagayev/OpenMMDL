import argparse

import pytest

from openmmdl.openmmdl_analysis.openmmdlanalysis import parse_bool_flag, parse_xyz


@pytest.mark.parametrize(
    "value, expected",
    [
        (True, True),
        (False, False),
    ],
)
def test_parse_bool_flag_returns_bool_inputs_unchanged(value, expected):
    assert parse_bool_flag(value) is expected


@pytest.mark.parametrize(
    "value",
    [
        "1",
        1,
        "true",
        "TRUE",
        " t ",
        "yes",
        "Y",
        " on ",
    ],
)
def test_parse_bool_flag_accepts_truthy_values(value):
    assert parse_bool_flag(value) is True


@pytest.mark.parametrize(
    "value",
    [
        "0",
        0,
        "false",
        "FALSE",
        " f ",
        "no",
        "N",
        " off ",
    ],
)
def test_parse_bool_flag_accepts_falsy_values(value):
    assert parse_bool_flag(value) is False


@pytest.mark.parametrize(
    "value",
    ["", "maybe", "2", "enable", None],
)
def test_parse_bool_flag_rejects_invalid_values(value):
    with pytest.raises(argparse.ArgumentTypeError) as exc_info:
        parse_bool_flag(value)

    assert str(exc_info.value) == (
        "Boolean value expected. Use one of: true/false, yes/no, y/n, 1/0, on/off."
    )


@pytest.mark.parametrize(
    "value, expected",
    [
        ("(1.0, 2.0, 3.0)", [1.0, 2.0, 3.0]),
        ("[1.0, 2.0, 3.0]", [1.0, 2.0, 3.0]),
        (
            "(np.float64(30.156), np.float64(41.233), np.float64(92.595))",
            [30.156, 41.233, 92.595],
        ),
    ],
)
def test_parse_xyz_parses_plain_and_numpy_float_coordinates(value, expected):
    assert parse_xyz(value) == expected


def test_parse_xyz_does_not_parse_float64_width_as_coordinate():
    parsed = parse_xyz(
        "(np.float64(30.156), np.float64(41.233), np.float64(92.595))"
    )

    assert parsed == [30.156, 41.233, 92.595]
    assert 64.0 not in parsed


@pytest.mark.parametrize("value", [None, 0, "0"])
def test_parse_xyz_returns_none_for_missing_coordinates(value):
    assert parse_xyz(value) is None