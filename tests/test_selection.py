"""mcpymol.selection: chain and residue identifiers as data, not syntax."""

import pytest

from mcpymol.selection import chain_clause, derived_name, quote


@pytest.mark.parametrize(
    "chains, expected",
    [
        ("A", 'chain "A"'),
        ("", 'chain ""'),  # a blank chain stays scoped to its object
        ("A+B", 'chain "A"+"B"'),  # a caller's list still means two chains
        (" A + B ", 'chain "A"+"B"'),
        (["", "A"], 'chain ""+"A"'),
        (("B", "A"), 'chain "B"+"A"'),
        ("-", 'chain "\\-"'),
    ],
)
def test_chain_clause(chains, expected):
    assert chain_clause(chains) == expected


def test_a_double_quote_is_refused_not_guessed_at():
    with pytest.raises(ValueError, match="double quote"):
        chain_clause('A"B')


def test_wiggles_uses_the_same_quote():
    from mcpymol.wiggles import atoms

    assert atoms.quote is quote


@pytest.mark.parametrize(
    ("base", "expected"),
    [
        ("1ubq", "1ubq_hbonds"),  # a plain name is unchanged
        ("AF_P04637-2", "AF_P04637-2_hbonds"),  # '-', '.', '+' survive in PyMOL names
        ("keep1 and chain A", "keep1_and_chain_A_hbonds"),  # one name, not a delete list
        ("1ubq or 4hhb", "1ubq_or_4hhb_hbonds"),
        ("a;b", "a_b_hbonds"),  # no second command in a `do` line
        ("*", "hbonds"),  # no wildcard, and not hidden behind a leading '_'
        ("(sele)", "sele__hbonds"),
        ("_hidden", "_hidden_hbonds"),  # a hidden object's helpers stay hidden
    ],
)
def test_derived_name_is_always_one_plain_name(base, expected):
    assert derived_name(base, "hbonds") == expected
