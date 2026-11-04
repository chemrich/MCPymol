"""mcpymol.selection: chain and residue identifiers as data, not syntax."""

import pytest

from mcpymol.selection import chain_clause, quote


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
