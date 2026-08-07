"""Scene conventions shared across the view presets.

Small, but worth having in one place: every preset that set a background was
setting it *almost* right, and the failure only shows up in a rendered image
rather than in the viewport.
"""

from mcpymol.bridge import send_request


def set_background(color: str = "black") -> None:
    """Set the viewport background, and make it survive a render.

    ``bg_color`` alone is not enough. Ray-traced output keeps an alpha channel
    unless ``opaque_background`` is set, so a preset that sets only
    ``bg_color`` looks correct in the viewport and renders transparent —
    which appears white wherever the PNG is actually used. Every preset except
    the ghost-heart style had that bug, including ``textbook_view``, whose
    background is deliberately white rather than black.
    """
    send_request("do", args=[f"bg_color {color}"])
    send_request("set", args=["opaque_background", "1"])


def black_background() -> None:
    """The house background for analytical views."""
    set_background("black")
