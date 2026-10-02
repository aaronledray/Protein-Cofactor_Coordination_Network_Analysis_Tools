"""Decide whether figures should be shown on screen.

The legacy plotting functions call ``plt.show()`` and plotly's ``fig.show()``.
Those open windows or a browser tab and block until they are dismissed or
fetched, which hangs headless, piped, cron, CI, and subprocess runs. By default
figures are shown only when both stdin and stdout are terminals; the files the
tools write are identical either way. ``set_show_figures(True/False)`` (the
CLI's ``--show`` / ``--no-show``) overrides the automatic choice.
"""

import sys
from typing import Any, Optional

_override: Optional[bool] = None


def set_show_figures(value: Optional[bool]) -> None:
    """Force showing (True), suppress it (False), or restore auto (None)."""
    global _override
    _override = value


def figures_should_show() -> bool:
    if _override is not None:
        return _override
    try:
        return bool(sys.stdin.isatty() and sys.stdout.isatty())
    except (AttributeError, ValueError):  # closed or replaced streams
        return False


def show_matplotlib(plt: Any) -> None:
    """``plt.show()`` when allowed; otherwise release the current figure."""
    if figures_should_show():
        plt.show()
    elif plt.get_fignums():
        plt.close(plt.gcf())


def show_plotly(figure: Any) -> None:
    """``figure.show()`` when allowed; otherwise do nothing."""
    if figures_should_show():
        figure.show()
