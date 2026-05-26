"""Package version lookup with a local fallback for editable checkouts."""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("pydft")
except PackageNotFoundError:
    __version__ = "0.0.0"  # optional fallback
