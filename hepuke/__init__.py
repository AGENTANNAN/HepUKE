"""HepUKE public API.

    from hepuke import HepUKE, Document, load_config

`Config` / `load_config` / `Document` import eagerly; `HepUKE` is resolved
lazily so that importing the package does not pull in the vector-store and
embedding stack until you actually build a client.
"""
from hepuke.config import Config, load_config
from hepuke.document import Document

__all__ = ["HepUKE", "Config", "load_config", "Document"]


def __getattr__(name):  # PEP 562
    if name == "HepUKE":
        from hepuke.core.index import HepUKE
        return HepUKE
    raise AttributeError(f"module 'hepuke' has no attribute {name!r}")
