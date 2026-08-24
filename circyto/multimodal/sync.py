from __future__ import annotations

from contextlib import contextmanager
from os import PathLike
from typing import Any, Iterator, Mapping


MUDATA_SYNC_POLICY = "pull_on_update=True"


@contextmanager
def _explicit_pull_on_update() -> Iterator[Any]:
    """Use MuData's historical automatic pull policy for one operation."""
    import mudata as mu

    settings = getattr(mu, "settings", None)
    override = getattr(settings, "override", None)
    if override is not None:
        with override(pull_on_update=True):
            yield mu
    else:
        with mu.set_options(pull_on_update=True):
            yield mu


def mudata_from_modalities(modalities: Mapping[str, Any], **kwargs: Any) -> Any:
    """Construct MuData while explicitly preserving automatic obs/var pulls."""
    with _explicit_pull_on_update() as mu:
        return mu.MuData(modalities, **kwargs)


def read_h5mu(path: str | PathLike[str], **kwargs: Any) -> Any:
    """Read H5MU while explicitly preserving automatic obs/var pulls."""
    with _explicit_pull_on_update() as mu:
        return mu.read_h5mu(path, **kwargs)


def write_h5mu(mdata: Any, path: str | PathLike[str], **kwargs: Any) -> None:
    """Write H5MU while explicitly preserving automatic obs/var pulls."""
    with _explicit_pull_on_update():
        mdata.write_h5mu(path, **kwargs)


def synchronize_mudata(mdata: Any) -> None:
    """Update MuData indices and pull obs/var columns using circyto's policy."""
    with _explicit_pull_on_update():
        mdata.update()


def refresh_prefixed_var_from_modality(mdata: Any, modality: str) -> None:
    """Replace stale modality-prefixed global var columns with explicit pulls."""
    modality_var = mdata.mod[modality].var
    prefixed_columns = [f"{modality}:{column}" for column in modality_var.columns]
    stale_columns = [column for column in prefixed_columns if column in mdata.var.columns]
    if stale_columns:
        mdata.var = mdata.var.drop(columns=stale_columns)
    if len(modality_var.columns) > 0:
        mdata.pull_var(mods=[modality])
