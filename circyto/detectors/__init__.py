from .ciri_full import run_ciri_full
from .ciri_long import run_ciri_long
from .find_circ import run_find_circ
from .circ_explorer2 import run_circ_explorer2
DETECTORS = {
    "ciri-full": run_ciri_full,
    "ciri-long": run_ciri_long,
    "find-circ": run_find_circ,
    "circ-explorer2": run_circ_explorer2,
}
def get_engine(name):
    if name not in DETECTORS:
        raise ValueError(f"Unknown engine: {name}. Available: {list(DETECTORS)}")
    return DETECTORS[name]
