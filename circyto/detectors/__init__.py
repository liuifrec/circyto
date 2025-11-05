from rich.console import Console

console = Console()

def get_engine(name: str):
    """
    Return a callable detector runner for the given engine name.
    Placeholder until individual engines are wrapped.
    """
    def _dummy_detector(**kwargs):
        console.print(f"[yellow]Detector '{name}' not yet implemented.[/yellow]")
        console.print(kwargs)
    return _dummy_detector
