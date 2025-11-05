from rich.console import Console
console = Console()

def get_engine(name: str):
    def _stub(**kwargs):
        console.print(f"[yellow]Detector '{name}' not installed yet — writing stub output if needed.[/yellow]")
    return _stub
