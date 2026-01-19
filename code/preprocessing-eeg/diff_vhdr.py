
import sys
from difflib import unified_diff
from rich.console import Console
from rich.table import Table
from rich.syntax import Syntax

"""
Script to compare two BrainVision header (.vhdr) files.

This script displays a rich, color-coded unified diff of two text files,
intended for comparing VHDR files.

Usage:
    python diff_vhdr.py <file1.vhdr> <file2.vhdr>
"""

def print_diff_rich(file1_path, file2_path):
    """
    Prints a color-coded diff of two files to the console.

    Args:
        file1_path (str): Path to the first file.
        file2_path (str): Path to the second file.
    """
    console = Console()

    try:
        with open(file1_path, 'r', encoding='utf-8') as f1, open(file2_path, 'r', encoding='utf-8') as f2:
            f1_lines = f1.readlines()
            f2_lines = f2.readlines()
    except FileNotFoundError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        return

    # Generate diff
    diff = list(unified_diff(
        f1_lines, f2_lines,
        fromfile=file1_path, tofile=file2_path,
        lineterm=''
    ))

    if not diff:
        console.print("[bold green]Files are identical.[/bold green]")
        return

    # Create a table for cleaner output
    table = Table(title=f"Diff: {file1_path} vs {file2_path}", show_lines=False)
    table.add_column("Type", style="bold", width=4)
    table.add_column("Content", style="dim")

    for line in diff:
        if line.startswith('---') or line.startswith('+++') or line.startswith('@@'):
            table.add_row("META", f"[cyan]{line}[/cyan]")
        elif line.startswith('+'):
            # Strip the + for the content, keep the color green
            table.add_row("+", f"[green]{line[1:].rstrip()}[/green]")
        elif line.startswith('-'):
            table.add_row("-", f"[red]{line[1:].rstrip()}[/red]")
        else:
            # Context lines
            table.add_row(" ", line.rstrip())

    console.print(table)

if __name__ == "__main__":
    # Check if the user provided enough arguments
    if len(sys.argv) < 3:
        print("Usage: python vdiff.py <file1.vhdr> <file2.vhdr>")
        sys.exit(1)

    # sys.argv[0] is the script name, [1] is first file, [2] is second file
    print_diff_rich(sys.argv[1], sys.argv[2])
