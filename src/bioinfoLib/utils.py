# Copyright 2025 Zhiyuan Yu (Heemskerk's lab, University of Michigan)

import subprocess

from rich import print
from rich.tree import Tree


def print_dep_tree(max_depth=1):
    result = subprocess.run(
        ["pipdeptree", f"--depth={max_depth}"], capture_output=True, text=True
    )
    tree = Tree("[bold #FFA500]Dependencies[/bold #FFA500]")

    for line in result.stdout.splitlines():
        tree.add(f"[bold green]{line.strip()}[/bold green]")

    print(tree)
