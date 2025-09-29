import subprocess
import os
import sys

def run_command(cmd, work_dir):
    """
    Run a command with live output. Stops if the command fails.
    """
    print(f"\nRunning: {' '.join(cmd)}\n{'='*50}")
    process = subprocess.Popen(
        cmd,
        cwd=work_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    # Print output line by line
    for line in process.stdout:
        print(line, end="")

    process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")

def build_paper():
    # Use the directory where this script is located
    work_dir = os.path.dirname(os.path.abspath(__file__))

    # List of commands
    commands = [
    ["quarto", "render", "paper.qmd"],  # stop at .tex
    ["python", "fix_theorems.py", "paper.tex"],
    ["xelatex","-interaction=nonstopmode", "fixed_paper.tex"],
    ["bibtex", "fixed_paper"],
    ["xelatex", "-interaction=nonstopmode", "fixed_paper.tex"]
]

    try:
        for cmd in commands:
            run_command(cmd, work_dir)
        print("\n✅ Build finished successfully!")
    except RuntimeError as e:
        print(f"\n❌ Build failed:\n{e}")
        sys.exit(1)

if __name__ == "__main__":
    build_paper()
