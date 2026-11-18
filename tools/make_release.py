"""Cut a release when notes were added since the last one.

Run by .github/workflows/release.yml on the 1st of each month (or by hand:
`python tools/make_release.py --dry-run` prints what it would do and changes nothing).

A release is what makes Zenodo archive a new version under the same DOI, and what
ReachPilot announces. No new notes since the last tag: nothing happens.
"""
import datetime
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, check=True, capture_output=True, text=True, encoding="utf-8").stdout.strip()


def next_version(tag: str) -> str:
    """v1.0.1 -> 1.1.0: every batch of new notes is a minor version."""
    major, minor, _patch = (int(part) for part in re.fullmatch(r"v(\d+)\.(\d+)\.(\d+)", tag).groups())
    return f"{major}.{minor + 1}.0"


def note_title(path: Path) -> str:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# "):
            return line[2:].strip()
    return path.stem


def new_notes(tag: str) -> list[Path]:
    added = git("diff", "--name-only", "--diff-filter=A", f"{tag}..HEAD", "--", "notes")
    return [ROOT / name for name in added.splitlines() if name.endswith(".md")]


def release_body(notes: list[Path], total: int) -> str:
    lines = [f"- New note: {note_title(path)}" for path in notes]
    lines.append(f"- The repository now holds {total} notes, each with its scripts and an AI referee re-run report")
    lines.append("")
    lines.append("Computed and re-run by AI agents; not reviewed by a human mathematician. No proofs, no prize claims.")
    return "\n".join(lines)


def bump_citation(version: str, today: str) -> None:
    path = ROOT / "CITATION.cff"
    text = path.read_text(encoding="utf-8")
    text = re.sub(r'^version: .*$', f'version: "{version}"', text, count=1, flags=re.M)
    text = re.sub(r"^date-released: .*$", f"date-released: {today}", text, count=1, flags=re.M)
    path.write_text(text, encoding="utf-8", newline="\n")


def main() -> None:
    dry_run = "--dry-run" in sys.argv
    tag = git("describe", "--tags", "--abbrev=0")
    notes = new_notes(tag)
    if not notes:
        print(f"No new notes since {tag}; no release.")
        return
    version = next_version(tag)
    total = len(list((ROOT / "notes").glob("*.md")))
    body = release_body(notes, total)
    print(f"Release v{version} ({len(notes)} new note(s) since {tag}):\n{body}")
    if dry_run:
        return
    bump_citation(version, datetime.date.today().isoformat())
    git("commit", "-am", f"Release v{version}")
    git("push", "origin", "HEAD:main")
    subprocess.run(
        ["gh", "release", "create", f"v{version}", "--target", git("rev-parse", "HEAD"),
         "--title", f"v{version} - {len(notes)} new note(s)", "--notes", body],
        cwd=ROOT, check=True,
    )


if __name__ == "__main__":
    main()
