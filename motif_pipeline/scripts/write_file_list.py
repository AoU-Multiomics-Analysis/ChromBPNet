"""Create a task-local list after WDL localizes and safely quotes each File."""
import argparse
from pathlib import Path

from motif_io import readable


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    parser.add_argument('--file', action='append', default=[])
    args = parser.parse_args()
    if len(args.file) != 5:
        raise ValueError('Expected exactly five localized fold files')
    # Validate every argument before serialization. CR/LF cannot enter a list.
    paths = [readable(value) for value in args.file]
    Path(args.output).write_text(''.join(str(path) + '\n' for path in paths))
    print(f'[files] Wrote five localized paths to {args.output}', flush=True)


if __name__ == '__main__':
    main()
