"""Create a task-local list from localized, individually quoted File arguments."""
import argparse
from pathlib import Path
from finemo_io import readable

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', action='append', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    files = [str(readable(value)) for value in args.file]
    Path(args.output).write_text(''.join(value + '\n' for value in files))

if __name__ == '__main__':
    main()
