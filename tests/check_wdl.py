"""Reject file-writing functions in all workflow expressions, including nested bodies."""
import sys
from pathlib import Path
import WDL

WRITERS = {'write_lines', 'write_tsv', 'write_map', 'write_json', 'write_objects', 'write_object'}


def workflow_writes(doc):
    found = []
    def visit(node):
        if isinstance(node, WDL.Expr.Apply) and node.function_name in WRITERS:
            found.append((node.function_name, node.pos.line))
        for child in node.children:
            visit(child)
    if doc.workflow:
        visit(doc.workflow)
    for imported in doc.imports:
        found.extend(workflow_writes(imported.doc))
    return found


if __name__ == '__main__':
    doc = WDL.load(sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parents[1] / 'workflows/score_variants.wdl'))
    violations = workflow_writes(doc)
    if violations:
        raise SystemExit(f'Workflow-scope file-writing functions: {violations}')
    print('PASS: no workflow-scope file-writing functions')
