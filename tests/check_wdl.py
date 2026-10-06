"""Reject engine file writes in workflow expressions and task command expressions."""
import sys
from pathlib import Path
import WDL

WRITERS = {'write_lines', 'write_tsv', 'write_map', 'write_json', 'write_objects', 'write_object'}


def expression_writes(root):
    found = []
    def visit(node):
        if isinstance(node, WDL.Expr.Apply) and node.function_name in WRITERS:
            found.append((node.function_name, node.pos.line))
        for child in node.children:
            visit(child)
    visit(root)
    return found


def workflow_writes(doc):
    found = expression_writes(doc.workflow) if doc.workflow else []
    for imported in doc.imports:
        found.extend(workflow_writes(imported.doc))
    return found


def command_writes(doc):
    found = []
    for task in doc.tasks:
        found.extend(expression_writes(task.command))
    for imported in doc.imports:
        found.extend(command_writes(imported.doc))
    return found


if __name__ == '__main__':
    doc = WDL.load(sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parents[1] / 'workflows/score_variants.wdl'))
    violations = workflow_writes(doc)
    if violations:
        raise SystemExit(f'Workflow-scope file-writing functions: {violations}')
    violations = command_writes(doc)
    if violations:
        raise SystemExit(f'Engine-generated files in task command expressions: {violations}')
    print('PASS: no workflow-scope file-writing functions')
    print('PASS: task commands create their own local files')
