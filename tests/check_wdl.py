"""Check Terra file creation and the WDL 1.0 standard library."""
import sys
from pathlib import Path
import WDL

WRITERS = {'write_lines', 'write_tsv', 'write_map', 'write_json', 'write_objects', 'write_object'}
# https://github.com/openwdl/wdl/blob/wdl-1.0/SPEC.md#standard-library
# miniwdl accepts some newer functions even when the document declares 1.0.
WDL_1_0_FUNCTIONS = WRITERS | {
    'stdout', 'stderr', 'glob', 'read_lines', 'read_tsv', 'read_map',
    'read_object', 'read_objects', 'read_json', 'read_int', 'read_string',
    'read_float', 'read_boolean', 'size', 'sub', 'range', 'transpose',
    'zip', 'cross', 'length', 'flatten', 'prefix', 'select_first',
    'select_all', 'defined', 'basename', 'floor', 'ceil', 'round',
}


def unsupported_functions(doc):
    found = []
    def visit(node):
        if isinstance(node, WDL.Expr.Apply):
            name = node.function_name
            # miniwdl represents operators with internal names such as _add.
            if not name.startswith('_') and name not in WDL_1_0_FUNCTIONS:
                found.append((name, node.pos.line))
        for child in node.children:
            visit(child)
    visit(doc)
    for imported in doc.imports:
        found.extend(unsupported_functions(imported.doc))
    return found


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
    violations = unsupported_functions(doc)
    if violations:
        raise SystemExit(f'Functions outside the WDL 1.0 standard library: {violations}')
    violations = workflow_writes(doc)
    if violations:
        raise SystemExit(f'Workflow-scope file-writing functions: {violations}')
    violations = command_writes(doc)
    if violations:
        raise SystemExit(f'Engine-generated files in task command expressions: {violations}')
    print('PASS: no workflow-scope file-writing functions')
    print('PASS: task commands create their own local files')
    print('PASS: all functions belong to the WDL 1.0 standard library')
