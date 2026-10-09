import unittest
from pathlib import Path
import WDL
from test_wdl import LocalStdLib
from check_wdl import command_writes, unsupported_functions, workflow_writes

ROOT = Path(__file__).resolve().parents[1]


class MergeRecoveryTests(unittest.TestCase):
    def test_merge_resources_and_retries_reach_task_runtime(self):
        doc = WDL.load(str(ROOT / 'workflows/score_variants.wdl'))
        call = next(node for node in doc.workflow.body
                    if isinstance(node, WDL.Tree.Call) and node.name == 'SummarizeVariants')
        env = WDL.Env.Bindings()
        values = [('merge_memory_gb', 96), ('merge_disk_gb', 500), ('merge_max_retries', 3),
                  ('summary_memory_gb', 32), ('summary_disk_gb', 500), ('summary_max_retries', 4)]
        for name, value in values:
            env = env.bind(name, WDL.Value.Int(value))
        shared_env = WDL.Env.Bindings()
        for name, value in values:
            shared_env = shared_env.bind(name, call.inputs[name].eval(env, LocalStdLib('1.0')))
        shared = call.callee
        scatter = next(node for node in shared.body if isinstance(node, WDL.Tree.Scatter))
        summary = next(node for node in scatter.body if isinstance(node, WDL.Tree.Call))
        combine = next(node for node in shared.body if isinstance(node, WDL.Tree.Call)
                       and node.name == 'CombineModelEffects')
        for task_call, memory, retries in [(summary, 32, 4), (combine, 96, 3)]:
            task_env = WDL.Env.Bindings()
            for name in ['memory_gb', 'disk_gb', 'max_retries']:
                task_env = task_env.bind(name, task_call.inputs[name].eval(shared_env, LocalStdLib('1.0')))
            runtime = {name: task_call.callee.runtime[name].eval(task_env, LocalStdLib('1.0')).value
                       for name in ['memory', 'disks', 'maxRetries']}
            self.assertEqual(runtime, {'memory': f'{memory}GB', 'disks': 'local-disk 500 SSD', 'maxRetries': retries})

    def test_merge_only_workflow_has_no_scoring_calls_and_preserves_file_inputs(self):
        doc = WDL.load(str(ROOT / 'workflows/merge_variants.wdl'))
        calls = [node for node in doc.workflow.body if isinstance(node, WDL.Tree.Call)]
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].callee.name, 'ChromBPNetSummarizeVariants')
        def callees(workflow):
            result = []
            def visit(node):
                if isinstance(node, WDL.Tree.Call):
                    result.append(node.callee.name)
                    if isinstance(node.callee, WDL.Tree.Workflow):
                        result.extend(callees(node.callee))
                elif isinstance(node, WDL.Tree.Scatter):
                    for child in node.body:
                        visit(child)
            for node in workflow.body:
                visit(node)
            return result
        self.assertNotIn('ScoreVariants', callees(doc.workflow))
        declarations = {decl.name: decl for decl in doc.workflow.inputs}
        self.assertEqual(str(declarations['score_files'].type), 'Array[File]')
        self.assertEqual(str(declarations['model_manifest'].type), 'File')
        for name, value in [('merge_memory_gb', 64), ('merge_disk_gb', 500), ('merge_max_retries', 2)]:
            self.assertEqual(declarations[name].expr.eval(WDL.Env.Bindings(), LocalStdLib('1.0')).value, value)
        cloud_files = WDL.Value.Array(WDL.Type.File(), [WDL.Value.File('gs://previous-run/fold0.tsv')])
        workflow_env = WDL.Env.Bindings().bind('score_files', cloud_files)
        files = calls[0].inputs['score_files'].eval(workflow_env, LocalStdLib('1.0'))
        self.assertIsInstance(files.value[0], WDL.Value.File)
        shared_doc = doc.imports[0].doc
        task = next(t for t in shared_doc.tasks if t.name == 'SummarizeModel')
        task_env = WDL.values_from_json(dict(score_files=['gs://previous-run/fold0.tsv'],
            peaks='gs://previous-run/peaks.bed', model_ids=['a_fold0'], model_group='a',
            cell_type='CD4', docker_image='test'), task.available_inputs, task.required_inputs)
        task_env = WDL.Value.rewrite_env_paths(task_env, lambda file: '/localized/fold0.tsv')
        command = task.command.eval(task_env, LocalStdLib('1.0')).value
        self.assertIn('/localized/fold0.tsv', command)
        self.assertNotIn('gs://', command)
        self.assertEqual(workflow_writes(doc), [])
        self.assertEqual(command_writes(doc), [])
        self.assertEqual(unsupported_functions(doc), [])


if __name__ == '__main__':
    unittest.main()
