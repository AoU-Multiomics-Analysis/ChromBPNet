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
                    if isinstance(node, WDL.Tree.Call) and node.name == 'MergeVariantEffects')
        env = WDL.Env.Bindings()
        for name, value in [('merge_memory_gb', 96), ('merge_disk_gb', 500), ('merge_max_retries', 3)]:
            env = env.bind(name, WDL.Value.Int(value))
        task_env = WDL.Env.Bindings()
        for name in ['memory_gb', 'disk_gb', 'max_retries']:
            task_env = task_env.bind(name, call.inputs[name].eval(env, LocalStdLib('1.0')))
        runtime = {name: call.callee.runtime[name].eval(task_env, LocalStdLib('1.0')).value
                   for name in ['memory', 'disks', 'maxRetries']}
        self.assertEqual(runtime, {'memory': '96GB', 'disks': 'local-disk 500 SSD', 'maxRetries': 3})

    def test_merge_only_workflow_has_no_scoring_calls_and_preserves_file_inputs(self):
        doc = WDL.load(str(ROOT / 'workflows/merge_variants.wdl'))
        calls = [node for node in doc.workflow.body if isinstance(node, WDL.Tree.Call)]
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0].callee.name, 'MergeVariantEffects')
        declarations = {decl.name: decl for decl in doc.workflow.inputs}
        self.assertEqual(str(declarations['score_files'].type), 'Array[File]')
        for name, value in [('merge_memory_gb', 64), ('merge_disk_gb', 500), ('merge_max_retries', 2)]:
            self.assertEqual(declarations[name].expr.eval(WDL.Env.Bindings(), LocalStdLib('1.0')).value, value)
        cloud_files = WDL.Value.Array(WDL.Type.File(), [WDL.Value.File('gs://previous-run/fold0.tsv')])
        workflow_env = WDL.Env.Bindings().bind('score_files', cloud_files)
        files = calls[0].inputs['score_files'].eval(workflow_env, LocalStdLib('1.0'))
        self.assertIsInstance(files.value[0], WDL.Value.File)
        task_env = WDL.Env.Bindings().bind('score_files', files)
        task_env = WDL.Value.rewrite_env_paths(task_env, lambda file: '/localized/fold0.tsv')
        command = calls[0].callee.command.eval(task_env, LocalStdLib('1.0')).value
        self.assertIn('/localized/fold0.tsv', command)
        self.assertNotIn('gs://', command)
        self.assertEqual(workflow_writes(doc), [])
        self.assertEqual(command_writes(doc), [])
        self.assertEqual(unsupported_functions(doc), [])


if __name__ == '__main__':
    unittest.main()
