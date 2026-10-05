import json
import shlex
import tempfile
import unittest
from pathlib import Path
import WDL
from check_wdl import workflow_writes

ROOT = Path(__file__).resolve().parents[1]


class LocalStdLib(WDL.StdLib.Base):
    def _devirtualize_filename(self, filename):
        return filename

    def _virtualize_filename(self, filename):
        return filename


class WDLTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = WDL.load(str(ROOT / 'workflows/score_variants.wdl'))

    def test_workflow_has_no_file_writes(self):
        self.assertEqual(workflow_writes(self.doc), [])

    def test_static_check_catches_nested_workflow_writes(self):
        source = 'version 1.0\nworkflow bad { input { Array[String] xs } scatter (x in xs) { File f = write_lines([x]) } output { Array[File] files = f } }'
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'bad.wdl'
            path.write_text(source)
            self.assertIn(('write_lines', 2), workflow_writes(WDL.load(str(path))))

    def task(self, name):
        return next(task for task in self.doc.tasks if task.name == name)

    def test_manifest_read_retains_file_types(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / 'models.json'
            manifest.write_text(json.dumps([dict(model_id='cd4', cell_type='CD4', model='gs://b/m.h5', peaks='gs://b/p.bed')]))
            env = WDL.Env.Bindings().bind('model_manifest', WDL.Value.File(str(manifest)))
            decl = next(item for item in self.doc.workflow.body if isinstance(item, WDL.Tree.Decl) and item.name == 'models')
            value = decl.expr.eval(env, LocalStdLib('1.0', tmp)).coerce(decl.type)
            for member in ['model', 'peaks']:
                self.assertIsInstance(value.value[0].value[member], WDL.Value.File)
                self.assertTrue(value.value[0].value[member].value.startswith('gs://'))

    def test_score_command_uses_localized_paths_and_safe_shell_quoting(self):
        task = self.task('ScoreVariants')
        suffix = "' $(touch NEVER) `touch NEVER` file"
        args = {name: '/localized/' + name + suffix for name in ['model', 'peaks', 'variants', 'genome', 'genome_index', 'chrom_sizes', 'manifest_validation']}
        args.update(model_id='cd4', cell_type="CD4 T ' cell $(touch NEVER)", docker_image='test', batch_size=2,
                    num_shuf=0, max_peaks=None, random_seed=1234, memory_gb=64, disk_gb=100, num_preempt=0)
        env = WDL.values_from_json(args, task.available_inputs, task.required_inputs)
        command = task.command.eval(env, LocalStdLib('1.0')).value
        # Extract only the Python invocation, without the log redirections.
        invocation = command[command.index('python /opt/'):command.index("echo '[score] Variant effects complete'")]
        words = shlex.split(invocation.replace('\\\n', ''))
        for option in ['model', 'peaks', 'variants', 'genome', 'genome-index', 'chrom-sizes', 'cell-type']:
            self.assertEqual(words[words.index('--' + option) + 1], args[option.replace('-', '_')])
        self.assertNotIn('--max-peaks', words)
        env = WDL.values_from_json(dict(args, max_peaks=100), task.available_inputs, task.required_inputs)
        command = task.command.eval(env, LocalStdLib('1.0')).value
        self.assertIn('--max-peaks 100', command)

    def test_merge_file_list_uses_command_time_localized_files(self):
        task = self.task('MergeVariantEffects')
        with tempfile.TemporaryDirectory() as tmp:
            paths = ['/localized/cd4 scores.tsv', "/localized/nk ' scores.tsv"]
            env = WDL.values_from_json(dict(score_files=paths, docker_image='test'), task.available_inputs, task.required_inputs)
            command = task.command.eval(env, LocalStdLib('1.0', tmp)).value
            invocation = command[command.index('python /opt/'):command.index("echo '[merge] Cell-type merge complete'")]
            words = shlex.split(invocation.replace('\\\n', ''))
            file_list = Path(words[words.index('--score-files') + 1])
            self.assertEqual(file_list.read_text().splitlines(), paths)
            self.assertNotIn('gs://', file_list.read_text())


if __name__ == '__main__':
    unittest.main()
