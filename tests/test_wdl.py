import shlex
import tempfile
import unittest
from pathlib import Path
import WDL
from check_wdl import command_writes, unsupported_functions, workflow_writes

ROOT = Path(__file__).resolve().parents[1]


class LocalStdLib(WDL.StdLib.Base):
    def _devirtualize_filename(self, filename):
        return filename

    def _virtualize_filename(self, filename):
        return filename


class CloudGeneratedFileStdLib(LocalStdLib):
    """Model the cloud path returned by Terra for an engine-generated file."""
    def _virtualize_filename(self, filename):
        return 'gs://test-bucket/call-MergeVariantEffects/' + Path(filename).name


class WDLTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = WDL.load(str(ROOT / 'workflows/score_variants.wdl'))

    def test_workflow_has_no_file_writes(self):
        self.assertEqual(workflow_writes(self.doc), [])

    def test_workflow_uses_only_wdl_1_0_functions(self):
        self.assertEqual(unsupported_functions(self.doc), [])

    def test_static_check_rejects_sep_function_but_accepts_separator_option(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'separator.wdl'
            path.write_text('''version 1.0
task bad { input { Array[File] xs } command <<< echo '~{sub(sep(" ", xs), "x", "y")}' >>> }
''')
            self.assertEqual(unsupported_functions(WDL.load(str(path))), [('sep', 2)])
            path.write_text('''version 1.0
task good { input { Array[File] xs } command <<< echo '~{sep=" " xs}' >>> }
''')
            self.assertEqual(unsupported_functions(WDL.load(str(path))), [])

    def test_static_check_catches_nested_workflow_writes(self):
        source = 'version 1.0\nworkflow bad { input { Array[String] xs } scatter (x in xs) { File f = write_lines([x]) } output { Array[File] files = f } }'
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'bad.wdl'
            path.write_text(source)
            self.assertIn(('write_lines', 2), workflow_writes(WDL.load(str(path))))

    def test_commands_have_no_engine_generated_file_writes(self):
        self.assertEqual(command_writes(self.doc), [])

    def test_static_check_catches_command_writes_inside_string_functions(self):
        source = '''version 1.0
task bad { input { Array[File] xs } command <<< cat '~{sub(write_lines(xs), "x", "y")}' >>> }
'''
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'bad.wdl'
            path.write_text(source)
            self.assertIn(('write_lines', 2), command_writes(WDL.load(str(path))))

    def task(self, name):
        return next(task for task in self.doc.tasks if task.name == name)

    def test_manifest_read_retains_file_types(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / 'models.normalized.tsv'
            manifest.write_text("cd4\tCD4 T cell\tgs://b/a ' model.h5\tgs://b/p.bed\n")
            env = WDL.Env.Bindings().bind('ValidateManifest.rows', WDL.Value.File(str(manifest)))
            decl = next(item for item in self.doc.workflow.body if isinstance(item, WDL.Tree.Decl) and item.name == 'manifest_rows')
            rows = decl.expr.eval(env, LocalStdLib('1.0', tmp)).coerce(decl.type)
            scatter = next(item for item in self.doc.workflow.body if isinstance(item, WDL.Tree.Scatter))
            entry = next(item for item in scatter.body if isinstance(item, WDL.Tree.Decl) and item.name == 'entry')
            value = entry.expr.eval(WDL.Env.Bindings().bind('row', rows.value[0]), LocalStdLib('1.0', tmp)).coerce(entry.type)
            for member in ['model', 'peaks']:
                self.assertIsInstance(value.value[member], WDL.Value.File)
                self.assertTrue(value.value[member].value.startswith('gs://'))
            self.assertEqual(value.value['model'].value, "gs://b/a ' model.h5")

    def test_score_command_uses_localized_paths_and_safe_shell_quoting(self):
        task = self.task('ScoreVariants')
        suffix = "' $(touch NEVER) `touch NEVER` file"
        args = {name: '/localized/' + name + suffix for name in ['model', 'peaks', 'variants', 'genome', 'genome_index', 'chrom_sizes', 'manifest_validation']}
        args.update(model_id='cd4', cell_type="CD4 T ' cell $(touch NEVER)", docker_image='test', batch_size=2,
                    num_shuf=0, max_peaks=None, random_seed=1234, memory_gb=64, disk_gb=100, num_preempt=0)
        env = WDL.values_from_json(args, task.available_inputs, task.required_inputs)
        command = task.command.eval(env, LocalStdLib('1.0')).value
        # Cromwell uses Java regex replacement rules. No replacement backslash is allowed.
        self.assertNotIn("\\'", command)
        # Extract only the Python invocation, without the log redirections.
        invocation = command[command.index('python /opt/'):command.index("echo '[score] Variant effects complete'")]
        words = shlex.split(invocation.replace('\\\n', ''))
        for option in ['model', 'peaks', 'variants', 'genome', 'genome-index', 'chrom-sizes', 'cell-type']:
            self.assertEqual(words[words.index('--' + option) + 1], args[option.replace('-', '_')])
        self.assertNotIn('--max-peaks', words)
        env = WDL.values_from_json(dict(args, max_peaks=100), task.available_inputs, task.required_inputs)
        command = task.command.eval(env, LocalStdLib('1.0')).value
        self.assertIn('--max-peaks 100', command)
        env = WDL.values_from_json(dict(args, genome_index=None), task.available_inputs, task.required_inputs)
        command = task.command.eval(env, LocalStdLib('1.0')).value
        self.assertNotIn('--genome-index', command)




if __name__ == '__main__':
    unittest.main()
