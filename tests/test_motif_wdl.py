import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import WDL
from check_wdl import command_writes, unsupported_functions, workflow_writes
from test_wdl import CloudGeneratedFileStdLib, LocalStdLib
import test_motif_pipeline

ROOT = Path(__file__).resolve().parents[1]


class MotifWDLTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = WDL.load(str(ROOT / 'workflows/discover_motifs.wdl'))

    def task(self, name):
        return next(task for task in self.doc.tasks if task.name == name)

    def test_terra_wdl_10_and_no_engine_file_writes(self):
        self.assertEqual(self.doc.wdl_version, '1.0')
        self.assertEqual(workflow_writes(self.doc), [])
        self.assertEqual(command_writes(self.doc), [])
        self.assertEqual(unsupported_functions(self.doc), [])

    def test_required_inputs_keep_file_types(self):
        inputs = {decl.name: decl.type for decl in self.doc.workflow.inputs}
        self.assertIsInstance(inputs['fold_models'], WDL.Type.Array)
        self.assertIsInstance(inputs['fold_models'].item_type, WDL.Type.File)
        for name in ['peaks', 'genome', 'genome_index', 'chrom_sizes', 'motif_database']:
            self.assertIsInstance(inputs[name], WDL.Type.File)

    def test_contribution_command_quotes_paths_and_has_no_json_argument_wrapper(self):
        task = self.task('FoldContributions')
        suffix = "' $(touch NEVER) `touch NEVER` file"
        args = dict(model='/local/model' + suffix, peaks='/local/peaks' + suffix,
                    genome='/local/genome' + suffix, genome_index='/local/index' + suffix,
                    preparation='/local/preparation' + suffix, fold_index=0,
                    cell_type="CD4 ' $(touch NEVER)", head='counts', random_seed=1234,
                    num_backgrounds=20, batch_size=16, threads=16,
                    docker_image='test', memory_gb=64, disk_gb=100, num_preempt=0)
        env = WDL.values_from_json(args, task.available_inputs, task.required_inputs)
        command = task.command.eval(env, LocalStdLib('1.0')).value
        invocation = command[command.index('python /opt/'):command.index("echo '[contributions] Complete'")]
        words = shlex.split(invocation.replace('\\\n', ''))
        for option in ['model', 'peaks', 'genome', 'genome-index', 'cell-type', 'preparation']:
            self.assertEqual(words[words.index('--' + option) + 1], args[option.replace('-', '_')])
        self.assertNotIn('--config', words)

    def test_average_runs_real_script_after_cloud_to_local_mapping(self):
        task = self.task('AverageContributions')
        fixture = test_motif_pipeline.MotifPipelineTests(methodName='test_average_is_arithmetic_mean_with_float32_output')
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        files, beds = zip(*(fixture.fixture(i) for i in range(5)))
        cloud_scores = [f'gs://bucket/shard-{i}/scores.h5' for i in range(5)]
        cloud_beds = [f'gs://bucket/shard-{i}/regions.bed' for i in range(5)]
        mapping = dict(zip(cloud_scores + cloud_beds, map(str, files + beds)))
        args = dict(score_files=cloud_scores, region_files=cloud_beds, docker_image='test',
                    expected_folds=5, chunk_rows=1, memory_gb=8, disk_gb=100)
        env = WDL.values_from_json(args, task.available_inputs, task.required_inputs)
        env = WDL.Value.rewrite_env_paths(env, lambda file: mapping[file.value])
        command = task.command.eval(env, CloudGeneratedFileStdLib('1.0', fixture.root)).value
        entrypoint = shlex.quote(sys.executable) + ' ' + shlex.quote(str(ROOT / 'motif_pipeline/scripts/average_contributions.py'))
        command = command.replace('python /opt/motif_pipeline/scripts/average_contributions.py', entrypoint)
        if sys.platform == 'darwin':
            command = command.replace('exec > >(tee average.log) 2>&1', '')
        result = subprocess.run(['bash', '-c', command], cwd=fixture.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((fixture.root / 'averaged.h5').is_file())
        self.assertEqual((fixture.root / 'scores.list').read_text().splitlines(), list(map(str, files)))
        self.assertNotIn('gs://', command)


if __name__ == '__main__':
    unittest.main()
