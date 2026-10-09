import json
import re
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

    def local_commands(self, command):
        for script in ['write_file_list.py', 'prepare_motif_inputs.py']:
            command = command.replace('python /opt/motif_pipeline/scripts/' + script,
                shlex.quote(sys.executable) + ' ' + shlex.quote(str(ROOT / 'motif_pipeline/scripts' / script)))
        if sys.platform == 'darwin':
            command = re.sub(r'exec > >\(tee [a-z]+\.log\) 2>&1', '', command)
        return command

    def render_cloud(self, task_name, args, root):
        task = self.task(task_name)
        mapping = {}
        cloud_args = dict(args)
        for binding in task.available_inputs:
            name, typ = binding.name, binding.value.type
            if name not in args or args[name] is None:
                continue
            value = args[name]
            if isinstance(typ, WDL.Type.File):
                cloud = 'gs://test-bucket/' + name
                mapping[cloud] = str(value)
                cloud_args[name] = cloud
            elif isinstance(typ, WDL.Type.Array) and isinstance(typ.item_type, WDL.Type.File):
                clouds = ['gs://test-bucket/' + name + str(i) for i in range(len(value))]
                mapping.update(zip(clouds, map(str, value)))
                cloud_args[name] = clouds
        env = WDL.values_from_json(cloud_args, task.available_inputs, task.required_inputs)
        env = WDL.Value.rewrite_env_paths(env, lambda file: mapping[file.value])
        return self.local_commands(task.command.eval(env, CloudGeneratedFileStdLib('1.0', root)).value)

    def test_preparation_cloud_localization_with_and_without_index(self):
        import h5py
        from pyfaidx import Fasta
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            models = []
            for fold in range(5):
                path = root / f"model ' $(touch NEVER) {fold}.h5"
                with h5py.File(path, 'w') as handle:
                    handle.attrs['model_config'] = json.dumps({'config': {'layers': [
                        {'class_name': 'InputLayer', 'config': {'batch_input_shape': [None, 64, 4]}}]}})
                    handle.create_dataset('weights', data=[fold])
                models.append(path)
            genome = root / "genome ' $(touch NEVER).fa"
            genome.write_text('>chr1\n' + 'ACGT' * 50 + '\n')
            with Fasta(str(genome)):
                pass
            peaks, sizes = root / 'peaks.bed', root / 'sizes.tsv'
            peaks.write_text('chr1\t60\t80\tpeak\n')
            sizes.write_text('chr1\t200\n')
            for use_index in [False, True]:
                with self.subTest(index=use_index):
                    args = dict(fold_models=models, peaks=peaks, genome=genome,
                        genome_index=str(genome) + '.fai' if use_index else None,
                        chrom_sizes=sizes, cell_type='CD4', head='counts', input_length=64,
                        discovery_window=40, random_seed=7, max_peaks=None, docker_image='test')
                    command = self.render_cloud('PrepareMotifInputs', args, root)
                    result = subprocess.run(['bash', '-c', command], cwd=root, capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertTrue((root / 'prepared/reference.fa.fai').is_file())
                    self.assertEqual((root / 'models.list').read_text().splitlines(), [str(path.resolve()) for path in models])
                    self.assertEqual('--genome-index' in command, use_index)
                    self.assertNotIn('gs://', command)
            self.assertFalse((root / 'NEVER').exists())

    def test_all_compute_tasks_receive_readable_localized_file_arguments(self):
        cases = {
            'FoldContributions': (['model', 'peaks', 'genome', 'genome_index', 'preparation'],
                dict(fold_index=0, cell_type="CD4 ' $(touch NEVER)", head='profile', random_seed=7,
                     num_backgrounds=2, batch_size=16, threads=2, memory_gb=64, disk_gb=100, num_preempt=0)),
            'DiscoverMotifs': (['contributions'], dict(window=400, max_seqlets=100, n_leiden=1,
                random_seed=7, cpu_count=2, memory_gb=8, disk_gb=100)),
            'ReportMotifs': (['motifs', 'motif_database'], dict(match_qvalue=0.05, n_matches=5)),
        }
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stub = root / 'capture.py'
            stub.write_text("import json, sys\nfrom pathlib import Path\nargs=sys.argv[1:]\n"
                "for option in " + repr(['--model', '--peaks', '--genome', '--genome-index',
                    '--preparation', '--contributions', '--motifs', '--motif-database']) + ":\n"
                "    if option in args:\n        path=args[args.index(option)+1]\n"
                "        assert '://' not in path and Path(path).is_file(), path\n"
                "Path('captured.json').write_text(json.dumps(args))\n")
            for task_name, (file_names, options) in cases.items():
                with self.subTest(task=task_name):
                    args = dict(options, docker_image='test')
                    for name in file_names:
                        path = root / (name + " ' $(touch NEVER) `touch NEVER` file")
                        path.write_text('input')
                        args[name] = path
                    command = self.render_cloud(task_name, args, root)
                    command = re.sub(r'python /opt/motif_pipeline/scripts/[a-z_]+\.py',
                        lambda match: shlex.quote(sys.executable) + ' ' + shlex.quote(str(stub)), command)
                    command = command.replace('nvidia-smi', 'true')
                    result = subprocess.run(['bash', '-c', command], cwd=root, capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    received = json.loads((root / 'captured.json').read_text())
                    for name in file_names:
                        flag = '--' + name.replace('_', '-')
                        self.assertEqual(received[received.index(flag) + 1], str(args[name]))
                    self.assertNotIn('gs://', command)
            self.assertFalse((root / 'NEVER').exists())

    def test_newline_delimiter_path_cannot_execute_shell_code(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            malicious = root / "bad\nMOTIF_SCORE_FILES\nprintf injected > NEVER\n#"
            malicious.parent.mkdir(parents=True, exist_ok=True)
            malicious.write_text('input')
            files = [malicious] + [root / f'fold{i}.h5' for i in range(4)]
            for path in files[1:]:
                path.write_text('input')
            args = dict(score_files=files, region_files=files, docker_image='test',
                        expected_folds=5, chunk_rows=1, memory_gb=8, disk_gb=100)
            command = self.render_cloud('AverageContributions', args, root)
            result = subprocess.run(['bash', '-c', command], cwd=root, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('invalid path', result.stdout + result.stderr)
            self.assertFalse((root / 'NEVER').exists())
            self.assertFalse((root / 'scores.list').exists())

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
        command = self.local_commands(command)
        if sys.platform == 'darwin':
            command = command.replace('exec > >(tee average.log) 2>&1', '')
        result = subprocess.run(['bash', '-c', command], cwd=fixture.root, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue((fixture.root / 'averaged.h5').is_file())
        self.assertEqual((fixture.root / 'scores.list').read_text().splitlines(), [str(path.resolve()) for path in files])
        self.assertNotIn('gs://', command)


if __name__ == '__main__':
    unittest.main()
