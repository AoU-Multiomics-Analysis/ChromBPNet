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
from test_wdl import LocalStdLib, CloudGeneratedFileStdLib

ROOT = Path(__file__).resolve().parents[1]


class FinemoWDLTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.docs = [WDL.load(str(ROOT / 'workflows' / name)) for name in ['call_motifs.wdl', 'annotate_variant_motifs.wdl']]

    def test_wdl_10_static_checks(self):
        for doc in self.docs:
            self.assertEqual(doc.wdl_version, '1.0')
            self.assertEqual(workflow_writes(doc), [])
            self.assertEqual(command_writes(doc), [])
            self.assertEqual(unsupported_functions(doc), [])

    def test_cloud_localization_for_every_new_task_and_optional_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stub = root / 'capture.py'
            stub.write_text('import sys,json\nfrom pathlib import Path\nargs=sys.argv[1:]\n'
                'flags=' + repr(['--contributions', '--regions', '--motifs', '--motif-provenance', '--coordinates',
                    '--hits', '--qc', '--preparation', '--call-metadata', '--tf-matches', '--model', '--sequences',
                    '--variants', '--genome', '--genome-index', '--chrom-sizes', '--models-list', '--hits-list',
                    '--changes-list', '--variants-list']) + '\n'
                'for flag in flags:\n'
                '    if flag in args:\n'
                '        value=args[args.index(flag)+1]\n'
                '        assert "://" not in value and Path(value).is_file(),value\n'
                'Path("captured.json").write_text(json.dumps(args))\n')
            for doc in self.docs:
                for task in doc.tasks:
                    for present in [False, True]:
                        with self.subTest(task=task.name, optional=present):
                            args, mapping, expected = {}, {}, {}
                            for decl in task.inputs:
                                typ = decl.type
                                if typ.optional and not present:
                                    args[decl.name] = None
                                    continue
                                if isinstance(typ, WDL.Type.File):
                                    path = root / (decl.name + " ' $(touch NEVER) `touch NEVER` file")
                                    path.write_text('input')
                                    cloud = 'gs://bucket/' + decl.name
                                    args[decl.name], mapping[cloud], expected[decl.name] = cloud, str(path), str(path)
                                elif isinstance(typ, WDL.Type.Array):
                                    clouds = []
                                    for i in range(5 if decl.name == 'fold_models' else 2):
                                        path = root / (decl.name + str(i) + " ' $(touch NEVER) file")
                                        path.write_text('input')
                                        cloud = 'gs://bucket/' + decl.name + str(i)
                                        clouds.append(cloud)
                                        mapping[cloud] = str(path)
                                    args[decl.name] = clouds
                                elif isinstance(typ, WDL.Type.String):
                                    args[decl.name] = {'head': 'counts', 'kind': 'variants', 'mode': 'pp'}.get(decl.name, "test ' $(touch NEVER)")
                                elif isinstance(typ, WDL.Type.Int):
                                    args[decl.name] = 5 if decl.name == 'input_length' else 2
                                elif isinstance(typ, WDL.Type.Float):
                                    args[decl.name] = 0.7
                            env = WDL.values_from_json(args, task.available_inputs, task.required_inputs)
                            env = WDL.Value.rewrite_env_paths(env, lambda value: mapping[value.value])
                            for decl in task.postinputs:
                                env = env.bind(decl.name, decl.expr.eval(env, LocalStdLib('1.0')).coerce(decl.type))
                            command = task.command.eval(env, CloudGeneratedFileStdLib('1.0', root)).value
                            self.assertNotIn('gs://', command)
                            command = re.sub(r'python /opt/(finemo_pipeline|motif_pipeline)/scripts/([a-z_]+\.py)',
                                lambda match: shlex.quote(sys.executable) + ' ' + shlex.quote(str(
                                    ROOT / match[1] / 'scripts' / match[2] if match[2] in ['write_local_files.py', 'write_file_list.py'] else stub)), command)
                            command = command.replace('nvidia-smi', 'true')
                            if sys.platform == 'darwin':
                                command = re.sub(r'exec > >\(tee [a-z_]+\.log\) 2>&1', '', command)
                            result = subprocess.run(['bash', '-c', command], cwd=root, capture_output=True, text=True)
                            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                            received = json.loads((root / 'captured.json').read_text())
                            for name, value in expected.items():
                                flag = '--' + name.replace('_', '-')
                                self.assertIn(flag, received)
                                self.assertEqual(received[received.index(flag)+1], value)
                            self.assertNotIn('--config', received)
                            if not present:
                                for decl in task.inputs:
                                    if isinstance(decl.type, WDL.Type.File) and decl.type.optional:
                                        self.assertNotIn('--' + decl.name.replace('_', '-'), received)
            self.assertFalse((root / 'NEVER').exists())

    def test_task_files_remain_files(self):
        names = {'contributions', 'regions', 'motifs', 'motif_provenance', 'coordinates', 'hits', 'qc',
                 'preparation', 'call_metadata', 'tf_matches', 'model', 'sequences', 'variants', 'genome',
                 'genome_index', 'chrom_sizes'}
        for doc in self.docs:
            for task in doc.tasks:
                for decl in task.inputs:
                    if decl.name in names:
                        self.assertIsInstance(decl.type, WDL.Type.File)
                    if decl.name in ['fold_models', 'hit_files', 'change_files', 'variant_files']:
                        self.assertIsInstance(decl.type.item_type, WDL.Type.File)


if __name__ == '__main__':
    unittest.main()
