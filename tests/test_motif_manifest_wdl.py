import csv
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

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / 'workflows/discover_motifs_manifest.wdl'
INPUT_FIELDS = ['model_id', 'cell_type', 'model', 'peaks']
RESULT_FIELDS = ['model_group', 'cell_type', 'head', 'modisco_motifs',
                 'averaged_contributions', 'interpreted_regions', 'report_bundle',
                 'discovered_motifs_meme', 'candidate_tf_matches', 'motif_inventory',
                 'tomtom_results']


class MotifManifestWDLTests(unittest.TestCase):
    def document(self):
        self.assertTrue(WORKFLOW.is_file(), 'Multi-model motif manifest workflow is missing')
        return WDL.load(str(WORKFLOW))

    def rows(self):
        # The two model groups share a cell type. Grouping by cell type is a bug.
        return [dict(model_id=f'{group}_fold{fold}', cell_type="GM12878 ' $(touch NEVER)",
                     model=f"gs://models/{group}/fold{fold}/model ' $(touch NEVER).h5",
                     peaks=f'gs://peaks/{group}.bed')
                for group, fold in [('A', 2), ('B', 4), ('A', 0), ('B', 0), ('A', 4),
                                    ('B', 1), ('A', 1), ('B', 3), ('A', 3), ('B', 2)]]

    def manifest(self, root, rows=None):
        path = root / "models ' $(touch NEVER).tsv"
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, INPUT_FIELDS, delimiter='\t', quoting=csv.QUOTE_NONE, quotechar=None)
            writer.writeheader()
            writer.writerows(self.rows() if rows is None else rows)
        return path

    def render(self, task_name, args, root):
        task = next(task for task in self.document().tasks if task.name == task_name)
        cloud_args, mapping = dict(args), {}
        for binding in task.available_inputs:
            if isinstance(binding.value.type, WDL.Type.File) and binding.name in args:
                uri = 'gs://inputs/' + binding.name
                mapping[uri] = str(args[binding.name])
                cloud_args[binding.name] = uri
        env = WDL.values_from_json(cloud_args, task.available_inputs, task.required_inputs)
        env = WDL.Value.rewrite_env_paths(env, lambda file: mapping[file.value])
        for decl in task.postinputs:
            env = env.bind(decl.name, decl.expr.eval(env, LocalStdLib('1.0')).coerce(decl.type))
        command = task.command.eval(env, CloudGeneratedFileStdLib('1.0', root)).value
        command = command.replace('python - ', shlex.quote(sys.executable) + ' - ')
        if sys.platform == 'darwin':
            command = re.sub(r'exec > >\(tee [a-z_]+\.log\) 2>&1', '', command)
        return command

    def run_task(self, task_name, args, root):
        command = self.render(task_name, args, root)
        return subprocess.run(['bash', '-c', command], cwd=root, text=True, capture_output=True)

    def plan(self, root, rows=None, head='counts'):
        manifest = self.manifest(root, rows)
        result = self.run_task('PlanMotifGroups', dict(model_manifest=manifest,
                               head=head, docker_image='existing-image'), root)
        return manifest, result

    def test_interleaved_folds_are_grouped_by_model_and_sorted(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, result = self.plan(root)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(json.loads((root / 'motif_groups.json').read_text()), [
                dict(model_group='A', cell_type="GM12878 ' $(touch NEVER)", indices=[2, 6, 0, 8, 4]),
                dict(model_group='B', cell_type="GM12878 ' $(touch NEVER)", indices=[3, 5, 9, 7, 1]),
            ])
            self.assertEqual((root / 'models.normalized.tsv').read_text().splitlines()[0],
                             "A_fold2\tGM12878 ' $(touch NEVER)\tgs://models/A/fold2/model ' $(touch NEVER).h5\tgs://peaks/A.bed")
            self.assertFalse((root / 'NEVER').exists())

    def test_invalid_groups_fail_before_any_model_is_opened(self):
        cases = ['missing', 'duplicate', 'extra_fold', 'unsuffixed', 'cell_type', 'peaks', 'model', 'bad_head']
        for case in cases:
            with self.subTest(case=case), tempfile.TemporaryDirectory() as tmp:
                rows = self.rows()
                head = 'counts'
                if case == 'missing':
                    rows.pop()
                elif case == 'duplicate':
                    rows.append(dict(rows[0]))
                elif case == 'extra_fold':
                    rows[0]['model_id'] = 'A_fold5'
                elif case == 'unsuffixed':
                    rows[0]['model_id'] = 'A'
                elif case == 'cell_type':
                    rows[0]['cell_type'] = 'NK'
                elif case == 'peaks':
                    rows[0]['peaks'] = 'gs://other/peaks.bed'
                elif case == 'model':
                    rows[0]['model'] = rows[2]['model']
                elif case == 'bad_head':
                    head = 'both'
                _, result = self.plan(Path(tmp), rows, head)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertFalse((Path(tmp) / 'motif_groups.json').exists())

    def test_workflow_selects_localizable_files_for_each_group(self):
        doc = self.document()
        scatter = next(node for node in doc.workflow.body if isinstance(node, WDL.Tree.Scatter))
        inner = next(node for node in scatter.body if isinstance(node, WDL.Tree.Scatter))
        model_decl = next(node for node in inner.body if isinstance(node, WDL.Tree.Decl))
        peaks_decl = next(node for node in scatter.body if isinstance(node, WDL.Tree.Decl) and node.name == 'group_peaks')
        group_type = WDL.Type.StructInstance('MotifGroup')
        group_type.members = doc.struct_typedefs['MotifGroup'].members
        rows = WDL.Value.Array(WDL.Type.Array(WDL.Type.String()), [
            WDL.Value.Array(WDL.Type.String(), [WDL.Value.String(row[name]) for name in INPUT_FIELDS])
            for row in self.rows()])
        group = WDL.Value.Struct(group_type, dict(model_group=WDL.Value.String('A'),
                 cell_type=WDL.Value.String('GM12878'), indices=WDL.Value.Array(WDL.Type.Int(),
                    [WDL.Value.Int(i) for i in [2, 6, 0, 8, 4]])))
        env = WDL.Env.Bindings().bind('manifest_rows', rows).bind('group', group)
        model = model_decl.expr.eval(env.bind('index', WDL.Value.Int(2)), LocalStdLib('1.0')).coerce(model_decl.type)
        peaks = peaks_decl.expr.eval(env, LocalStdLib('1.0')).coerce(peaks_decl.type)
        self.assertIsInstance(model, WDL.Value.File)
        self.assertIsInstance(peaks, WDL.Value.File)
        self.assertEqual(model.value, "gs://models/A/fold0/model ' $(touch NEVER).h5")
        self.assertEqual(peaks.value, 'gs://peaks/A.bed')

    def result_argument(self, group, head='counts', bad_uri=False, cell_type=None, double_quotes=False,
                        uri_suffix=" ' $(touch NEVER) `touch NEVER`"):
        # Evaluate the real workflow metadata declarations, including shell quoting.
        doc = self.document()
        scatter = next(node for node in doc.workflow.body if isinstance(node, WDL.Tree.Scatter))
        group_type = WDL.Type.StructInstance('MotifGroup')
        group_type.members = doc.struct_typedefs['MotifGroup'].members
        value = WDL.Value.Struct(group_type, dict(model_group=WDL.Value.String(group),
                  cell_type=WDL.Value.String(cell_type or "GM12878 ' $(touch NEVER)"),
                  indices=WDL.Value.Array(WDL.Type.Int(), [WDL.Value.Int(0)])))
        env = WDL.Env.Bindings().bind('group', value).bind('head', WDL.Value.String(head))
        for name in RESULT_FIELDS[3:]:
            uri = f'gs://outputs/{group}/{name}' + uri_suffix
            if double_quotes:
                uri += ' "quoted"'
            if bad_uri and name == 'modisco_motifs':
                uri += '\nINJECTED\nprintf bad > NEVER\n'
            env = env.bind('ModelMotifs.' + name, WDL.Value.File(uri))
        for decl in scatter.body:
            if isinstance(decl, WDL.Tree.Decl) and decl.name != 'group_peaks':
                env = env.bind(decl.name, decl.expr.eval(env, LocalStdLib('1.0')).coerce(decl.type))
        return env['result_row_arg'].value

    def test_output_manifests_keep_cloud_uris_and_original_fold_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest, planned = self.plan(root)
            self.assertEqual(planned.returncode, 0, planned.stderr)
            args = dict(model_manifest=manifest, group_plan=root / 'motif_groups.json', head='counts',
                        result_row_args=[self.result_argument('B'), self.result_argument('A')],
                        docker_image='existing-image')
            result = self.run_task('WriteMotifManifests', args, root)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (root / 'modisco_manifest.tsv').open() as stream:
                reader = csv.DictReader(stream, delimiter='\t', quoting=csv.QUOTE_NONE)
                self.assertEqual(reader.fieldnames, RESULT_FIELDS)
                results = list(reader)
            self.assertEqual([row['model_group'] for row in results], ['A', 'B'])
            self.assertEqual(results[0]['modisco_motifs'], "gs://outputs/A/modisco_motifs ' $(touch NEVER) `touch NEVER`")
            with (root / 'models.with_modisco.tsv').open() as stream:
                enriched = list(csv.DictReader(stream, delimiter='\t', quoting=csv.QUOTE_NONE))
            self.assertEqual([{name: row[name] for name in INPUT_FIELDS} for row in enriched], self.rows())
            self.assertEqual(enriched[0]['modisco_motifs'], results[0]['modisco_motifs'])
            self.assertEqual(enriched[1]['modisco_motifs'], results[1]['modisco_motifs'])
            self.assertFalse((root / 'NEVER').exists())

    def test_double_quotes_are_preserved_in_both_manifests(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = self.rows()
            for row in rows:
                row['cell_type'] = 'GM12878 "reference"'
                row['model'] += ' "model"'
                row['peaks'] += ' "peaks"'
            manifest, planned = self.plan(root, rows)
            self.assertEqual(planned.returncode, 0, planned.stdout + planned.stderr)
            arguments = [self.result_argument(group, cell_type='GM12878 "reference"', double_quotes=True)
                         for group in ['A', 'B']]
            result = self.run_task('WriteMotifManifests', dict(model_manifest=manifest,
                group_plan=root / 'motif_groups.json', head='counts', result_row_args=arguments,
                docker_image='existing-image'), root)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (root / 'models.with_modisco.tsv').open() as stream:
                enriched = list(csv.DictReader(stream, delimiter='\t', quoting=csv.QUOTE_NONE))
            self.assertEqual([{field: row[field] for field in INPUT_FIELDS} for row in enriched], rows)
            with (root / 'modisco_manifest.tsv').open() as stream:
                results = list(csv.DictReader(stream, delimiter='\t', quoting=csv.QUOTE_NONE))
            self.assertEqual(results[0]['cell_type'], 'GM12878 "reference"')
            self.assertEqual(results[0]['modisco_motifs'],
                             'gs://outputs/A/modisco_motifs \' $(touch NEVER) `touch NEVER` "quoted"')

    def test_existing_model_ids_with_leading_hyphens_or_underscores_are_accepted(self):
        sys.path.insert(0, str(ROOT / 'scripts'))
        from validate_manifest import read_manifest
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = self.rows()
            for row in rows:
                row['model_id'] = ('-' if row['model_id'].startswith('A') else '_') + row['model_id']
                row['cell_type'] = 'GM12878'
            manifest = self.manifest(root, rows)
            self.assertEqual(len(read_manifest(manifest)), 10)
            manifest, planned = self.plan(root, rows)
            self.assertEqual(planned.returncode, 0, planned.stdout + planned.stderr)
            result = self.run_task('WriteMotifManifests', dict(model_manifest=manifest,
                group_plan=root / 'motif_groups.json', head='counts',
                result_row_args=[self.result_argument(group, cell_type='GM12878', uri_suffix='')
                                 for group in ['-A', '_B']],
                docker_image='existing-image'), root)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with (root / 'modisco_manifest.tsv').open() as stream:
                results = list(csv.DictReader(stream, delimiter='\t', quoting=csv.QUOTE_NONE))
            self.assertEqual([row['model_group'] for row in results], ['-A', '_B'])

    def test_output_writer_rejects_incomplete_duplicate_or_invalid_results(self):
        for case in ['missing', 'duplicate', 'wrong_group', 'wrong_head', 'newline']:
            with self.subTest(case=case), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                manifest, planned = self.plan(root)
                self.assertEqual(planned.returncode, 0, planned.stderr)
                rows = [self.result_argument('A'), self.result_argument('B')]
                if case == 'missing':
                    rows.pop()
                elif case == 'duplicate':
                    rows[1] = rows[0]
                elif case == 'wrong_group':
                    rows[1] = self.result_argument('C')
                elif case == 'wrong_head':
                    rows[1] = self.result_argument('B', head='profile')
                elif case == 'newline':
                    rows[0] = self.result_argument('A', bad_uri=True)
                result = self.run_task('WriteMotifManifests', dict(model_manifest=manifest,
                    group_plan=root / 'motif_groups.json', head='counts', result_row_args=rows,
                    docker_image='existing-image'), root)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse((root / 'modisco_manifest.tsv').exists())
                self.assertFalse((root / 'models.with_modisco.tsv').exists())
                self.assertFalse((root / 'NEVER').exists())

    def test_terra_wdl_10_and_no_engine_file_writes(self):
        doc = self.document()
        self.assertEqual(doc.wdl_version, '1.0')
        self.assertEqual(workflow_writes(doc), [])
        self.assertEqual(command_writes(doc), [])
        self.assertEqual(unsupported_functions(doc), [])
        task = next(task for task in doc.tasks if task.name == 'WriteMotifManifests')
        inputs = {decl.name: decl.type for decl in task.inputs}
        self.assertIsInstance(inputs['model_manifest'], WDL.Type.File)
        self.assertIsInstance(inputs['group_plan'], WDL.Type.File)
        self.assertIsInstance(inputs['result_row_args'].item_type, WDL.Type.String)


if __name__ == '__main__':
    unittest.main()
