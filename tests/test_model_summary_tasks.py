import csv
import gzip
import json
import shlex
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import WDL
from check_wdl import command_writes, unsupported_functions, workflow_writes
from test_wdl import CloudGeneratedFileStdLib

ROOT = Path(__file__).resolve().parents[1]


class ModelSummaryTaskTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def task(self, name):
        source = ROOT / 'workflows/summarize_variants.wdl'
        self.assertTrue(source.is_file(), 'Shared model summary workflow is missing')
        return next(t for t in WDL.load(str(source)).tasks if t.name == name)

    def run_task(self, name, args, files, directory, success=True):
        task = self.task(name)
        directory.mkdir()
        env = WDL.values_from_json(dict(docker_image='existing-image', **args),
                                   task.available_inputs, task.required_inputs)
        localized = WDL.Value.rewrite_env_paths(env, lambda value: str(files[value.value]))
        command = task.command.eval(localized, CloudGeneratedFileStdLib('1.0', str(directory))).value
        self.assertNotIn('gs://', command)
        # Execute the actual task code with local copies of unchanged image scripts.
        command = command.replace('/opt/chrombpnet/scripts', str(ROOT / 'scripts'))
        command = command.replace('python ', shlex.quote(sys.executable) + ' ')
        if sys.platform == 'darwin':
            command = command.replace('exec > >(tee summary.log) 2>&1', '')
            command = command.replace('exec > >(tee combine.log) 2>&1', '')
            command = command.replace('exec > >(tee groups.log) 2>&1', '')
        result = subprocess.run(['bash', '-c', command], cwd=directory, text=True,
                                capture_output=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertFalse((directory / 'NEVER').exists())
        return result

    def scores(self, name, model, cell, effect, reverse=False):
        path = self.root / name
        fields = ['model_id', 'cell_type', 'chr', 'pos', 'allele1', 'allele2',
                  'variant_id', 'logfc', 'jsd', 'active_allele_quantile', 'logfc.pval']
        variants = [('001', 10, 'A', 'T'), ('NA', 11, 'A', 'C'),
                    ('deletion', 8, 'AAA', 'A'), ('insertion', 20, '-', 'A'),
                    ('other', 10, 'A', 'T')]
        if reverse:
            variants.reverse()
        with path.open('w') as stream:
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
            writer.writerow(fields)
            for ident, pos, ref, alt in variants:
                writer.writerow([model, cell, 'chr2' if ident == 'other' else 'chr1',
                                 pos, ref, alt, ident, effect, 0.1, 0.8, 0.01])
        return path

    def summary(self, prefix='a', cell='CD4', effects=(1, -0.5), peaks=None,
                reverse=False, expected=None, success=True):
        paths = [self.scores(f"{prefix} ' $(touch NEVER) `touch NEVER` {i} scores.tsv", f'{prefix}_fold{i}', cell,
                             effect, reverse and i == 0) for i, effect in enumerate(effects)]
        if peaks is None:
            peaks = self.root / "peaks ' $(touch NEVER) `touch NEVER`.bed.gz"
            with gzip.open(peaks, 'wt') as stream:
                stream.write('# comment\ntrack name=test\nchr1\t9\t10\nchr1\t19\t20\n')
        uris = [f'gs://scores/{prefix}/{i}' for i in range(len(paths))]
        output = self.root / f'summary_{prefix}_{len(list(self.root.glob("summary_*")))}'
        result = self.run_task('SummarizeModel', dict(score_files=uris,
            model_ids=expected or [f'{prefix}_fold{i}' for i in range(len(paths))],
            model_group=prefix, cell_type=cell, peaks='gs://peaks/input'),
            dict(zip(uris, paths), **{'gs://peaks/input': peaks}), output, success)
        return output, result

    def table(self, directory, name):
        with (directory / name).open() as stream:
            return list(csv.DictReader(stream, delimiter='\t'))

    def test_one_model_summary_and_peak_overlap_boundaries(self):
        output, _ = self.summary()
        rows = self.table(output, 'variant_effects.fold_summary.tsv')
        self.assertEqual(len(rows), 5)
        self.assertEqual({r['model_group'] for r in rows}, {'a'})
        self.assertEqual({r['n_folds'] for r in rows}, {'2'})
        self.assertEqual({r['model_ids'] for r in rows}, {'a_fold0,a_fold1'})
        self.assertEqual({r['variant_id']: r['in_peak'] for r in rows},
                         {'001': 'true', 'NA': 'false', 'deletion': 'true',
                          'insertion': 'true', 'other': 'false'})
        self.assertAlmostEqual(float(rows[0]['logfc.mean']), 0.25)
        self.assertNotIn('logfc.pval.mean', rows[0])
        self.assertEqual(len(self.table(output, 'variant_effects.all_models.tsv')), 10)

    def test_wrong_model_assignment_fails_before_summary(self):
        _, result = self.summary(expected=['b_fold0', 'b_fold1'], success=False)
        self.assertIn('model', (result.stdout + result.stderr).lower())

    def test_cell_type_cli_argument_is_safely_quoted(self):
        cell = "CD4 T ' cell $(touch NEVER) `touch NEVER`"
        output, _ = self.summary(cell=cell)
        self.assertEqual({r['cell_type'] for r in self.table(output, 'variant_effects.fold_summary.tsv')}, {cell})

    def test_valid_unsuffixed_model_id_cannot_close_shell_file_block(self):
        model = 'CHROMBPNET_MODEL_IDS'
        scores = self.scores('delimiter_model.tsv', model, 'CD4', 1)
        peaks = self.root / 'delimiter_peaks.bed'
        peaks.write_text('chr1\t9\t10\n')
        output = self.root / 'delimiter_summary'
        self.run_task('SummarizeModel', dict(score_files=['gs://scores/0'],
            model_ids=[model], model_group=model, cell_type='CD4', peaks='gs://peaks/0'),
            {'gs://scores/0': scores, 'gs://peaks/0': peaks}, output)
        self.assertEqual(self.table(output, 'variant_effects.fold_summary.tsv')[0]['model_group'], model)

    def test_nonlocalized_peak_is_rejected(self):
        output = self.root / 'bad_localization'
        paths = [self.scores('input.tsv', 'a_fold0', 'CD4', 1)]
        result = self.run_task('SummarizeModel', dict(score_files=['gs://scores/0'],
            model_ids=['a_fold0'], model_group='a', cell_type='CD4', peaks='gs://peaks/0'),
            {'gs://scores/0': paths[0], 'gs://peaks/0': Path('s3://unresolved/peaks')},
            output, success=False)
        self.assertIn('localization', (result.stdout + result.stderr).lower())

    def test_peak_nested_intervals_and_invalid_coordinates(self):
        peaks = self.root / 'nested.bed'
        peaks.write_text('chr1\t0\t30\nchr1\t1\t2\nchr1\t4\t5\n')
        output, _ = self.summary(peaks=peaks)
        self.assertTrue(all(r['in_peak'] == 'true' for r in
                            self.table(output, 'variant_effects.fold_summary.tsv') if r['chr'] == 'chr1'))
        peaks.write_text('chr1\t10\t9\n')
        _, result = self.summary(peaks=peaks, success=False)
        self.assertIn('peak', (result.stdout + result.stderr).lower())

    def manifest(self, rows):
        path = self.root / "manifest ' $(touch NEVER).tsv"
        with path.open('w') as stream:
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n', quoting=csv.QUOTE_NONE)
            writer.writerow(['model_id', 'cell_type', 'model', 'peaks'])
            writer.writerows(rows)
        return path

    def test_plan_groups_interleaved_manifest_without_serializing_file_inputs(self):
        manifest = self.manifest([['a_fold1', 'CD4', 'gs://m/a1', 'gs://p/a'],
                                  ['b_fold0', 'NK', 'gs://m/b0', 'gs://p/b'],
                                  ['a_fold0', 'CD4', 'gs://m/a0', 'gs://p/a']])
        output = self.root / 'plan'
        self.run_task('PlanModelSummaries', dict(model_manifest='gs://manifest', n_score_files=3),
                      {'gs://manifest': manifest}, output)
        plan = json.loads((output / 'summary_groups.json').read_text())
        self.assertEqual(plan, [dict(model_group='a', cell_type='CD4', indices=[2, 0]),
                                dict(model_group='b', cell_type='NK', indices=[1])])
        self.assertNotIn('gs://', (output / 'summary_groups.json').read_text())
        self.assertIn('gs://p/a', (output / 'models.normalized.tsv').read_text())

    def test_plan_rejects_ambiguous_group_and_wrong_file_count(self):
        cases = [([['a_fold0', 'CD4', 'gs://m/0', 'gs://p/a'],
                   ['a_fold1', 'CD4', 'gs://m/1', 'gs://p/b']], 2, 'peak'),
                 ([['a_fold0', 'CD4', 'gs://m/0', 'gs://p/a'],
                   ['a_fold1', 'NK', 'gs://m/1', 'gs://p/a']], 2, 'cell'),
                 ([['a_fold0', 'CD4', 'gs://m/0', 'gs://p/a'],
                   ['a_fold00', 'CD4', 'gs://m/1', 'gs://p/a']], 2, 'fold'),
                 ([['a_fold0', 'CD4', 'gs://m/0', 'gs://p/a']], 2, 'count')]
        for index, (rows, count, message) in enumerate(cases):
            with self.subTest(message=message):
                manifest = self.manifest(rows)
                result = self.run_task('PlanModelSummaries',
                    dict(model_manifest='gs://manifest', n_score_files=count),
                    {'gs://manifest': manifest}, self.root / f'bad_plan_{index}', success=False)
                self.assertIn(message, (result.stdout + result.stderr).lower())

    def test_group_plan_selects_only_its_localized_scores_and_peak_file(self):
        manifest = self.manifest([['a_fold1', 'CD4', 'gs://m/a1', 'gs://p/a'],
                                  ['b_fold0', 'NK', 'gs://m/b0', 'gs://p/b'],
                                  ['a_fold0', 'CD4', 'gs://m/a0', 'gs://p/a']])
        output = self.root / 'selection_plan'
        self.run_task('PlanModelSummaries', dict(model_manifest='gs://manifest', n_score_files=3),
                      {'gs://manifest': manifest}, output)
        doc = WDL.load(str(ROOT / 'workflows/summarize_variants.wdl'))
        cloud_scores = ['gs://scores/a1', 'gs://scores/b0', 'gs://scores/a0']
        env = WDL.values_from_json(dict(model_manifest='gs://manifest',
            score_files=cloud_scores, docker_image='existing-image'),
            doc.workflow.available_inputs, doc.workflow.required_inputs)
        env = env.bind('PlanModelSummaries.groups', WDL.Value.File(str(output / 'summary_groups.json')))
        env = env.bind('PlanModelSummaries.rows', WDL.Value.File(str(output / 'models.normalized.tsv')))
        stdlib = CloudGeneratedFileStdLib('1.0', str(output))
        for node in doc.workflow.body:
            if isinstance(node, WDL.Tree.Decl) and node.name in ['groups', 'manifest_rows']:
                env = env.bind(node.name, node.expr.eval(env, stdlib).coerce(node.type))
        scatter = next(node for node in doc.workflow.body if isinstance(node, WDL.Tree.Scatter))
        env = env.bind('group', env['groups'].value[0])
        for node in scatter.body:
            if isinstance(node, WDL.Tree.Decl) and node.name in ['group_indices', 'peak_index']:
                env = env.bind(node.name, node.expr.eval(env, stdlib).coerce(node.type))
        inner = next(node for node in scatter.body if isinstance(node, WDL.Tree.Scatter))
        selected = {'group_scores': [], 'group_model_ids': []}
        for index in inner.expr.eval(env, stdlib).value:
            member_env = env.bind('index', index)
            for decl in inner.body:
                value = decl.expr.eval(member_env, stdlib).coerce(decl.type)
                selected[decl.name].append(value)
        self.assertTrue(all(isinstance(value, WDL.Value.File) for value in selected['group_scores']))
        self.assertEqual([value.value for value in selected['group_scores']],
                         ['gs://scores/a0', 'gs://scores/a1'])
        peak_decl = next(node for node in scatter.body if isinstance(node, WDL.Tree.Decl) and node.name == 'group_peaks')
        peak = peak_decl.expr.eval(env, stdlib).coerce(peak_decl.type)
        self.assertIsInstance(peak, WDL.Value.File)
        self.assertEqual(peak.value, 'gs://p/a')
        peaks = self.root / 'selected.bed'
        peaks.write_text('chr1\t9\t10\n')
        files = {'gs://scores/a0': self.scores('selected0.tsv', 'a_fold0', 'CD4', 1),
                 'gs://scores/a1': self.scores('selected1.tsv', 'a_fold1', 'CD4', -0.5),
                 'gs://p/a': peaks}
        # A foreign cell's file is not even available to this task.
        args = dict(score_files=[value.value for value in selected['group_scores']],
            model_ids=[value.value for value in selected['group_model_ids']],
            model_group='a', cell_type='CD4', peaks=peak.value)
        directory = self.root / 'selected_summary'
        self.run_task('SummarizeModel', args, files, directory)
        self.assertEqual({r['model_group'] for r in self.table(directory, 'variant_effects.fold_summary.tsv')}, {'a'})

    def combine(self, directories, success=True):
        args, paths = {}, {}
        for key, name in [('long_files', 'variant_effects.all_models.tsv'),
                          ('wide_files', 'variant_effects.wide.tsv'),
                          ('summary_files', 'variant_effects.fold_summary.tsv')]:
            args[key] = [f'gs://results/{key}/{i}' for i in range(len(directories))]
            paths.update({uri: directory / name for uri, directory in zip(args[key], directories)})
        output = self.root / 'combined'
        result = self.run_task('CombineModelEffects', args, paths, output, success)
        return output, result

    def test_combination_keeps_models_separate_and_handles_variant_order(self):
        a, _ = self.summary()
        peaks = self.root / 'nk.bed'
        peaks.write_text('chr1\t10\t11\n')
        b, _ = self.summary(prefix='b', cell='NK', effects=(2,), peaks=peaks, reverse=True)
        output, _ = self.combine([a, b])
        summary = self.table(output, 'variant_effects.fold_summary.tsv')
        self.assertEqual(len(summary), 10)
        self.assertEqual({r['model_group'] for r in summary}, {'a', 'b'})
        self.assertEqual({r['logfc.mean'] for r in summary if r['model_group'] == 'b'}, {'2.0'})
        self.assertEqual(next(r for r in summary if r['model_group'] == 'b' and r['variant_id'] == '001')['in_peak'], 'false')
        self.assertEqual(len(self.table(output, 'variant_effects.all_models.tsv')), 15)
        wide = self.table(output, 'variant_effects.wide.tsv')
        self.assertEqual(len(wide), 5)
        self.assertEqual(wide[0]['variant_id'], '001')
        self.assertEqual(wide[0]['b_fold0.logfc'], '2')
        self.assertIn('a_fold1.logfc', wide[0])

    def test_combination_rejects_mismatched_variant_alleles(self):
        a, _ = self.summary()
        b, _ = self.summary(prefix='b', cell='NK', effects=(2,))
        path = b / 'variant_effects.wide.tsv'
        path.write_text(path.read_text().replace('chr1\t10\tA\tT', 'chr1\t10\tG\tT'))
        _, result = self.combine([a, b], success=False)
        self.assertIn('identity', (result.stdout + result.stderr).lower())

    def test_combination_rejects_missing_variant(self):
        a, _ = self.summary()
        b, _ = self.summary(prefix='b', cell='NK', effects=(2,))
        path = b / 'variant_effects.wide.tsv'
        lines = path.read_text().splitlines()
        path.write_text('\n'.join(lines[:-1]) + '\n')
        _, result = self.combine([a, b], success=False)
        self.assertIn('variant set', (result.stdout + result.stderr).lower())

    def test_workflow_scatter_selects_typed_files_and_has_no_scoring_calls(self):
        source = ROOT / 'workflows/summarize_variants.wdl'
        self.assertTrue(source.is_file(), 'Shared model summary workflow is missing')
        doc = WDL.load(str(source))
        self.assertEqual(workflow_writes(doc), [])
        self.assertEqual(command_writes(doc), [])
        self.assertEqual(unsupported_functions(doc), [])
        scatter = next(node for node in doc.workflow.body if isinstance(node, WDL.Tree.Scatter))
        summary = next(node for node in scatter.body if isinstance(node, WDL.Tree.Call))
        self.assertEqual(summary.callee.name, 'SummarizeModel')
        self.assertEqual(str(summary.callee.available_inputs['score_files'].type), 'Array[File]')
        self.assertEqual(str(summary.callee.available_inputs['peaks'].type), 'File')
        final = next(node for node in doc.workflow.body if isinstance(node, WDL.Tree.Call)
                     and node.callee.name == 'CombineModelEffects')
        self.assertIn('SummarizeModel', str(final.inputs['summary_files']))


if __name__ == '__main__':
    unittest.main()
