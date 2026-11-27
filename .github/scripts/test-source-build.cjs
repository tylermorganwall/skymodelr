#!/usr/bin/env node
'use strict';
const assert = require('node:assert/strict');
const crypto = require('node:crypto');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const { spawnSync } = require('node:child_process');
const { test } = require('node:test');

// Test the actual YAML blocks, rather than a second implementation of them.
const workflow = path.resolve(process.env.SOURCE_BUILD_WORKFLOW || '.github/workflows/source-build.yaml');
const yaml = fs.readFileSync(workflow, 'utf8');
function jobBlock(name) {
  const block = yaml.split(`\n  ${name}:\n`)[1];
  assert.ok(block, `Missing job: ${name}`);
  return block.split(/\n  [a-z][a-z0-9_-]*:\n/)[0];
}
function stepBlock(job, name) {
  const block = job.split(`      - name: ${name}\n`)[1];
  assert.ok(block, `Missing step: ${name}`);
  return block.split(/\n      - /)[0];
}
function runBlock(name) {
  const lines = yaml.split('\n');
  const start = lines.findIndex(line => line.trim() === `- name: ${name}`);
  assert.notEqual(start, -1, `Missing step: ${name}`);
  const begin = lines.findIndex((line, i) => i > start && line.trim() === 'run: |');
  assert.notEqual(begin, -1, `Missing run block: ${name}`);
  const end = lines.findIndex((line, i) => i > begin && line.trim() && !line.startsWith('          '));
  return lines.slice(begin + 1, end === -1 ? undefined : end).map(line => line.slice(10)).join('\n');
}
const build = runBlock('Build the triggering commit');
const label = runBlock('Verify and label the source package');
const verify = runBlock('Verify the downloaded tarball and metadata');
const verifyPublisherMetadata = runBlock('Verify downloaded metadata against the build job');
const sha256 = bytes => crypto.createHash('sha256').update(bytes).digest('hex');
const description = `Package: testpkg
Version: 1.2.3-4
Title: Source Build Test Fixture
Description: Provides a minimal package for workflow regression tests.
Authors@R: person("Test", "Runner", email = "test@example.com", role = c("aut", "cre"))
License: MIT
Encoding: UTF-8
`;
function command(program, args, options = {}) {
  const result = spawnSync(program, args, { encoding: 'utf8', ...options });
  if (result.error) throw result.error;
  return result;
}
function passed(result) {
  assert.equal(result.status, 0, result.stdout + result.stderr);
  return result;
}
function fixture(t, branch = 'main') {
  const root = fs.mkdtempSync(path.join(os.tmpdir(), 'source-workflow-test-'));
  t.after(() => fs.rmSync(root, { recursive: true, force: true }));
  const source = path.join(root, 'source');
  const runner = path.join(root, 'runner');
  fs.mkdirSync(path.join(source, 'R'), { recursive: true });
  fs.mkdirSync(path.join(runner, 'source-build'), { recursive: true });
  fs.mkdirSync(path.join(runner, 'source-output'));
  fs.writeFileSync(path.join(source, 'DESCRIPTION'), description);
  fs.writeFileSync(path.join(source, 'NAMESPACE'), 'export(test_value)\n');
  fs.writeFileSync(path.join(source, 'R', 'value.R'), 'test_value <- function() 42L\n');
  const env = { ...process.env, GITHUB_WORKSPACE: source, RUNNER_TEMP: runner,
    GITHUB_REF_NAME: branch, GITHUB_SHA: 'a'.repeat(40), GITHUB_RUN_ID: '123456',
    GITHUB_RUN_NUMBER: '987', GITHUB_RUN_ATTEMPT: '2', GITHUB_OUTPUT: path.join(root, 'github-output') };
  const f = { root, source, runner, env, output: path.join(runner, 'source-output'),
    download: path.join(source, 'source-output'), archive: path.join(runner, 'source-build', 'testpkg_1.2.3-4.tar.gz') };
  f.run = block => command('bash', ['--noprofile', '--norc', '-c', block], { env, cwd: source });
  f.tar = (contents = description, member = 'testpkg/DESCRIPTION') => {
    const staging = path.join(root, 'staging');
    fs.rmSync(staging, { recursive: true, force: true });
    fs.mkdirSync(path.dirname(path.join(staging, member)), { recursive: true });
    fs.writeFileSync(path.join(staging, member), contents);
    passed(command('tar', ['-czf', f.archive, '-C', staging, member]));
  };
  f.label = () => {
    passed(f.run(label));
    const output = fs.readFileSync(env.GITHUB_OUTPUT, 'utf8');
    assert.match(output, /^metadata-sha256=[0-9a-f]{64}\n$/);
    env.EXPECTED_METADATA_SHA256 = output.trim().split('=')[1];
    assert.equal(env.EXPECTED_METADATA_SHA256, sha256(fs.readFileSync(path.join(f.output, 'source-build.json'))));
  };
  f.metadata = () => JSON.parse(fs.readFileSync(path.join(f.output, 'source-build.json'), 'utf8'));
  f.downloadArtifact = () => fs.cpSync(f.output, f.download, { recursive: true });
  return f;
}

for (const [name, block] of [['build', build], ['label', label], ['verify', verify], ['publisher metadata', verifyPublisherMetadata]]) {
  test(`Bash syntax: ${name}`, () => passed(command('bash', ['-n'], { input: block })));
}
test('no Python helper calls and read-only build', () => {
  assert.doesNotMatch(yaml, /build-source\.py/);
  assert.match(yaml, /permissions:\n  contents: read/);
});
test('the build exports the upload ID and metadata digest without downloading', () => {
  const job = jobBlock('build');
  assert.match(job, /artifact-id: \$\{\{ steps\.artifact\.outputs\.artifact-id \}\}/);
  assert.match(job, /metadata-sha256: \$\{\{ steps\.source\.outputs\.metadata-sha256 \}\}/);
  assert.match(stepBlock(job, 'Verify and label the source package'), /id: source/);
  assert.match(stepBlock(job, 'Save verified source package'), /id: artifact/);
  assert.doesNotMatch(job, /download-artifact/);
  assert.doesNotMatch(job, /permissions:/);
});
test('every branch verifies the artifact before main-only publication', () => {
  const verification = jobBlock('verify');
  assert.match(verification, /^    needs: build$/m);
  assert.doesNotMatch(verification, /^    (if|permissions|continue-on-error):/m);
  const publication = jobBlock('publish');
  assert.match(publication, /^    needs: \[build, verify\]$/m);
  assert.match(publication, /if: >-\n      github\.ref == 'refs\/heads\/main' &&\n      \(github\.event_name == 'push' \|\| github\.event_name == 'workflow_dispatch'\)/);
  assert.match(publication, /permissions:\n      contents: write/);
  assert.match(publication, /concurrency:\n      group: source-release-main\n      cancel-in-progress: false/);
  assert.match(publication, /ref: \$\{\{ github\.sha \}\}\n          persist-credentials: false\n          sparse-checkout: \.github\/scripts/);
  assert.doesNotMatch(publication, /setup-r|R CMD|Rscript/);
});
for (const [jobName, verificationName] of [
  ['verify', 'Verify the downloaded tarball and metadata'],
  ['publish', 'Verify downloaded metadata against the build job'],
]) {
  test(`${jobName} downloads the exact build artifact with one failure-only retry before verification`, () => {
    const job = jobBlock(jobName);
    const download = stepBlock(job, 'Download the exact source artifact');
    const wait = stepBlock(job, 'Wait for artifact visibility');
    const retry = stepBlock(job, 'Retry the source artifact download');
    const verification = stepBlock(job, verificationName);
    assert.equal((job.match(/uses: actions\/download-artifact@v4/g) || []).length, 2);
    assert.match(download, /id: download\n        continue-on-error: true/);
    assert.doesNotMatch(download, /\n        if:/);
    for (const step of [download, retry]) {
      assert.match(step, /artifact-ids: \$\{\{ needs\.build\.outputs\.artifact-id \}\}\n          merge-multiple: true\n          path: source-output/);
      assert.doesNotMatch(step, /\n          (name|pattern|run-id):/);
    }
    for (const step of [wait, retry]) assert.match(step, /if: steps\.download\.outcome == 'failure'/);
    assert.match(wait, /run: sleep 20/);
    assert.doesNotMatch(retry + verification, /continue-on-error:/);
    assert.doesNotMatch(verification, /\n        if:/);
    assert.match(verification, /EXPECTED_METADATA_SHA256: \$\{\{ needs\.build\.outputs\.metadata-sha256 \}\}/);
    assert.ok(job.indexOf(download) < job.indexOf(wait));
    assert.ok(job.indexOf(wait) < job.indexOf(retry));
    assert.ok(job.indexOf(retry) < job.indexOf(verification));
    if (jobName === 'publish') assert.ok(job.indexOf(verification) < job.indexOf('uses: actions/github-script@v7'));
  });
}
test('real R CMD build, complete metadata and downloaded integrity', t => {
  const f = fixture(t);
  passed(f.run(build));
  f.label();
  const m = f.metadata();
  const filename = `testpkg_1.2.3-4_main-${sha256('main').slice(0, 12)}_${'a'.repeat(40)}_123456-2.tar.gz`;
  const bytes = fs.readFileSync(path.join(f.output, filename));
  assert.deepEqual(m, { package: 'testpkg', version: '1.2.3-4', branch: 'main', sha: 'a'.repeat(40), run_id: 123456, run_number: 987, run_attempt: 2, filename, sha256: sha256(bytes) });
  assert.equal(fs.readFileSync(path.join(f.output, filename + '.sha256'), 'utf8'), `${m.sha256}  ${filename}\n`);
  f.downloadArtifact();
  // Downstream jobs have only the artifact and the build's independent digest.
  fs.rmSync(f.output, { recursive: true });
  passed(f.run(verify));
  passed(f.run(verifyPublisherMetadata));
});
for (const branch of ['feature/topic', 'feature-topic', '42/merge', 'feat/\u00e9\u2603"$(touch pwned)', 'feature/' + 'a'.repeat(100)]) {
  test(`safe deterministic filename and JSON for branch ${JSON.stringify(branch)}`, t => {
    const f = fixture(t, branch);
    f.tar(); f.label();
    const m = f.metadata();
    const slug = Array.from(branch, ch => /^[A-Za-z0-9._-]$/.test(ch) ? ch : '-').join('').slice(0, 60);
    assert.equal(m.branch, branch);
    assert.equal(m.filename, `testpkg_1.2.3-4_${slug}-${sha256(branch).slice(0, 12)}_${'a'.repeat(40)}_123456-2.tar.gz`);
    assert.match(m.filename, /^[A-Za-z0-9._+-]+\.tar\.gz$/);
    assert.equal(fs.existsSync(path.join(f.source, 'pwned')), false);
    f.downloadArtifact(); passed(f.run(verify));
  });
}
for (const [name, mutate] of [
  ['no tarball', f => fs.unlinkSync(f.archive)],
  ['multiple tarballs', f => fs.copyFileSync(f.archive, f.archive + '.extra.tar.gz')],
  ['hidden extra tarball', f => fs.copyFileSync(f.archive, path.join(path.dirname(f.archive), '.extra.tar.gz'))],
  ['mismatched version', f => f.tar(description.replace('1.2.3-4', '9.0'))],
  ['mismatched package', f => f.tar(description.replace('testpkg', 'otherpkg'))],
  ['missing archive DESCRIPTION', f => f.tar(description, 'testpkg/OTHER')],
  ['multiple archived DESCRIPTION records', f => f.tar(description + '\n' + description)],
  ['missing source version', f => fs.writeFileSync(path.join(f.source, 'DESCRIPTION'), description.replace('Version: 1.2.3-4\n', ''))],
  ['unsafe package', f => fs.writeFileSync(path.join(f.source, 'DESCRIPTION'), description.replace('Package: testpkg', 'Package: ../testpkg'))],
  ['unsafe version', f => fs.writeFileSync(path.join(f.source, 'DESCRIPTION'), description.replace('Version: 1.2.3-4', 'Version: ../1'))],
  ['invalid commit SHA', f => { f.env.GITHUB_SHA = '../oops'; }],
  ['invalid run ID', f => { f.env.GITHUB_RUN_ID = '{"oops":1}'; }],
  ['zero run number', f => { f.env.GITHUB_RUN_NUMBER = '0'; }],
  ['invalid attempt', f => { f.env.GITHUB_RUN_ATTEMPT = 'no'; }],
]) {
  test(`rejects ${name}`, t => {
    const f = fixture(t); f.tar(); mutate(f);
    const result = f.run(label);
    assert.notEqual(result.status, 0, `Unexpected success: ${result.stdout}`);
  });
}

for (const [name, block] of [['artifact verification', verify], ['publisher metadata verification', verifyPublisherMetadata]]) {
  for (const expected of [undefined, '', 'a'.repeat(63), 'G'.repeat(64), 'b'.repeat(64), 'a'.repeat(64) + '\nsource-build.json']) {
    test(`${name} rejects a missing, malformed or incorrect expected digest: ${JSON.stringify(expected)}`, t => {
      const f = fixture(t); f.tar(); f.label(); f.downloadArtifact();
      if (expected === undefined) delete f.env.EXPECTED_METADATA_SHA256;
      else f.env.EXPECTED_METADATA_SHA256 = expected;
      assert.notEqual(f.run(block).status, 0);
    });
  }
  test(`${name} rejects coordinated tampering against the original build digest`, t => {
    const f = fixture(t); f.tar(); f.label(); f.downloadArtifact();
    const m = f.metadata();
    const bytes = Buffer.from('coordinated forgery');
    const digest = sha256(bytes);
    fs.writeFileSync(path.join(f.download, m.filename), bytes);
    fs.writeFileSync(path.join(f.download, m.filename + '.sha256'), `${digest}  ${m.filename}\n`);
    fs.writeFileSync(path.join(f.download, 'source-build.json'), JSON.stringify({ ...m, sha256: digest }));
    fs.rmSync(f.output, { recursive: true });
    assert.notEqual(f.run(block).status, 0);
  });
}
for (const [name, mutate] of [
  ['changed tarball', (f, m) => fs.appendFileSync(path.join(f.download, m.filename), 'tampered')],
  ['changed metadata', f => fs.appendFileSync(path.join(f.download, 'source-build.json'), '\n')],
  ['changed checksum', (f, m) => fs.writeFileSync(path.join(f.download, m.filename + '.sha256'), `${'b'.repeat(64)}  ${m.filename}\n`)],
  ['renamed tarball', (f, m) => fs.renameSync(path.join(f.download, m.filename), path.join(f.download, 'renamed.tar.gz'))],
  ['missing tarball', (f, m) => fs.unlinkSync(path.join(f.download, m.filename))],
  ['extra tarball', (f, m) => fs.copyFileSync(path.join(f.download, m.filename), path.join(f.download, 'extra.tar.gz'))],
  ['hidden extra tarball', (f, m) => fs.copyFileSync(path.join(f.download, m.filename), path.join(f.download, '.extra.tar.gz'))],
  ['missing checksum', (f, m) => fs.unlinkSync(path.join(f.download, m.filename + '.sha256'))],
  ['extra checksum', (f, m) => fs.copyFileSync(path.join(f.download, m.filename + '.sha256'), path.join(f.download, 'extra.sha256'))],
  ['hidden extra checksum', (f, m) => fs.copyFileSync(path.join(f.download, m.filename + '.sha256'), path.join(f.download, '.extra.sha256'))],
  ['extra file', f => fs.writeFileSync(path.join(f.download, 'unexpected.txt'), 'extra')],
  ['extra directory', f => fs.mkdirSync(path.join(f.download, 'unexpected'))],
  ['checksum with wrong filename', (f, m) => fs.writeFileSync(path.join(f.download, m.filename + '.sha256'), `${m.sha256}  renamed.tar.gz\n`)],
  ['checksum with extra line', (f, m) => fs.appendFileSync(path.join(f.download, m.filename + '.sha256'), `${m.sha256}  ${m.filename}\n`)],
  ['missing metadata', f => fs.unlinkSync(path.join(f.download, 'source-build.json'))],
  ['forged tarball and checksum', (f, m) => {
    const bytes = Buffer.from('forged'); const digest = sha256(bytes);
    fs.writeFileSync(path.join(f.download, m.filename), bytes);
    fs.writeFileSync(path.join(f.download, m.filename + '.sha256'), `${digest}  ${m.filename}\n`);
  }],
  ['forged tarball, metadata and checksum', (f, m) => {
    const bytes = Buffer.from('forged'); const digest = sha256(bytes);
    fs.writeFileSync(path.join(f.download, m.filename), bytes);
    fs.writeFileSync(path.join(f.download, m.filename + '.sha256'), `${digest}  ${m.filename}\n`);
    fs.writeFileSync(path.join(f.download, 'source-build.json'), JSON.stringify({ ...m, sha256: digest }, null, 2) + '\n');
  }],
]) {
  test(`rejects downloaded artifact with ${name}`, t => {
    const f = fixture(t); f.tar(); f.label(); f.downloadArtifact(); mutate(f, f.metadata());
    const result = f.run(verify);
    assert.notEqual(result.status, 0, `Unexpected success: ${result.stdout}`);
  });
}
