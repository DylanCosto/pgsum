#!/usr/bin/env python3
"""Build or verify deterministic, self-contained pgsum binary release archives.

Python 3.8+; standard library only. Nothing in this script publishes a release.
Identical payload bytes, version, target, SOURCE_DATE_EPOCH, and packaging/compression
runtime produce identical archives. This does not assert reproducible compilation.
"""

import argparse
import gzip
import hashlib
import html
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import tarfile
import tempfile
from urllib.parse import unquote, urlsplit


SCHEMA = 'pgsum-release-package-v1'
MANIFEST = 'package-manifest.json'
TOP_FILES = ['README.md', 'DESIGN.md', 'RELEASE_NOTES.md', 'LICENSE-MIT', 'LICENSE-APACHE']
BENCH_FILES = ['README.md', 'performance.py', 'missingness.py', 'compare_plink.py',
               'public_chr22.py', 'public_chr22.requirements.txt', 'robustness_public.py',
               'robustness_genomewide.py', 'robustness_genomewide.requirements.txt',
               'robustness_scale.py', 'test_robustness_scale.py', 'check_study_summary.py', 'test_study_summary.py']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def version_at(root):
    cargo = (root / 'Cargo.toml').read_text(encoding='utf-8')
    section = re.search(r'^\[package\]\s*\n(.*?)(?=^\[|\Z)', cargo, re.M | re.S)
    require(section is not None, 'Cargo.toml has no [package] section')
    match = re.search(r'^version\s*=\s*"([0-9]+\.[0-9]+\.[0-9]+(?:-[A-Za-z0-9.-]+)?)"', section.group(1), re.M)
    require(match is not None, 'could not read package version')
    return match.group(1)


def safe_path(name):
    path = PurePosixPath(name)
    require(name and not path.is_absolute() and '\\' not in name,
            'package member must be a portable relative path: ' + name)
    require(all(part not in ('', '.', '..') for part in name.split('/')),
            'unsafe package member: ' + name)
    return path


def without_fences(text):
    result, inside = [], None
    for line in text.splitlines():
        stripped = line.lstrip()
        if stripped.startswith(('```', '~~~')):
            marker = stripped[:3]
            if inside is None:
                inside = marker
            elif inside == marker:
                inside = None
            continue
        if inside is None:
            result.append(line)
    return '\n'.join(result)


def markdown_links(text):
    text = without_fences(text)
    inline = re.findall(r'!?\[[^\]]*\]\(\s*(<[^>]+>|[^\s)]+)', text)
    defined = re.findall(r'^\s*\[[^\]]+\]:\s*(<[^>]+>|\S+)', text, re.M)
    return [html.unescape(value.strip('<>')) for value in inline + defined]


def anchors(text):
    used, result = {}, set()
    for line in without_fences(text).splitlines():
        match = re.match(r'^\s{0,3}#{1,6}\s+(.+?)\s*#*\s*$', line)
        if match is None:
            continue
        heading = re.sub(r'<[^>]+>', '', html.unescape(match.group(1))).lower()
        slug = ''.join(char for char in heading if char.isalnum() or char in '_- ')
        slug = slug.replace(' ', '-')
        count = used.get(slug, 0)
        used[slug] = count + 1
        result.add(slug if count == 0 else '%s-%d' % (slug, count))
    return result


def check_markdown(payload):
    checked = 0
    for name, data in payload.items():
        if not name.endswith('.md'):
            continue
        text = data.decode('utf-8')
        for link in markdown_links(text):
            url = urlsplit(link)
            if url.scheme or url.netloc:
                continue
            require(not url.path.startswith('/'), 'nonportable absolute Markdown link in %s: %s' % (name, link))
            parts = list(PurePosixPath(name).parent.parts)
            for part in unquote(url.path).split('/') if url.path else []:
                if part in ('', '.'):
                    continue
                if part == '..':
                    require(parts, 'Markdown link escapes archive in %s: %s' % (name, link))
                    parts.pop()
                else:
                    parts.append(part)
            target = '/'.join(parts) if url.path else name
            exists = target in payload or any(item.startswith(target.rstrip('/') + '/') for item in payload)
            require(exists, 'missing packaged Markdown target in %s: %s' % (name, link))
            if url.fragment and target.endswith('.md'):
                require(unquote(url.fragment) in anchors(payload[target].decode('utf-8')),
                        'missing packaged Markdown anchor in %s: %s' % (name, link))
            checked += 1
    return checked


def collect(root, binary):
    payload = {}

    def add(path, name=None):
        require(path.is_file() and not path.is_symlink(), 'expected regular input file: ' + str(path))
        name = name or path.relative_to(root).as_posix()
        safe_path(name)
        payload[name] = path.read_bytes()

    add(binary, 'pgsum')
    for name in TOP_FILES:
        add(root / name)
    add(root / 'scripts/package_release.py')
    for directory in ['docs', 'examples/robustness']:
        for path in sorted((root / directory).rglob('*')):
            if path.is_file():
                require('__pycache__' not in path.parts and path.suffix != '.pyc',
                        'remove generated Python cache from package input: ' + str(path))
                add(path)
    for name in BENCH_FILES:
        add(root / 'bench' / name)
    # These small checked records support the links in the shipped documentation;
    # upstream VCFs and scoring files are deliberately not release payloads.
    for path in sorted((root / 'bench/results').glob('*.json')):
        add(path)
    for name in ['Cargo.toml', 'Cargo.lock', 'src/main.rs']:
        add(root / 'bench/plink2-export' / name)
    require('examples/robustness/run.py' in payload, 'offline example runner is missing')
    check_markdown(payload)
    return payload


def mode_for(name):
    return 0o755 if name == 'pgsum' or name.endswith('.py') else 0o644


def archive_info(name, mode, epoch, size=0, directory=False):
    entry = tarfile.TarInfo(name)
    entry.uid = entry.gid = 0
    entry.uname = entry.gname = ''
    entry.mode, entry.mtime, entry.size = mode, epoch, size
    if directory:
        entry.type = tarfile.DIRTYPE
    return entry


def build(args):
    root = args.source_root.resolve()
    binary = args.binary.resolve()
    version = version_at(root)
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', args.target) is not None, 'invalid target triple')
    epoch = int(args.source_date_epoch)
    require(0 <= epoch <= 0xffffffff, 'SOURCE_DATE_EPOCH must fit an unsigned 32-bit timestamp')
    payload = collect(root, binary)
    archive_root = 'pgsum-v%s-%s' % (version, args.target)
    manifest = {
        'schema': SCHEMA, 'version': version, 'target': args.target,
        'archive_root': archive_root, 'source_date_epoch': epoch,
        'files': {name: {'sha256': sha(data), 'bytes': len(data), 'mode': '%04o' % mode_for(name)}
                  for name, data in sorted(payload.items())},
        'hash_scope': 'Every payload file except package-manifest.json itself; the adjacent archive checksum covers the manifest too.',
        'reproducibility': 'Identical payload bytes, version, target, epoch, and packaging/compression runtime yield identical archives. Binary compilation reproducibility is not asserted.',
    }
    payload[MANIFEST] = (json.dumps(manifest, indent=2, sort_keys=True) + '\n').encode('utf-8')
    args.output.mkdir(parents=True, exist_ok=True)
    archive = args.output / (archive_root + '.tar.gz')
    checksum = archive.with_name(archive.name + '.sha256')
    require(not archive.exists() and not checksum.exists(), 'refusing to replace an existing archive/checksum')
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=str(args.output), prefix='.package-', delete=False) as raw:
            temporary = Path(raw.name)
            with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=epoch, compresslevel=9) as compressed:
                with tarfile.open(fileobj=compressed, mode='w', format=tarfile.PAX_FORMAT) as output:
                    directories = {archive_root}
                    for name in payload:
                        directories.update(str(parent) for parent in PurePosixPath(archive_root, name).parents
                                           if str(parent) != '.')
                    for directory in sorted(directories):
                        output.addfile(archive_info(directory, 0o755, epoch, directory=True))
                    for name, data in sorted(payload.items()):
                        output.addfile(archive_info(archive_root + '/' + name, mode_for(name), epoch, len(data)), io.BytesIO(data))
        os.link(str(temporary), str(archive))  # Atomic publication without overwriting another build.
        with checksum.open('x', encoding='utf-8') as out:
            out.write(sha(archive.read_bytes()) + '  ' + archive.name + '\n')
    finally:
        if temporary is not None:
            temporary.unlink()
    verify(archive, args.run_example, args.run_study_example)
    print('Package: ' + str(archive.resolve()))
    print('Checksum: ' + str(checksum.resolve()))


def read_verified(archive):
    checksum = archive.with_name(archive.name + '.sha256')
    parts = checksum.read_text(encoding='utf-8').strip().split()
    require(len(parts) == 2 and parts[1] == archive.name, 'invalid adjacent archive checksum')
    require(parts[0] == sha(archive.read_bytes()), 'archive SHA-256 does not match')
    files, modes, names = {}, {}, set()
    root = None
    with tarfile.open(str(archive), mode='r:gz') as source:
        for member in source:
            path = safe_path(member.name.rstrip('/') if member.isdir() else member.name)
            require(member.name not in names, 'duplicate archive member: ' + member.name)
            names.add(member.name)
            root = root or path.parts[0]
            require(path.parts[0] == root, 'archive contains multiple top-level directories')
            require(member.uid == 0 and member.gid == 0, 'archive contains host user/group identities')
            require(member.mode in (0o644, 0o755), 'unsafe or unexpected archive mode')
            if member.isdir():
                continue
            require(member.isfile() and len(path.parts) > 1, 'archive contains a link or special file')
            name = '/'.join(path.parts[1:])
            data = source.extractfile(member)
            require(data is not None, 'archive file cannot be read: ' + name)
            files[name], modes[name] = data.read(), member.mode
    require(MANIFEST in files, 'package manifest is missing')
    manifest = json.loads(files.pop(MANIFEST).decode('utf-8'))
    require(manifest['schema'] == SCHEMA and manifest['archive_root'] == root, 'unexpected package manifest identity')
    require(archive.name == root + '.tar.gz', 'archive filename and manifest root differ')
    require(set(files) == set(manifest['files']), 'archive file inventory differs from manifest')
    for name, data in files.items():
        expected = manifest['files'][name]
        require(sha(data) == expected['sha256'] and len(data) == expected['bytes'], 'payload hash/size differs: ' + name)
        require('%04o' % modes[name] == expected['mode'], 'payload mode differs: ' + name)
    require(modes.get('pgsum') == 0o755, 'packaged binary is not executable')
    require('examples/robustness/run.py' in files and 'docs/robustness.md' in files, 'required offline example/docs are missing')
    return manifest, files, modes


def verify(archive, run_example=False, run_study_example=False):
    archive = archive.resolve()
    manifest, files, modes = read_verified(archive)
    link_count = check_markdown(files)
    print('Verified %d payload hashes and %d local Markdown links: %s' % (len(files), link_count, archive.name))
    if run_example or run_study_example:
        # Both spaces in the extraction path and an unrelated cwd catch accidental
        # dependence on the build checkout or shell working directory.
        with tempfile.TemporaryDirectory(prefix='pgsum package check ') as temporary:
            directory = Path(temporary)
            extracted = directory / manifest['archive_root']
            for name, data in files.items():
                path = extracted.joinpath(*safe_path(name).parts)
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(data)
                path.chmod(modes[name])
            unrelated = directory / 'different working directory'
            unrelated.mkdir()
            binary = extracted / 'pgsum'
            version = subprocess.check_output([str(binary), '--version'], cwd=str(unrelated), text=True).strip()
            require(version == 'pgsum ' + manifest['version'], 'binary version differs from package manifest')
            subprocess.run([sys.executable, str(extracted / 'examples/robustness/run.py'),
                            '--binary', str(binary), '--out', str(directory / 'offline example output')] + (['--study'] if run_study_example else []),
                           cwd=str(unrelated), check=True)
        print('Extracted-binary offline example passed.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest='command', required=True)
    package = subcommands.add_parser('build', help='create archive and checksum, then verify them')
    package.add_argument('--binary', type=Path, required=True)
    package.add_argument('--target', required=True)
    package.add_argument('--output', type=Path, default=Path('dist'))
    package.add_argument('--source-root', type=Path, default=Path(__file__).resolve().parents[1])
    package.add_argument('--source-date-epoch', default=os.environ.get('SOURCE_DATE_EPOCH', '0'))
    package.add_argument('--run-example', action='store_true', help='execute archive binary; native targets only')
    package.add_argument('--run-study-example', action='store_true', help='also run the paired study example (native targets, Python 3.9+)')
    inspection = subcommands.add_parser('verify', help='check archive, manifest, paths, modes, and Markdown links')
    inspection.add_argument('archive', type=Path)
    inspection.add_argument('--run-example', action='store_true', help='execute archive binary; native targets only')
    inspection.add_argument('--run-study-example', action='store_true', help='also run the paired study example (native targets, Python 3.9+)')
    args = parser.parse_args()
    if args.command == 'build':
        build(args)
    else:
        verify(args.archive, args.run_example, args.run_study_example)


if __name__ == '__main__':
    main()
