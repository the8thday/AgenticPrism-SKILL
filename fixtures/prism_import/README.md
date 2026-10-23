# Public Prism extraction fixtures

25 unchanged public files: 1 `.prism` + 24 `.pzfx`. `sources.json` pins the
repository, commit, source path and SHA-256 of every file. Both upstream license
notices are retained in each subdirectory. These are test data; no upstream
parser code is copied into the runtime.

- [Biomiha/prism2R](https://github.com/Biomiha/prism2R/tree/3bef7980c9abf4d194d89ca15f796246c3b84a47):
  paired project examples (their save dates differ); native format `1-4-1`,
  created in Prism 10.4.0.527 and modified in 10.4.1.532 on Mac.
- [Yue-Jiang/pzfx](https://github.com/Yue-Jiang/pzfx/tree/ceb5dbb5331d702d0197b26f773816c339dbe835):
  22 test cases and the exponential-decay example.

Tests assert source hashes, every stored cell's text/position, flags, roles,
classification and unknown-data behavior. Synthetic mutations in tests are
engineering checks, never evidence of a native Prism version or GUI behavior.
Licenses are reproduced as supplied, including the upstream LICENSE/NOTICE
copyright-holder differences. Do not treat filenames as unique project titles.
