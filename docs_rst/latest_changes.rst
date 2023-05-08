Change log
==========

v2023.5.8
---------
* PR #2959 from @materialsproject dependabot/github_actions/pypa/cibuildwheel-2.12.3
    Bumps [pypa/cibuildwheel](https://github.com/pypa/cibuildwheel) from 2.11.4 to 2.12.3.
    <details>
    <summary>Release notes</summary>
    <p><em>Sourced from <a href="https://github.com/pypa/cibuildwheel/releases">pypa/cibuildwheel's releases</a>.</em></p>
    <blockquote>
    <h2>v2.12.3</h2>
    <ul>
    <li>🐛 Fix an import error when running on Python 3.7. (<a href="https://redirect.github.com/pypa/cibuildwheel/issues/1479">#1479</a>)</li>
    </ul>
    <h2>v2.12.2</h2>
    <ul>
    <li>🐛 Fix a bug that caused an extra empty config-setting to be passed to the backend when CIBW_BUILD_FRONTEND is set to <code>build</code>. (<a href="https://redirect.github.com/pypa/cibuildwheel/issues/1474">#1474</a>)</li>
    <li>🐛 Fix a crash that occurred when overwriting an existing wheel on Windows. (<a href="https://redirect.github.com/pypa/cibuildwheel/issues/1464">#1464</a>)</li>
    <li>🛠 Pinned version updates, including CPython 3.10.11, 3.11.3, pip 23.1 and wheel 0.40.0.</li>
    </ul>
    <h2>v2.12.1</h2>
    <ul>
    <li>🐛 Fix a bug that prevented the use of CIBW_CONFIG_SETTINGS with the 'pip' build backend. (<a href="https://redirect.github.com/pypa/cibuildwheel/issues/1430">#1430</a>)</li>
    </ul>
    <h2>v2.12.0</h2>
    <ul>
    <li>✨ Adds support for PyPy arm64 wheels. This means that you can build PyPy wheels for Apple Silicon machines. Cross-compilation is not supported for these wheels, so you'll have to build on an Apple Silicon machine. (<a href="https://redirect.github.com/pypa/cibuildwheel/issues/1372">#1372</a>)</li>
    <li>🛠 Pinned version updates, including PyPy to v7.3.11 and setuptools to 66.0.0.</li>
    </ul>
    </blockquote>
    </details>
    <details>
    <summary>Changelog</summary>
    <p><em>Sourced from <a href="https://github.com/pypa/cibuildwheel/blob/main/docs/changelog.md">pypa/cibuildwheel's changelog</a>.</em></p>
    <blockquote>
    <h3>v2.12.3</h3>
    <p><em>19 April 2023</em></p>
    <ul>
    <li>🐛 Fix an import error when running on Python 3.7. (<a href="https://redirect.github.com/pypa/cibuildwheel/issues/1479">#1479</a>)</li>
    </ul>
    <h3>v2.12.2</h3>
    <p><em>18 April 2023</em></p>
    <ul>
    <li>🐛 Fix a bug that caused an extra empty config-setting to be passed to the backend when CIBW_BUILD_FRONTEND is set to <code>build</code>. (<a href="https://redirect.github.com/pypa/cibuildwheel/issues/1474">#1474</a>)</li>
    <li>🐛 Fix a crash that occurred when overwriting an existing wheel on Windows. (<a href="https://redirect.github.com/pypa/cibuildwheel/issues/1464">#1464</a>)</li>
    <li>🛠 Pinned version updates, including CPython 3.10.11, 3.11.3, pip 23.1 and wheel 0.40.0.</li>
    </ul>
    <h3>v2.12.1</h3>
    <p><em>11 March 2023</em></p>
    <ul>
    <li>🐛 Fix a bug that prevented the use of CIBW_CONFIG_SETTINGS with the 'pip' build backend. (<a href="https://redirect.github.com/pypa/cibuildwheel/issues/1430">#1430</a>)</li>
    </ul>
    <h3>v2.12.0</h3>
    <p><em>16 Jan 2023</em></p>
    <ul>
    <li>✨ Adds support for PyPy arm64 wheels. This means that you can build PyPy wheels for Apple Silicon machines. Cross-compilation is not supported for these wheels, so you'll have to build on an Apple Silicon machine. (<a href="https://redirect.github.com/pypa/cibuildwheel/issues/1372">#1372</a>)</li>
    <li>🛠 Pinned version updates, including PyPy to v7.3.11 and setuptools to 66.0.0.</li>
    </ul>
    </blockquote>
    </details>
    <details>
    <summary>Commits</summary>
    <ul>
    <li><a href="https://github.com/pypa/cibuildwheel/commit/5e15bb25b428e1bf2daf2215f173d2b40135f56f"><code>5e15bb2</code></a> Bump version: v2.12.3</li>
    <li><a href="https://github.com/pypa/cibuildwheel/commit/52572ccdd01fcde11603283bcfda090086a7d11c"><code>52572cc</code></a> Merge pull request <a href="https://redirect.github.com/pypa/cibuildwheel/issues/1479">#1479</a> from henryiii/henryiii/fix/py37</li>
    <li><a href="https://github.com/pypa/cibuildwheel/commit/c6027e4fc72eee01851e06b1cd173aabeb3abfe0"><code>c6027e4</code></a> docs: fix options.md typos (<a href="https://redirect.github.com/pypa/cibuildwheel/issues/1477">#1477</a>)</li>
    <li><a href="https://github.com/pypa/cibuildwheel/commit/aff6dd5adc49114f1ee49c82bcb386a142fc614a"><code>aff6dd5</code></a> refactor: tomllib in _compat</li>
    <li><a href="https://github.com/pypa/cibuildwheel/commit/31bd9c91743e09f5b48c4def38c255c14c00d1af"><code>31bd9c9</code></a> refactor: restore typing for non-backports</li>
    <li><a href="https://github.com/pypa/cibuildwheel/commit/d996af554ae21435028e538e417e403e8db7433e"><code>d996af5</code></a> fix: restore Python 3.7 support</li>
    <li><a href="https://github.com/pypa/cibuildwheel/commit/e42f243b523a93352f0686ba2f4a62b7bf5f7d34"><code>e42f243</code></a> ci: add Python 3.7</li>
    <li><a href="https://github.com/pypa/cibuildwheel/commit/9f18d385cbb38b9c24f46a2c798f361ecea8fd5c"><code>9f18d38</code></a> Bump version: v2.12.2</li>
    <li><a href="https://github.com/pypa/cibuildwheel/commit/dfbc6c3902b689c250f9d9dd10ce8f88db8e2b7b"><code>dfbc6c3</code></a> Merge pull request <a href="https://redirect.github.com/pypa/cibuildwheel/issues/1476">#1476</a> from henryiii/henryiii/chore/ruff_ex</li>
    <li><a href="https://github.com/pypa/cibuildwheel/commit/8c5f89c035ba6ac81794ed5a5d6cd2110957f181"><code>8c5f89c</code></a> chore: update typing to be generic on function args</li>
    <li>Additional commits viewable in <a href="https://github.com/pypa/cibuildwheel/compare/v2.11.4...v2.12.3">compare view</a></li>
    </ul>
    </details>
    <br />
    [![Dependabot compatibility score](https://dependabot-badges.githubapp.com/badges/compatibility_score?dependency-name=pypa/cibuildwheel&package-manager=github_actions&previous-version=2.11.4&new-version=2.12.3)](https://docs.github.com/en/github/managing-security-vulnerabilities/about-dependabot-security-updates#about-compatibility-scores)
    Dependabot will resolve any conflicts with this PR as long as you don't alter it yourself. You can also trigger a rebase manually by commenting `@dependabot rebase`.
    [//]: # (dependabot-automerge-start)