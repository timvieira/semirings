"""
The documentation page content/semirings.md is part of the test suite:

* every `>>>` example on the page must run and print exactly what the page says;
* every link of the form [`Name`](.../blob/main/<path>) must point at a file
  that defines `Name`, so links cannot rot when code moves.
"""
import doctest
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAGE = ROOT / 'content' / 'semirings.md'
REPO = 'https://github.com/timvieira/semirings/blob/main/'


def python_blocks(text):
    return re.findall(r'```python\n(.*?)```', text, re.DOTALL)


def test_examples():
    text = PAGE.read_text()
    blocks = python_blocks(text)
    assert blocks, 'no examples found'
    parser = doctest.DocTestParser()
    runner = doctest.DocTestRunner(optionflags=doctest.NORMALIZE_WHITESPACE)
    globs = {}                           # examples share state, top to bottom
    for i, block in enumerate(blocks):
        test = parser.get_doctest(block, globs, f'semirings.md[block {i}]', str(PAGE), 0)
        runner.run(test)
    result = runner.summarize(verbose=False)
    assert result.failed == 0, f'{result.failed} of {result.attempted} examples failed'


def test_links_point_at_definitions():
    text = PAGE.read_text()
    links = re.findall(r'\[`(\w+)`\]\(' + re.escape(REPO) + r'([^)#]+)\)', text)
    assert links, 'no library links found'
    for name, path in links:
        f = ROOT / path
        assert f.exists(), f'{name}: {path} does not exist'
        src = f.read_text()
        defined = re.search(rf'^(class|def)\s+{name}\b|^{name}\s*=', src, re.M)
        assert defined, f'{name} is not defined in {path}'


def test_file_links_exist():
    text = PAGE.read_text()
    for path in re.findall(re.escape(REPO) + r'([^)#\s]+)', text):
        assert (ROOT / path).exists(), path
