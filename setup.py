from setuptools import setup

setup(name='semirings',
      version='0.3.2',
      description='Semirings are a powerful abstraction for dynamic programming.',
      author='Tim Vieira',
      project_url = 'https://github.com/timvieira/semirings',
      install_requires = [
          'numpy',
          'scipy',       # for convex hull support
          'graphviz',
          'pytest',
          'arsenal @ git+https://github.com/timvieira/arsenal',
          'wfsa @ git+https://github.com/timvieira/wfsa',
          'fsa @ git+https://github.com/timvieira/fsa',
      ],
      long_description=open('README.md').read(),
      packages = ['semirings'],
)
