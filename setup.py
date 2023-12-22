from setuptools import setup

setup(name='semirings',
      version='1.0',
      description='Semirings are a powerful abstraction for dynamic programming.',
      author='Tim Vieira',
      project_url = 'https://github.com/timvieira/semirings',
      install_requires = [
          'numpy',
          'IPython',
          'graphviz',
          'pytest',
          'arsenal',
      ],
      readme=open('README.md').read(),
      packages = ['semirings'],
)
