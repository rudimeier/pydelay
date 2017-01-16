from distutils.core import setup
setup(name='pydelay',
      version='0.1.1',
      author='Valentin Flunkert',
      author_email='flunkert@gmail.com',
      license='MIT',
      packages=['pydelay'],
      package_dir={'pydelay': 'pydelay'},
      package_data={'pydelay': ['doc/pyplots/*', 
          'doc/sphinxext/*', 
          'doc/Makefile', 
          'doc/conf.py', 
          'doc/README', 
          'doc/index.rst', 
          'doc/pydelay.pdf', 
          'examples/*.py']},
      )

