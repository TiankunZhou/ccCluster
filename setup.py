__author__ = "Rita Giordano, Gianluca Santoni, Tiankun Zhou"
__copyright__ = "Copyright 2015-2026"
__credits__ = ["Rita Giordano, Gianluca Santoni, Tiankun Zhou, Alexander Popov"]
__license__ = ""
__version__ = "2.0"
__maintainer__ = "Tiankun Zhou"
__email__ = "tiankun.zhou@esrf.fr"
__status__ = "Beta"


from setuptools import setup, find_packages

setup(name="ccCluster",
      version="2.0",
      description='Hierarchical cluster analysis for MX',
      author='Gianluca Santoni, Tiankun Zhou',
      packages=["ccCluster"],
      setup_requires=['setuptools'],
      install_requires=[
          'numpy',
          'scipy',
          'pathos',
          'matplotlib',
          'PySide6'],  # TODO: iotbx
      entry_points={
          'console_scripts': [
              'ccCalc = ccCluster.ccCalc:main',
              'HCAPipeline = ccCluster.HCAPipeline:main',
              'ccCluster = ccCluster.ccCluster:main',
	      'ccCluster-gui = ccCluster.ccCluster_gui:main'],
          },
      )
