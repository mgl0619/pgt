try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
      'description':'Tandem Mass Spectrometry Analysis Tool',
      'author':'Gang Liu',
      'url':'to be determined',
      'download_url': 'to be de determined',
      'author_email': 'gangliu2001@gmail.com',
      'version':'0.01',
      'install_requires':['nose'],
      'packages':['pgnat'],
      'scripts':[],
      'name':'pgnat'
    }
setup(**config)
   
