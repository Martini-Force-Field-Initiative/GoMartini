from distutils.core import setup

setup(
  name = 'YOURPACKAGENAME',         # How you named your package folder (MyLib)
  packages = ['YOURPACKAGENAME'],   # Chose the same as "name"
  version = '0.1',
  license='Apache 2.0',
  description = 'TYPE YOUR DESCRIPTION HERE',   # Give a short description about your library
  author = 'L. Borges-Araújo',
  author_email = 'lpborara@gmail.com',
  url = 'https://github.com/user/reponame'   # Provide either the link to your github or to your website
  download_url = 'https://github.com/user/reponame/archive/v_01.tar.gz',    # I explain this later on
  keywords = ['martini', 'MD', 'Martinize', 'go model'],   
  install_requires=['numpy',],
  classifiers=[
    'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Environment :: Console',
    'Intended Audience :: Developers',      
    'Intended Audience :: Science/Research'
    'Topic :: Software Development :: Build Tools',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry',
    'License :: OSI Approved :: Apache Software License',   
    'Programming Language :: Python :: 3',     
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
  ],
  scripts=['',]
)