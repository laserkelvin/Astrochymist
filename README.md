# Astrochymist Parsing

This is a script to pull information from the Astrochymist, a valuable source
of discoveries of molecules in space.

You can install the requirements needed by the script by running:

`pip install -r requirements.txt`

The script is written in Python 3, parsing with `BeautifulSoup` and outputting
with `Pandas`. The script can be called simply by running:

`python astro_parser.py`

or alternatively, by importing it as a function (if you are using Jupyter
notebooks).

The enclosed Jupyter notebook is an example.

The current version does not separate the chemical formulas from the names.
This is mainly due to the HTML, although in theory it could be done perhaps
with some regex.

