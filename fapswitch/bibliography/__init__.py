"""
fapswitch.references

A simple bibliography that has a single dictionary with all the included
references.

"""

class Reference(object):
    """Store information about a citation/reference."""
    def __init__(self, key, doi, title, author, journal, year):
         self.key = key
         self.doi = doi
         self.title = title
         self.author = author
         self.journal = journal
         self.year = year


references = {
    'Kadantsev2013' : Reference(
        'Kadantsev2013', '10.1021/jz401479k',
        'Fast and Accurate Electrostatics in Metal Organic '
        'Frameworks with a Robust Charge Equilibration Parameterization '
        'for High-Throughput Virtual Screening of Gas Adsorption',
        'Eugene S. Kadantsev, Peter G. Boyd, Thomas D. Daff, '
        'and Tom K. Woo', 'The Journal of Physical Chemistry Letters',
        '2013'),
    'Blank2015' : Reference(
        'Blank2015', '0.1021/jz123456z',
        'Some uninteresitng paper about things',
        'A. Guy, A. Student',
        'The Journal of Things',
        '2015'
    ),
    'Blank2017' : Reference(
        'Blank2017', '0.1055/jz123456z',
        'Synthesis of materials to make the best MOFS',
        'M. Post Doc, A. Student, T. E. Knetitian',
        'The Journal of Synthesis',
        '2017'
    )
}