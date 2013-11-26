"""
fapswitch.references

A simple bibliography that has a single dictionary with all the included
references.

"""

from collections import namedtuple

Reference = namedtuple('Reference',
                       ['key', 'doi', 'title', 'author', 'journal', 'year'])


references = {
    'Kadantsev2013': Reference(
        key='Kadantsev2013',
        doi='10.1021/jz401479k',
        title='Fast and Accurate Electrostatics in Metal Organic '
        'Frameworks with a Robust Charge Equilibration Parameterization '
        'for High-Throughput Virtual Screening of Gas Adsorption',
        author='Eugene S. Kadantsev, Peter G. Boyd, Thomas D. Daff, '
        'and Tom K. Woo',
        journal='The Journal of Physical Chemistry Letters',
        year='2013'),
    'Ertl2009': Reference(
        key='Ertl2009',
        doi='10.1186/1758-2946-1-8',
        title='Estimation of synthetic accessibility score of drug-like '
        'molecules based on molecular complexity and fragment contributions',
        author='Peter Ertl and Ansgar Schuffenhauer',
        journal='Journal of Cheminformatics',
        year='2009'),
    'Daff2013': Reference(
        key='Daff2013',
        doi='10.xxxxxx/xxxxxxxxx',
        title='unpublished fapswitch method',
        author='Thomas D Daff and Tom K Woo',
        journal='In Preparation',
        year='2009')
}
