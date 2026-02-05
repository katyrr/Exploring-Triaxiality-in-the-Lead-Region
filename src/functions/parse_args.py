from argparse import ArgumentParser

parser = ArgumentParser(add_help=False)

parser.add_argument("data_subfolder")
parser.add_argument('-v', '--verbose', action='store_true')

parser.add_argument('-d', '--display-figures', action='store_true')

args = parser.parse_args()