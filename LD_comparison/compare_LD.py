
import argparse

parser = argparse.ArgumentParser(description='Compare LD between metabrain and eqtlgen')
parser.add_argument('metaBrain', help='File with the metaBrain results')
parser.add_argument('eqtlGen', help='File with the eqtlGen results')

args = parser.parse_args()
