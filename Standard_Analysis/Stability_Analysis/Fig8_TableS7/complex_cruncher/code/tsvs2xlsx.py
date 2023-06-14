#!/usr/bin/env python3

# Packages
import argparse
import pandas
import os


# Paramaters parser
parser = argparse.ArgumentParser(description='Combine multiple tsv files into a single xlsx workbook file.',)

parser.add_argument('-i', metavar='paths', dest='fins', nargs='+', required=True, action='store', help=('List of tsv input files. Example: $HOME/myprojects/*.tsv'))
parser.add_argument('-o', metavar='xlsx',  dest='fout', required=True, action='store', help=('Output xlsx file name'))


# Get script parameters
args = parser.parse_args()

# Create the xlsx workbook
writer = pandas.ExcelWriter(args.fout, engine='openpyxl')

for finput in args.fins:

    print(finput)
    
    # Get tsv file name without the extension
    (f_path, f_name) = os.path.split(finput)
    (f_short_name, f_extension) = os.path.splitext(f_name)

    # Read the tsv file
    df = pandas.read_table(finput)
        
    # Convert the table into a sheet
    df.to_excel(writer, index=False, sheet_name=f_short_name)


# Save the xlsx workbook into a file
writer.save()

