#!/usr/bin/env python
#==============================================================================
#
#    File: 
#    correlate.1.py
#        
#    Usage: 
#    ./correlate.1.py -1 <infile 1> -2 <infile 2> -o <outfile prefix>
#
#    Description:
#    Reads two files with the same number of columns and rows. It is assumed
#    the first row is a header and that the first column contains ids that are
#    present in both files. The numbers will be correlated using R's cor.test
#    function. The results of the correlation are written to outfiles with the
#    specified prefix.
#
#    Date:
#    v1:2014-09-17
#
#    Author:
#    Johannes Asplund Samuelsson
#
#==============================================================================


from optparse import OptionParser
import sys

parser = OptionParser()
parser.add_option("-1", "--infile1", dest="in1",
        action="store", type="string",
        help="Read from infile 1.")
parser.add_option("-2", "--infile2", dest="in2",
        action="store", type="string",
        help="Read from infile 2.")
parser.add_option("-o", "--outfile", dest="outfile",
        action="store", type="string",
        help="Write to outfile.")

(options, args) = parser.parse_args()

if not options.in1 or not options.in2 or not options.outfile: sys.exit(parser.print_help())


# Read both infiles


print "Reading infile 1."

n = 0

data1 = {}

for line in open(options.in1):
    n += 1
    if n == 1: continue
    line = line.rstrip()
    line = line.split("\t")
    data1[line[0]] = line[1:]


print "Reading infile 2."

n = 0

data2 = {}

for line in open(options.in2):
    n += 1
    if n == 1: continue
    line = line.rstrip()
    line = line.split("\t")
    data2[line[0]] = line[1:]


# Confirm that ids match

if set(data1.keys()) != set(data2.keys()):
    sys.exit("Error: Row identifiers do not match.")


# Perform correlation tests

import rpy2.robjects as R
import pandas as pd
from ggplot import *


cor_test = R.r['cor.test']

outfile = open(options.outfile+".tab", 'w')

header = "id\tp (Pearson)\tcor (Pearson)\tp (Spearman)\trho (Spearman)\n"
outfile.write(header)

out_df_list = []

for identifier in data1.keys():

    # Make R FloatVectors
    vector1 = R.FloatVector(data1[identifier])
    vector2 = R.FloatVector(data2[identifier])

    # Create a pandas DataFrame
    df = pd.DataFrame(data = zip(vector1, vector2), columns = ['x','y'])
   
    # Determine min and max values for y and x axis
    ymax = float(max(df['y'].tolist()))
    ymin = float(min(df['y'].tolist()))
    xmax = float(max(df['x'].tolist()))
    xmin = float(min(df['x'].tolist()))

    # Perform Spearman and Pearson correlation tests
    result_spearman = cor_test(vector1, vector2, method='spearman')
    result_pearson = cor_test(vector1, vector2, method='pearson')

    p_value_spearman = result_spearman.rx(3)[0][0]
    p_value_pearson = result_pearson.rx(3)[0][0]

    rho_value_spearman = result_spearman.rx(4)[0][0]
    cor_value_pearson = result_pearson.rx(4)[0][0]
    
    # Assemble a title with the identifier and the test statistics
    title = '%s\np=%.4f, cor=%.4f (Pearson)\np=%.4f, rho=%.4f (Spearman)' % (identifier, p_value_pearson, cor_value_pearson, p_value_spearman, rho_value_spearman)
    print title
    print ""
    
    # Gather data for facet_wrap output of all correlations in the end
    for value1, value2 in zip(vector1, vector2):
        out_df_tuple = (value1, value2, title)
        out_df_list.append(out_df_tuple)
    
    # Create ggplot for the current correlation
    plot = ggplot(aes(x='x', y='y'), data=df) + geom_point(colour='steelblue', size = 200) + geom_smooth(method='lm', colour='steelblue', size = 5) + xlim(xmin,xmax) + ylim(ymin,ymax) + ggtitle(title) + scale_x_continuous(title='') + scale_y_continuous(title='')

    # Save plot to pdf file
    ggsave(filename = options.outfile+identifier+".pdf", plot = plot)

    # Write data to outfile
    text_output = '\t'.join([identifier, str(p_value_pearson), str(cor_value_pearson), str(p_value_spearman), str(rho_value_spearman)])+'\n'
    outfile.write(text_output)

outfile.close()


# Create a pandas DataFrame for the facet_wrap output
out_df = pd.DataFrame(data = out_df_list, columns = ['x', 'y', 'title'])
    

# Create facet_wrap plot (currently not functional)
plot = ggplot(aes('x', 'y'), data=out_df) + geom_point(colour='steelblue') + geom_smooth(method='lm', colour='steelblue') + facet_wrap('title', ncol=7)

