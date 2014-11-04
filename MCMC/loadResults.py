from AXLmodel import *
import bayessb
from bayessb import multichain
from bayessb import convergence
from bayessb.report import *
from bayessb.report import reporters

if __name__ == '__main__':

    Afiles = list(["ASPATHSPATUSPATESPAT2SPATY"])
    #BTfiles = list(["12RAQ5M"])
    
    fullSet = dict( [('A549', Afiles)])#, ('BT549', BTfiles)] )
    Alist = list((
        reporters.marginals, reporters.maximum_likelihood,
        reporters.convergence_criterion, reporters.maximum_posterior))

    Areport = Report(fullSet, Alist, None, 1)

    Areport.write_html_table_with_links('out.html')
    
    

    print('Done.')
