#!/usr/bin/python
"""
one-off macro for Bernd that plots PDFs from Werner
"""
import math
import sys
sys.path.append('/sw/lib/root')
import ROOT
ROOT.gStyle.SetOptStat(0)
from array import array

def integrate(hinput, logx=True):
    ## style stuff
    ROOT.gStyle.SetCanvasColor(10)
    ROOT.gStyle.SetFillColor(10)
    ROOT.gStyle.SetStatColor(0)
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetOptDate(1)
    
    xaxis = hinput.GetXaxis()
    xbins = [ xaxis.GetBinUpEdge(bin) for bin in range(xaxis.GetNbins()+1) ]
    if logx:
        xbins = [ math.pow(10, b) for b in xbins ]
    x = [ (xbins[j+1] + xbins[j])/2 for j in range(len(xbins)-1) ]
    xdg = [ hinput.GetBinContent(bin+1) for bin in range(xaxis.GetNbins()) ]
    dg = [ xdg[i]/x[i] for i in range(len(xdg)) ]
    rootbins = array('f', xbins)
    hnew = ROOT.TH1D('hnew', '', len(x), rootbins)
    hintegral = ROOT.TH1D('hintegral', '', len(x), rootbins)
    hnorm_integral = ROOT.TH1D('hnorm_integral', '', len(x), rootbins)
    [ hnew.SetBinContent(bin+1, dg[bin]) for bin in range(len(dg)) ]

    total_integral = hnew.Integral('width')
    for bin, xmin in enumerate(xbins[:-1]):
        integral = hnew.Integral(bin+1, len(xbins)-1, 'width')
        hintegral.SetBinContent(bin+1, integral)
        hnorm_integral.SetBinContent(bin+1, integral/total_integral)
    
    f = open('%s.out' % hinput.GetName(), 'w')
    f.write('   x      x*dg(x)      dg(x)       dG      dG_frac\n')
    f.write('--------------------------------------------------\n')
    
    for bin in range(len(dg)):
        f.write('%.5f   %.5f   % 10.5f   %.5f   %.5f\n' % \
            (x[bin], xdg[bin], dg[bin], hintegral.GetBinContent(bin+1), 
            hnorm_integral.GetBinContent(bin+1)) )
    f.write('\n')
    f.close()
    
    c = ROOT.TCanvas('c', '', 700, 700)
    c.Divide(2,2)
    pad = c.cd(1)
    hinput.Draw()
    hinput.SetXTitle('log(x)')
    hinput.SetYTitle('x * #Delta g(x)')
    pad = c.cd(2)
    hnew.Draw()
    hnew.SetXTitle('x')
    hnew.SetYTitle('#Delta g(x)')
    pad.SetLogx()
    pad = c.cd(3)
    hintegral.SetXTitle('x_{min}')
    hintegral.SetYTitle('#DeltaG')
    hintegral.Draw()
    pad.SetLogx()
    pad = c.cd(4)
    hnorm_integral.SetXTitle('x_{min}')
    hnorm_integral.SetYTitle('#DeltaG(x_{min} .. 1) / #DeltaG(10^-4 .. 1)')
    hnorm_integral.Draw()
    pad.SetLogx()
    raw_input('press enter to continue:')

if __name__ == '__main__':
    if len(sys.argv) == 3:
        f = ROOT.TFile(sys.argv[1])
        integrate(f.Get(sys.argv[2]))
    else:
        print 'Usage: ./pdfs.py filename histogram_name'

