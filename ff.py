from scipy.special import beta as Beta
from math import pow
mu0 = 1
## mu0charm = 1.43
## mu0bottom = 4.3

## this is for pi+
N = {
    'u+ubar':   0.345,
    'd+dbar':   0.380,
    'ubar'  :   0.115,
    's+sbar':   0.190,
    'c+cbar':   0.271,
    'b+bbar':   0.501,
    'g':        0.279
}

alpha = {
    'u+ubar':   -0.015,
    'd+dbar':   -0.015,
    'ubar'  :    0.520,
    's+sbar':    0.520,
    'c+cbar':   -0.905,
    'b+bbar':   -1.305,
    'g':         0.899
}

beta = {
    'u+ubar':   1.20,
    'd+dbar':   1.20,
    'ubar'  :   3.27,
    's+sbar':   3.27,
    'c+cbar':   3.23,
    'b+bbar':   5.67,
    'g':        1.57
}

gamma = {
    'u+ubar':   11.06,
    'd+dbar':   11.06,
    'ubar'  :   16.26,
    's+sbar':   16.26,
    'c+cbar':   0.00,
    'b+bbar':   0.00,
    'g':        20.00
}

delta = {
    'u+ubar':   4.23,
    'd+dbar':   4.23,
    'ubar'  :   8.46,
    's+sbar':   8.46,
    'c+cbar':   0.00,
    'b+bbar':   0.00,
    'g':        4.91
}

def dss(z, flavor, charge):
    i = flavor
    num = N[i] * pow(z,alpha[i]) * pow(1-z,beta[i]) * \
        (1+gamma[i]*pow(1-z, delta[i]))
    denom = Beta(2+alpha[i], beta[i]+1) + \
        gamma[i]*Beta(2+alpha[i], beta[i]+delta[i]+1)
    return num/denom


def pythia(ckMin, ckMax, nevents):
    """
    runs standalone Pythia to measure charged pion FFs
    """
    import ROOT
    pythia = ROOT.TPythia6()
    
    # CDF Tune A for STAR
    pythia.SetMSEL(1)
    pythia.SetMSTP(51, 7)
    pythia.SetMSTP(81, 1)
    pythia.SetMSTP(82, 4)
    pythia.SetPARP(67, 4.0)
    pythia.SetPARP(83, 0.5)
    pythia.SetPARP(84, 0.4)
    pythia.SetPARP(85, 0.9)
    pythia.SetPARP(86, 0.95)
    pythia.SetPARP(89, 1800)
    pythia.SetPARP(90, 0.25)
    pythia.SetPARP(91, 1.0)
    
    pythia.SetCKIN(3, ckMin)
    pythia.SetCKIN(4, ckMax)
    
    pythia.Initialize('CMS', 'p', 'p', 200)
    
    h = {
        211: {
            1:  ROOT.TH1D('d_plus', '', 200, 0., 1.),
            2:  ROOT.TH1D('u_plus', '', 200, 0., 1.),
            3:  ROOT.TH1D('s_plus', '', 200, 0., 1.),
            4:  ROOT.TH1D('c_plus', '', 200, 0., 1.),
            5:  ROOT.TH1D('b_plus', '', 200, 0., 1.),
            -1: ROOT.TH1D('dbar_plus', '', 200, 0., 1.),
            -2: ROOT.TH1D('ubar_plus', '', 200, 0., 1.),
            -3: ROOT.TH1D('sbar_plus', '', 200, 0., 1.),
            -4: ROOT.TH1D('cbar_plus', '', 200, 0., 1.),
            -5: ROOT.TH1D('bbar_plus', '', 200, 0., 1.),
            21: ROOT.TH1D('g_plus', '', 200, 0., 1.)
        },
        -211: {
            1:  ROOT.TH1D('d_minus', '', 200, 0., 1.),
            2:  ROOT.TH1D('u_minus', '', 200, 0., 1.),
            3:  ROOT.TH1D('s_minus', '', 200, 0., 1.),
            4:  ROOT.TH1D('c_minus', '', 200, 0., 1.),
            5:  ROOT.TH1D('b_minus', '', 200, 0., 1.),
            -1: ROOT.TH1D('dbar_minus', '', 200, 0., 1.),
            -2: ROOT.TH1D('ubar_minus', '', 200, 0., 1.),
            -3: ROOT.TH1D('sbar_minus', '', 200, 0., 1.),
            -4: ROOT.TH1D('cbar_minus', '', 200, 0., 1.),
            -5: ROOT.TH1D('bbar_minus', '', 200, 0., 1.),
            21: ROOT.TH1D('g_minus', '', 200, 0., 1.)
        }
    }
    
    def makeVector(pythia, i):
        return ROOT.TLorentzVector( pythia.GetP(i,1), pythia.GetP(i,2),
            pythia.GetP(i,3), pythia.GetP(i,4) )
    
    
    ## Particle numbers are as follows:
    ## 1,2 protons
    ## 3,4 incoming partons before ISR
    ## 5,6 incoming partons after ISR
    ## 7,8 outgoing partons  ... use these?
    for i in range(nevents):
        if i % 1000 == 0: print 'generating event', i
        pythia.GenerateEvent()
        nparticles = pythia.GetN()
        
        for i in range(9,nparticles+1):
            pid = pythia.GetK(i, 2)
            
            ## select only charged pions
            if abs(pid) != 211: continue
            
            status = pythia.GetK(i, 1)
            assert(status == 1)
            
            parent = pythia.GetK(i, 3)            
            while parent > 8:
                parent = pythia.GetK(parent, 3)
            
            ## skip pions from beam remnants, ISR
            if parent not in (7,8): continue
            
            mom = makeVector(pythia, i)
            
            parton = makeVector(pythia, parent)
            flavor = pythia.GetK(parent, 2)
            
            z = mom.Vect().Dot(parton.Vect()) / parton.P()**2
            
            h[pid][flavor].Fill(z)
        
    pythia.Pystat(1)
    return h

    
