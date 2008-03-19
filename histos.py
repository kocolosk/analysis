import math
import os
##import uuid 
import time
from array import array

import ROOT
import minimc

pidCalibration = {
6988 : ( 0.066608, 0.889596),
6990 : ( 0.028513, 0.872612),
6992 : ( 0.016085, 0.856032),
6994 : ( 0.076727, 0.867838),
6995 : ( 0.167452, 0.855377),
6997 : ( 0.206903, 0.882462),
6998 : (-0.018938, 0.844754),
7001 : (-0.006526, 0.863949),
7002 : ( 0.014921, 0.865907),
7032 : (-0.059741, 0.856267),
7034 : (-0.064317, 0.870471),
7035 : (-0.090822, 0.871124),
7048 : (-0.387454, 0.866186),
7049 : (-0.393791, 0.865776),
7051 : (-0.137894, 0.792171),
7055 : (-0.420703, 0.858151),
7064 : (-0.283141, 0.878024),
7067 : (-0.087372, 0.873845),
7068 : (-0.077816, 0.840146),
7069 : ( 0.010756, 0.861486),
7070 : (-0.054644, 0.879342),
7072 : (-0.022931, 0.863193),
7075 : (-0.065000, 0.814292),
7079 : (-0.261988, 0.896589),
7085 : (-0.207958, 0.879845),
7087 : (-0.062486, 0.818310),
7088 : (-0.184964, 0.869739),
7092 : (-0.095016, 0.879631),
7102 : (-0.161471, 0.899666),
7103 : (-0.180760, 0.900305),
7110 : (-0.194342, 0.894236),
7112 : (-0.116899, 0.878893),
7114 : (-0.140619, 0.871448),
7118 : (-0.206615, 0.883038),
7120 : (-0.096027, 0.884061),
7122 : (-0.187495, 0.903031),
7123 : (-0.172261, 0.885190),
7124 : (-0.175057, 0.906155),
7125 : (-0.134531, 0.887071),
#7127 : (-0.149092, 0.897838), ## RunLog_onl problem
7128 : (-0.149092, 0.897838),
7131 : (-0.197370, 0.889249),
7133 : (-0.207133, 0.888365),
#7134 : (-0.130486, 0.893232), ## RunLog_onl problem
7134 : (-0.119566, 0.901472),
7136 : (-0.026902, 0.621362),
7138 : (-0.150951, 0.892991),
7151 : ( 0.033011, 0.851425),
7153 : (-0.094305, 0.892003),
7154 : ( 0.010605, 0.824362),
7161 : (-0.167815, 0.894401),
7162 : (-0.152793, 0.906894),
7164 : (-0.125891, 0.894419),
7165 : (-0.162527, 0.870979),
7166 : (-0.133315, 0.884772),
7172 : (-0.207923, 0.894637),
7232 : (-0.110842, 0.864740),
7237 : (-0.269559, 0.916554),
7238 : (-0.259372, 0.891293),
7249 : (-0.247810, 0.894708),
7250 : (-0.316081, 0.909442),
7253 : (-0.271044, 0.845922),
7255 : (-0.249242, 0.857201),
7265 : (-0.190056, 0.865755),
7266 : (-0.243196, 0.869884),
7269 : (-0.242802, 0.892759),
7270 : (-0.314700, 0.886310),
7271 : (-0.389984, 0.910892),
7272 : (-0.194932, 0.883244),
7274 : (-0.344731, 0.903323),
7276 : (-0.198418, 0.897573),
7278 : (-0.273830, 0.875383),
7279 : (-0.219441, 0.895625),
7293 : (-0.346773, 0.803589), ## transverse
7295 : (-0.268542, 0.859643), ## transverse
7296 : (-0.328538, 0.901968), ## transverse
7300 : (-0.283700, 0.876021),
7301 : (-0.325708, 0.891317),
7302 : (-0.333164, 0.890520),
7303 : (-0.274641, 0.894173),
7304 : (-0.306279, 0.890682),
7305 : ( 0.726845, 0.934432),
7308 : (-0.297895, 0.913868),
7311 : (-0.232436, 0.905442),
7317 : (-0.301910, 0.886113),
7320 : (-0.229533, 0.893824),
7325 : (-0.251390, 0.898416),
7327 : (-0.251943, 0.895118),
7722 : (-0.059290, 0.867884), ## transverse
7724 : (-0.098639, 0.876645), ## transverse
7725 : (-0.076397, 0.857092), ## transverse
7729 : (-0.174965, 0.901068), ## transverse
7740 : (-0.075899, 0.890939), ## transverse
7847 : (-0.027059, 0.874814),
7850 : (-0.123179, 0.900492),
7851 : (-0.036289, 0.889293),
7852 : (-0.045614, 0.891709),
7853 : (-0.119185, 0.898291),
7855 : (-0.134738, 0.908692),
7856 : (-0.089265, 0.895362),
7858 : (-0.102251, 0.897319),
7863 : (-0.083709, 0.911717),
7864 : (-0.113995, 0.904617),
7865 : ( 0.009516, 0.813364),
7871 : (-0.109983, 0.906486),
7872 : (-0.138412, 0.882469),
7883 : (-0.182051, 0.885964),
7886 : (-0.185474, 0.889406),
7887 : (-0.155607, 0.886540),
7889 : (-0.100630, 0.863594),
7890 : (-0.053211, 0.921421),
7891 : (-0.099581, 0.861152),
7892 : (-0.077996, 0.892110),
7893 : (-0.123162, 0.907638),
7896 : (-0.094806, 0.894019),
7898 : (-0.057477, 0.896248),
7901 : (-0.014140, 0.865671),
7908 : (-0.040456, 0.883598),
7909 : (-0.077470, 0.880354),
7911 : (-0.048150, 0.861430),
7913 : (-0.076636, 0.896822),
7916 : (-0.078441, 0.894539),
7918 : ( 0.043434, 0.888361),
7921 : (-0.010669, 0.876733),
7922 : ( 0.022409, 0.907525),
7926 : (-0.033154, 0.886284),
7944 : (-0.068744, 0.905098),
7949 : (-0.039570, 0.876990),
7951 : ( 0.056927, 0.936951),
7952 : ( 0.025446, 0.885525),
7954 : ( 0.205871, 0.864685),
7957 : ( 0.014065, 0.899964)
}

class EventCuts:
    def __init__(self, event=None):
        if event is None:
            self.vertex  = False
            self.bbc         = False
            self.all         = False
        else:
            self.set(event)
    
    
    def set(self, event):
        self.vertex = event.nVertices() > 0
        
        timeDiff = event.bbcTimeBin()
        #if timeDiff % 32 != 0:
        #    self.bbc = False
        #else:
        #    bin = timeDiff / 32
        #    self.bbc = bin in (7,8,9) or (event.runId() > 7000000 and bin == 6)
        
        ## dropping the discrete timebin cut -- APK 2008-01-09
        bin = timeDiff / 32
        self.bbc = bin in (7,8,9) or (event.runId() > 7000000 and bin == 6)
        
        self.all = self.vertex and self.bbc
    
    


class TrackCuts:
    def __init__(self, fill):
        self.eta = False
        self.dca = False
        self.fit = False
        self.pid = False
        self.all = False
        
        ## need this try..except to bootstrap PID calibration
        try:
            pidFit = pidCalibration[fill]
        except KeyError:
            pidFit = (0.0, 1.0)
        #pidFit = pidCalibration[fill]
        self.pidMin = pidFit[0] - 1.0*pidFit[1]
        self.pidMax = pidFit[0] + 2.0*pidFit[1]
        
        ## try a "cleaner" bg sample at least 2 sigma from pion mean
        self.pid_bg = False
        self.pidBgMin = pidFit[0] - 2.0*pidFit[1]
        self.pidBgMax = self.pidMax
    
    
    def set(self, track):
        self.eta = math.fabs( track.eta() ) < 1.0
        self.dca = math.fabs( track.globalDca().mag() ) < 1.0
        self.fit = track.nHitsFit() > 25
        
        self.pid = self.pidMin < track.nSigmaPion() < self.pidMax
        self.pid_bg = (track.nSigmaPion() < self.pidBgMin) or \
                      (track.nSigmaPion() > self.pidBgMax)
        
        self.all = self.eta and self.dca and self.pid and self.fit
        


class JetCuts:
    triggerThresholds = { 96201:13, 96211:17, 96221:66, 96233:83, 137221:58, 137222:60 }
    minPhi2005 = 40
    maxPhi2005 = 320
    patchPhi2005 = [90., 30., -30., -90., -150., 150.]
    minPhi2006 = 36
    maxPhi2006 = 324
    patchPhi2006 = [150., 90., 30., -30., -90., -150., 150., 90., 30., -30., -90., -150.]
    
    
    def __init__(self, jet=None, event=None):
        if jet is None:
            self.eta = False
            self.rt = False
            self.trig = []
        else:
            self.set(jet, event)
    
    
    def set(self, jet, event):
        if event.runId() < 7000000:
            self.eta = 0.2 < jet.detectorEta() < 0.8
            self.rt = 0.1 < (jet.tpcEtSum() / jet.Et()) < 0.9
            self.trig = []
            
            ## HT geometric trigger condition
            for particle in jet.particles():
                if particle.detectorId() == ROOT.kBarrelEmcTowerId:
                    adc = event.highTowerAdc(particle.index())
                    if adc > self.triggerThresholds[96211]:
                        self.trig.append(96201)
                        self.trig.append(96211)
                    elif adc > self.triggerThresholds[96201]:
                        self.trig.append(96201)
            
            ## JP geometric trigger condition
            for patchId in range(6):
                adc = event.jetPatchAdc(patchId)
                if adc > self.triggerThresholds[96221]:
                    dPhi = math.fabs( math.degrees(jet.Phi()) - self.patchPhi2005[patchId] )
                    if dPhi < self.minPhi2005 or dPhi > self.maxPhi2005:
                        self.trig.append(96221)
                        if adc > self.triggerThresholds[96233]:
                            self.trig.append(96233)
        else:
            self.eta = -0.7 < jet.detectorEta() < 0.9
            self.rt = (jet.tpcEtSum() / jet.Et()) < 0.85
            self.trig = []
        
            ## JP geometric trigger condition
            for patchId in range(12):
                adc = event.jetPatchAdc(patchId)
                if adc > self.triggerThresholds[137221]:
                    dPhi = math.fabs( math.degrees(jet.Phi()) - self.patchPhi2006[patchId] )
                    if dPhi < self.minPhi2006 or dPhi > self.maxPhi2006:
                        self.trig.append(137221)
                        if adc > self.triggerThresholds[137222]:
                            self.trig.append(137222)
        
    


class TrackHistogramCollection(dict):
    """histograms filled for each track"""
    mcPtBins             = minimc.MiniMcHistos.ptBins
    mcEtaBins        = minimc.MiniMcHistos.etaBins
    mcPhiBins        = minimc.MiniMcHistos.phiBins
    mcNFitBins       = minimc.MiniMcHistos.nFitBins
    mcPBins          = minimc.MiniMcHistos.pBins
    mcDEdxBins       = minimc.MiniMcHistos.dEdxBins
    mcDcaGBins       = minimc.MiniMcHistos.dcaGBins
    
    allKeys = ['pt', 'eta', 'phi', 'nHitsFit', 'dEdx', 'dcaG', 'nSigmaPion',
        'dphi_deta', 'z', 'z_away', 'pt_near', 'pt_away', 'pt_bg']
    
    def __init__(self, name, tfile=None, keys=None):
        if tfile is not None:
            for key in self.allKeys:
                if keys is None or key in keys:
                    self[key] = tfile.Get('%s_%s' % (name,key))
        else:
            self['pt'] = ROOT.TH1D('%s_pt' % (name,), '', self.mcPtBins[0], self.mcPtBins[1], self.mcPtBins[2])
            self['eta'] = ROOT.TH1D('%s_eta' % (name,), '', self.mcEtaBins[0], self.mcEtaBins[1], self.mcEtaBins[2])
            self['phi'] = ROOT.TH1D('%s_phi' % (name,), '', self.mcPhiBins[0], self.mcPhiBins[1], self.mcPhiBins[2])
            self['nHitsFit'] = ROOT.TH1D('%s_nHitsFit' % (name,), '', self.mcNFitBins[0], self.mcNFitBins[1], self.mcNFitBins[2])
            self['dEdx'] = ROOT.TH1D('%s_dEdx' % (name,), '', self.mcDEdxBins[0], self.mcDEdxBins[1], self.mcDEdxBins[2])
            self['dcaG'] = ROOT.TH1D('%s_dcaG' % (name,), '', self.mcDcaGBins[0], self.mcDcaGBins[1], self.mcDcaGBins[2])
            
            self['nSigmaPion'] = ROOT.TH1D('%s_nSigmaPion' % (name,), '', 240, -6.0, 6.0)
            
            self['dphi_deta'] = ROOT.TH2D('%s_dphi_deta' % (name,), '', \
                50, -1.5*math.pi, 0.5*math.pi, 50, -2.0, 2.0)
            z_xbins = [2.0, 3.0, 4.0, 5.5, 7.0, 10.0]
            ar = array('d',z_xbins)
            self['z'] = ROOT.TH2D('%s_z' % (name,), '', len(z_xbins)-1, ar, 50, 0., 1.)
            self['z_away'] = ROOT.TH2D('%s_z_away' % (name,), '', len(z_xbins)-1, ar, 50, 0., 1.)
            
            self['pt_away'] = ROOT.TH1D('%s_pt_away' % (name,), '', self.mcPtBins[0], self.mcPtBins[1], self.mcPtBins[2])
            self['pt_near'] = ROOT.TH1D('%s_pt_near' % (name,), '', self.mcPtBins[0], self.mcPtBins[1], self.mcPtBins[2])
            
            self['pt_bg'] = ROOT.TH1D('%s_pt_bg' % name, '', self.mcPtBins[0], self.mcPtBins[1], self.mcPtBins[2])
    
    
    def fillTrack(self, track, tcuts):
        if tcuts.eta and tcuts.dca and tcuts.fit and tcuts.pid:
            self['pt'].Fill(track.pt())
            self['phi'].Fill(track.phi())
            self['dEdx'].Fill(track.dEdx() * 1e7)
        
        if tcuts.eta and tcuts.dca and tcuts.fit and tcuts.pid_bg:
            self['pt_bg'].Fill(track.pt())
        
        if tcuts.dca and tcuts.fit and tcuts.pid:
            self['eta'].Fill(track.eta())
        
        if tcuts.eta and tcuts.fit and tcuts.pid:
            self['dcaG'].Fill(track.globalDca().mag())
        
        if tcuts.eta and tcuts.dca and tcuts.pid:
            self['nHitsFit'].Fill(track.nHitsFit())
            
        if tcuts.eta and tcuts.dca and tcuts.fit:
            self['nSigmaPion'].Fill(track.nSigmaPion())
    
    
    def fillTrackJetPair(self, track, tcuts, jet, jcuts):
        """all jet-pion correlation studies go here"""
        if tcuts.eta and tcuts.dca and tcuts.fit and tcuts.pid and jcuts.rt and jcuts.eta:
            deta = track.Eta() - jet.Eta()
            dphi = track.Phi() - jet.Phi()
            
            #normalize dphi
            if dphi < -1.5*math.pi:
                dphi = dphi + 2*math.pi
            if dphi > 0.5*math.pi:
                dphi = dphi - 2*math.pi
                  
            self['dphi_deta'].Fill(dphi, deta)
            
            #renormalize
            if dphi < -math.pi:
                dphi = dphi + 2*math.pi
            if dphi > math.pi:
                dphi = dphi - 2*math.pi
            
            dR = math.sqrt(deta**2 + dphi**2)
            z = track.Pt() / jet.Pt()
            if dR < 0.4: 
                self['z'].Fill(track.Pt(), z)
                self['pt_near'].Fill(track.Pt())
            elif dR > 1.5:
                self['z_away'].Fill(track.Pt(), z)
                self['pt_away'].Fill(track.Pt())
    
    
    def Add(self, other):
        [ h.Add(other[key]) for key, h in self.items() ]
    
    
    def Write(self):
        """make it persistent"""
        [ h.Write() for h in self.values() ]
    


class HistogramCollection(dict):
    """convenience class for dealing with multiple histos all filled the same way"""
    mcVzBins             = minimc.MiniMcHistos.vzBins
    allKeys = ['nVertices', 'vx_vy', 'vz', 'vzBBC', 'spinBit', 'bx7', 'bbc']
    
    def __init__(self, name, tfile=None, keys=None):
        super(HistogramCollection, self).__init__()
        if tfile is not None:
            for key in self.allKeys:
                if keys is None or key in keys:
                    self[key] = tfile.Get('%s_%s' % (name,key))
        else:
            self['nVertices'] = ROOT.TH1D('%s_nVertices' % (name,),'',15,-0.5,14.5)
            self['vx_vy'] = ROOT.TH2D('%s_vx_vy' % (name,),'',400,-2.0,2.0, 400,-2.0,2.0)
            self['vz'] = ROOT.TH1D('%s_vz' % (name,),'',self.mcVzBins[0], self.mcVzBins[1], self.mcVzBins[2])
            self['vzBBC'] = ROOT.TH1D('%s_vzBBC' % (name,),'',self.mcVzBins[0], self.mcVzBins[1], self.mcVzBins[2])
            self['spinBit'] = ROOT.TH1D('%s_spinBit' % (name,),'',17,0.5,16.5)
            self['bx7'] = ROOT.TH1D('%s_bx7' % (name,),'',128,-0.5,127.5)
            self['bbc'] = ROOT.TH1D('%s_bbc' % (name,),'',400,-0.5,399.5)
        
        self.tracks_plus    = TrackHistogramCollection('%s_plus' % (name,), tfile, keys)
        self.tracks_minus = TrackHistogramCollection('%s_minus' % (name,), tfile, keys)
        self.tracks_sum = TrackHistogramCollection('%s_sum' % (name,), tfile, keys)
    
    def fillEvent(self, event, ecuts):
        vertex = event.vertex(0)
        self['nVertices'].Fill(event.nVertices())
        self['bbc'].Fill(event.bbcTimeBin())
        
        if ecuts.vertex:
            self['vz'].Fill(vertex.z())
            
        if ecuts.vertex and ecuts.bbc:
            self['vx_vy'].Fill(vertex.x(), vertex.y())
            self['vzBBC'].Fill(vertex.z())
            self['spinBit'].Fill(event.spinBit())
            self['bx7'].Fill(event.bx7())
    
    def trackHistograms(self,charge):
        if charge == 1:  return self.tracks_plus
        if charge == -1: return self.tracks_minus
        if charge == 0:  return self.tracks_sum
        return None
    
    def Add(self, other):
        [ h.Add(other[key]) for key, h in self.items() ]
        self.tracks_plus.Add( other.tracks_plus )
        self.tracks_minus.Add( other.tracks_minus )
        self.tracks_sum.Add( other.tracks_sum )
    
    
    def Write(self):
        """make it persistent"""
        [ h.Write() for h in self.values() ]
        self.tracks_plus.Write()
        self.tracks_minus.Write()
        self.tracks_sum.Write()
    


class HistogramManager(dict):
    """generates histograms from StChargedPionEvent trees and reads them back from disk"""
    spinKeys = {5:'uu', 6:'du', 9:'ud', 10:'dd'}
    trigSetups = ('96011','96201','96211','96221','96233','hightower','jetpatch', 'alltrigs',\
          '117001','137221','137222','137611','137622')
    
    def __init__(self, tfile=None, keys=None):
        super(HistogramManager, self).__init__()
        
        if tfile is not None:
            self.name = os.path.basename(tfile.GetName())
        
        self.fill = 0
        
        self.allHistos = []
        
        for spin in ('uu','ud','du','dd','other','anyspin'):
            self[spin] = {}
            setattr(self, spin, self[spin])
            for trig in self.trigSetups:
                self.allHistos.append( HistogramCollection('_%s_%s' % (trig, spin), tfile, keys) )
                self[spin][trig] = self.allHistos[-1]
    
    
    def processEvent(self, event):
        """fill histograms with info from event"""
        ## spin sorting
        spin = 'other'
        if event.isSpinValid():
            try: spin = self.spinKeys[event.spinBit()]
            except KeyError: pass
        
        ## trigger selection -- convert transverse trigger IDs to longitudinal
        triggerOk = {}
        if event.runId() < 7000000:
            if event.isPolLong():
                for trigId in (96011, 96201, 96211, 96221, 96233):
                    triggerOk[str(trigId)] = event.isTrigger(trigId) and event.isSimuTrigger(trigId)
                    triggerOk['%d_hw' % (trigId,)] = event.isTrigger(trigId)
            elif event.isPolTrans():
                for trigId in (106011, 106201, 106211, 106221, 106233):
                    triggerOk[str(trigId-10000)] = event.isTrigger(trigId) and event.isSimuTrigger(trigId)
                    triggerOk['%d_hw' % (trigId-10000,)] = event.isTrigger(trigId)
            triggerOk['hightower'] = triggerOk['96011'] or triggerOk['96201'] or triggerOk['96211']
            triggerOk['jetpatch']  = triggerOk['96011'] or triggerOk['96221'] or triggerOk['96233']
            triggerOk['alltrigs']  = triggerOk['96011'] or triggerOk['96201'] or triggerOk['96211'] or triggerOk['96221'] or triggerOk['96233']
            activeTriggers = ('96011','96201','96211','96221','96233','hightower','jetpatch', 'alltrigs')
        else:
            if event.isPolLong():
                triggerOk['137611'] = event.isTrigger(137611) and event.isSimuTrigger(137611)
                triggerOk['137622'] = event.isTrigger(137622) and event.isSimuTrigger(137622)
                triggerOk['jetpatch'] = (event.isTrigger(137221) and event.isSimuTrigger(137221)) or \
                                        (event.isTrigger(137222) and event.isSimuTrigger(137222))
            elif event.isPolTrans():
                triggerOk['137611'] = event.isTrigger(127611) and event.isSimuTrigger(127611)
                triggerOk['137622'] = event.isTrigger(127622) and event.isSimuTrigger(127622)
                triggerOk['jetpatch'] = event.isTrigger(127221) and event.isSimuTrigger(127221)
            triggerOk['alltrigs'] = triggerOk['jetpatch'] or triggerOk['137611'] or triggerOk['137622']
            activeTriggers = ('137611', '137622', 'jetpatch', 'alltrigs')
          
        ## event-wise histograms
        ecuts = EventCuts(event)
        for trig in activeTriggers:
            if not triggerOk[trig]: continue
            self[spin][trig].fillEvent(event, ecuts)
        
        ## track-wise histograms
        if not ecuts.all: return
        tcuts = TrackCuts(self.fill)
        jcuts = [JetCuts(jet, event) for jet in event.jets()]
        for track in event.tracks():
            tcuts.set(track)
            for trig in activeTriggers:
                if not triggerOk[trig]: continue
                tcoll = self[spin][trig].trackHistograms(track.charge())
                tcoll.fillTrack(track, tcuts)
                
                ## track-jet correlations
                for i,jet in enumerate(event.jets()):
                    if 137222 in jcuts[i].trig:
                        if trig in ('137221', '137222', 'jetpatch', 'alltrigs'):
                            tcoll.fillTrackJetPair(track, tcuts, jet, jcuts[i])
                    elif 137221 in jcuts[i].trig:
                        if trig in ('137221', 'jetpatch', 'alltrigs'):
                            tcoll.fillTrackJetPair(track, tcuts, jet, jcuts[i])                                   
    
    
    def drawQA(self,):
        """creates and returns a canvas with a QA summary for the run"""
        c = ROOT.TCanvas('c', self.name, 600, 800)
        
        writer_config = { 'title':ROOT.TText(), 'text':ROOT.TText() }
        write = {}
        for key,value in writer_config.items():
            write[key] = value.DrawTextNDC
            
        pads = []
        
        writer_config['title'].SetTextSize(0.035)
        writer_config['title'].SetTextAlign(22)
        writer_config['text'].SetTextSize(0.02)
        
        margin = 0.02
        
        ## here's the title
        write['title'](0.5, 1.0 - margin, self.name)
        
        ## Basic run information
        #write['text'](margin, 0.95, '# events %d' % (self['anyspin'][]) )
        write['text'](margin, 0.85, 'minbias')
        
        #write['text'](0.2,0.2,'my text')
        t = ROOT.TText()
        t.DrawText(0.1,0.1,'this is so stupid')
        #pads.append( ROOT.TPad(str(uuid.uuid1()),'', margin, 0.74, 1.0-margin, 0.84) )
        pads.append( ROOT.TPad('some_random_name','', margin, 0.74, 1.0-margin, 0.84) )
        pad = pads[-1]
        pad.SetFillColor(2)
        pad.Divide(5,1)
        pad.Draw()
        pad.cd(1)
        ROOT.gPad.SetLogy()
        self['anyspin']['96011'].tracks_sum['pt'].Draw()
        
        
        #c.Draw()
        raw_input('wait here')
        return c
    
    
    def Write(self,):
        """make it persistent in the current tfile"""
        ## create charge-summed histos
        for spin in ('uu','ud','du','dd','other'):
            for trig in self.trigSetups:
                self[spin][trig].trackHistograms(0).Add( self[spin][trig].trackHistograms(1)    )
                self[spin][trig].trackHistograms(0).Add( self[spin][trig].trackHistograms(-1) )
        
        ## create spin-integrated histos
        for trig in self.trigSetups:
            self['anyspin'][trig].Add( self['uu'][trig] )
            self['anyspin'][trig].Add( self['ud'][trig] )
            self['anyspin'][trig].Add( self['du'][trig] )
            self['anyspin'][trig].Add( self['dd'][trig] )
            self['anyspin'][trig].Add( self['other'][trig] )
        
        ## write everything to the current TFile
        [ collection.Write() for collection in self.allHistos ]
    


def writeHistograms(treeDir='~/data/run5/tree', globber='*'):
    import analysis
    chain = ROOT.TChain('tree')
    chain.Add(treeDir + '/chargedPions_' + globber + '.tree.root')
    
    entries = chain.GetEntries()
    
    chain.GetEntry(0)
    fname = chain.GetCurrentFile().GetName()
    outname = os.path.basename(fname).replace('.tree.','.hist.')
    outFile = ROOT.TFile(outname, 'recreate')
    h = HistogramManager()
    h.fill = analysis.getFill(analysis.getRun(fname))
    
    for i in xrange(entries):
        chain.GetEntry(i)
        
        ## found a new runnumber
        if fname != chain.GetCurrentFile().GetName():
            outFile.cd()
            print 'saving', outname
            h.Write()
            outFile.Close()
            fname = chain.GetCurrentFile().GetName()
            outname = os.path.basename(fname).replace('.tree.','.hist.')
            outFile = ROOT.TFile(outname, 'recreate')
            h = HistogramManager()
            h.fill = analysis.getFill(analysis.getRun(fname))
            
        h.processEvent(chain.event)
    
    print 'saving', outname
    outFile.cd()
    h.Write()
    outFile.Close()


def condenseIntoFills(histDir='/Users/kocolosk/data/run5/hist',useLSF=False,fillList=None):
    """uses hadd to make histogram files for each fill instead of each run"""
    import analysis
    allFiles = os.listdir(histDir)
    fill_runlists = {}
    runlist = []
    for fname in allFiles:
        if not fname.endswith('.root'): continue
        run = analysis.getRun(fname)
        runlist.append(run)
    
    tuples = analysis.getAllFills(runlist)
    for run, fill in tuples:
        if run not in runlist: continue
        ## temporary hacks
        ## http://www.star.bnl.gov/HyperNews-star/protected/get/starspin/3324.html
        if 6144002 <= run <= 6144029: fill = 7128
        if 6144041 <= run <= 6144042: fill = 7129        
        if 6145067 <= run <= 6145068: fill = 7136
        if 6146001 <= run <= 6146026: fill = 7138
        try:
             fill_runlists[fill].append(run)
        except KeyError:
             fill_runlists[fill] = [run]
    
    for fill, runlist in fill_runlists.items():
        #if fill in (7048, 7055, 7327): continue  ## no final polarizations
        #if fill not in (7127, 7128, 7129, 7134, 7136, 7138): continue
        if fillList is None or fill in fillList:
            cmd = 'hadd chargedPions_%d.hist.root ' % (fill,)
            if len(runlist) == 1:
                print 'no need for hadd here!', fill, runlist[0]
                os.system('cp %s/chargedPions_%d.hist.root chargedPions_%d.hist.root' %
                    (histDir, runlist[0], fill))
            else:
                for run in runlist:
                    cmd += '%s/chargedPions_%d.hist.root ' % (histDir,run)
                if useLSF:
                    cmd = 'bsub -q star_cas_short -e err/%d.err -o out/%d.out ' % (fill,fill) + cmd
                os.system(cmd)

def bsub(treeDir, runlist=None):
    import analysis
    """submits a single writeHistograms job to LSF for each tree.root file in treeDir"""
    allfiles = os.listdir(treeDir)
    for fname in allfiles:
        if not fname.endswith('tree.root'): continue
        run = analysis.getRun(fname)
        if runlist is None or run in runlist:
            os.system('bsub -q star_cas_short -e err/%d.err -o out/%d.out python -c \
                "import analysis; analysis.histos.writeHistograms(\'%s\',globber=\'*%d*\')"' \
                % (run, run, treeDir, run))
            time.sleep(0.2)


if __name__ == '__main__':
    directory = sys.argv[1]
    run = sys.argv[2]
    writeHistograms(directory, run)
