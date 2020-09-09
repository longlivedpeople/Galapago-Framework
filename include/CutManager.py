import math

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      ############################
      ######  Event cuts   #######
      ############################
      self.nTrack = self.brackets('RefittedPV_nPFTrack + RefittedPV_nLostTrack + RefittedPV_nExcludedTrack > 5')
      self.highPU = self.brackets('nPUTrue > 35')
      self.lowPU = self.brackets('nPUTrue < 20')


      #################################
      ######  Basic Lepton Cuts  ######
      #################################
      self.EEChannel = self.brackets('Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15 == 1 && Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 == 0')
      self.haveEEBase = self.brackets('nEEBase > 0')
      self.EESR_dPhi = self.brackets('fabs(EEBase_dPhi[EEBase_maxIxy])< 3.14/2.0')
      self.EECR_dPhi = self.brackets('fabs(EEBase_dPhi[EEBase_maxIxy]) > 3.14/2.0')
      self.EE_highpT = self.brackets('EEBase_leadingPt[EEBase_maxIxy] > 50 && EEBase_subleadingPt[EEBase_maxIxy] > 40')
      self.EE_highET = self.brackets('EEBase_leadingEt[EEBase_maxIxy] > 50 && EEBase_subleadingEt[EEBase_maxIxy] > 40')
      self.EEtailRegime = self.brackets('fabs(EEBase_trackIxy[EEBase_maxIxy])> 5.0')
      self.EEpromptRegime = self.brackets('fabs(EEBase_trackIxy[EEBase_maxIxy])< 5.0')


      self.MMChannel = self.brackets('Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 == 1')
      self.haveMM = self.brackets('nDMDMBase > 0')
      self.MMSR_dPhi = self.brackets('fabs(DMDMBase_dPhi[DMDMBase_maxIxy])< 3.14/2.0')
      self.MMCR_dPhi = self.brackets('fabs(DMDMBase_dPhi[DMDMBase_maxIxy])> 3.14/2.0')
      self.MM_highpT = self.brackets('DMDMBase_leadingPt[DMDMBase_maxIxy] > 40 && DMDMBase_subleadingPt[DMDMBase_maxIxy] > 40')
      self.MM_SScharge = self.brackets('DGM_charge[DMDMBase_idxA[DMDMBase_maxIxy]]*DGM_charge[DMDMBase_idxB[DMDMBase_maxIxy]] > 0')
      self.MM_OScharge = self.brackets('DGM_charge[DMDMBase_idxA[DMDMBase_maxIxy]]*DGM_charge[DMDMBase_idxB[DMDMBase_maxIxy]] < 0')
      self.cosmicRejection = self.brackets('DMDMBase_cosAlpha[DMDMBase_maxIxy] > -0.80')

      #########################
      ######  Loop Cuts  ######
      #########################
      self.EEChannel = self.brackets('Flag_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15 == 1 && Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 == 0')
      self.LoopHaveEEBase = self.brackets('ev.nEEBase > 0')
      self.LoopEESR_dPhi = self.brackets('abs(ev.EEBase_dPhi[ee_maxIxy])< 3.14/2.0')
      self.LoopEECR_dPhi = self.brackets('abs(ev.EEBase_dPhi[ee_maxIxy]) > 3.14/2.0')
      self.LoopEE_pT = self.brackets('ev.EEBase_leadingPt[ee_maxIxy] > 50 && ev.EEBase_subleadingPt[ee_maxIxy] > 40')
      self.LoopEE_ET = self.brackets('ev.EEBase_leadingEt[ee_maxIxy] > 50 && ev.EEBase_subleadingEt[ee_maxIxy] > 40')
      self.LoopEE_SScharge = self.brackets('ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxA[ee_maxIxy]]]*ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxB[ee_maxIxy]]] > 0')
      self.LoopEE_OScharge = self.brackets('ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxA[ee_maxIxy]]]*ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EEBase_idxB[ee_maxIxy]]] < 0')


      self.LoopMMChannel = self.brackets('ev.Flag_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10 == 1')
      self.LoopHaveMM = self.brackets('ev.nDMDMBase > 0')
      self.LoopMMSR_dPhi = self.brackets('abs(ev.DMDMBase_dPhi[mm_maxIxy])< 3.14/2.0')
      self.LoopMMCR_dPhi = self.brackets('abs(ev.DMDMBase_dPhi[mm_maxIxy])> 3.14/2.0')
      self.LoopMM_SScharge = self.brackets('ev.DGM_charge[ev.DMDMBase_idxA[mm_maxIxy]]*ev.DGM_charge[ev.DMDMBase_idxB[mm_maxIxy]] > 0')
      self.LoopMM_OScharge = self.brackets('ev.DGM_charge[ev.DMDMBase_idxA[mm_maxIxy]]*ev.DGM_charge[ev.DMDMBase_idxB[mm_maxIxy]] < 0')
      self.LoopCosmicRejection = self.brackets('ev.DMDMBase_cosAlpha[mm_maxIxy] > -0.80')

      


   def donotB(self, cut):
     return '(!' + self.brackets(cut) + ')'

   def brackets(self, cut):
      if cut != '':
          return '('+cut+')'

   def AddListB(self, cutlist):
      returncut = ''
      for cut in cutlist:
          if cut != '':
              returncut += cut
              if not cutlist.index(cut) == len(cutlist)-1:
                  returncut += ' && '
      return self.brackets(returncut)
  
   def AddB(self, cut1, cut2):
      if cut1 == '':
          return cut2 
      if cut2 == '':
          return cut1
      return self.brackets(cut1 + " && " + cut2 )
  
   def ORB(self, cut1, cut2):

      return self.brackets(cut1 + " || " + cut2 )

   def donot(self, cut):
     return '(not ' + self.brackets(cut) + ')'

   def AddList(self, cutlist):
      returncut = ''
      for cut in cutlist:
          if cut != '':
              returncut += cut
              if not cutlist.index(cut) == len(cutlist)-1:
                  returncut += ' and '
      return self.brackets(returncut)
  
   def Add(self, cut1, cut2):
      if cut1 == '':
          return cut2 
      if cut2 == '':
          return cut1
      return self.brackets(cut1 + " and " + cut2 )
  
   def OR(self, cut1, cut2):

      return self.brackets(cut1 + " or " + cut2 )


      
