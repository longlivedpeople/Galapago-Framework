import math
import ROOT as r

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      self.initEventCuts()
      self.initEECuts()
      self.initMMCuts()


   #########################################
   #######  General cut definition   #######
   #########################################

   def initEventCuts(self):

      self.highPU = self.brackets('nPV > 22')
      self.medPU = self.brackets('nPV > 14 and nPV < 23')
      self.lowPU = self.brackets('nPV < 15')
      self.goodPV = self.brackets('PV_passAcceptance')


   ####################################
   #######  EE cut definition   #######
   ####################################

   def initEECuts(self): 

      self.EE_OnZ = self.brackets('abs(91 - ev.EE_mass[iee]) <= 10')
      self.EE_OffZ = self.brackets('abs(91 - ev.EE_mass[iee]) > 10')
      self.EE_OS = self.brackets('ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[iee]]]*ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[iee]]] < 0')
      self.EE_SS = self.brackets('ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxA[iee]]]*ev.IsoTrackSel_charge[ev.ElectronCandidate_isotrackIdx[ev.EE_idxB[iee]]] > 0')
      self.EE_iso0l = self.brackets('ev.EE_relisoA[iee] > 0.1 and ev.EE_relisoB[iee] > 0.1')
      self.EE_iso1l = self.brackets('(ev.EE_relisoA[iee] < 0.1 and ev.EE_relisoB[iee] > 0.1) or (ev.EE_relisoA[iee] > 0.1 and ev.EE_relisoB[iee] < 0.1)')
      self.EE_iso2l = self.brackets('ev.EE_relisoA[iee] < 0.1 and ev.EE_relisoB[iee] < 0.1')
      self.EE_leadingPt45 = self.brackets('ev.EE_leadingPt[iee] > 45')
      self.EE_subleadingPt28 = self.brackets('ev.EE_leadingPt[iee] > 28')
      self.EE_leadingEt45 = self.brackets('ev.EE_leadingEt[iee] > 45')
      self.EE_subleadingEt28 = self.brackets('ev.EE_leadingEt[iee] > 28')
      self.EE_normChi2_10 = self.brackets('ev.EE_normalizedChi2[iee] < 10')
      self.EE_mass15 = self.brackets('ev.EE_mass[iee] > 15')
      self.EE_eta2 = self.brackets('abs(ev.ElectronCandidate_eta[ev.EE_idxA[iee]]) < 2.0 and abs(ev.ElectronCandidate_eta[ev.EE_idxB[iee]]) < 2.0')
      self.EE_dPhibackward = self.brackets('ev.EE_dPhi[iee] > 3.14/2.0')
      self.EE_dPhiforward = self.brackets('ev.EE_dPhi[iee] < 3.14/2.0')
      self.EE_Ixy6high = self.brackets('ev.EE_trackIxy_PV[iee] > 6')
      self.EE_Ixy6prompt = self.brackets('ev.EE_trackIxy_PV[iee] < 6')


   ####################################
   #######  MM cut definition   #######
   ####################################

   def initMMCuts(self):

      self.MM_OnZ = self.brackets('abs(91 - ev.DMDM_mass[imm]) <= 10')
      self.MM_OffZ = self.brackets('abs(91 - ev.DMDM_mass[imm]) > 10')
      self.MM_OS = self.brackets('ev.DGM_charge[ev.DMDM_idxA[imm]]*ev.DGM_charge[ev.DMDM_idxB[imm]] < 0')
      self.MM_SS = self.brackets('ev.DGM_charge[ev.DMDM_idxA[imm]]*ev.DGM_charge[ev.DMDM_idxB[imm]] > 0')
      self.MM_iso0l = self.brackets('ev.DGM_relPFiso[ev.DMDM_idxA[imm]] > 0.2 and ev.DGM_relPFiso[ev.DMDM_idxB[imm]] > 0.2')
      self.MM_iso1l = self.brackets('(ev.DGM_relPFiso[ev.DMDM_idxA[imm]] < 0.2 and ev.DGM_relPFiso[ev.DMDM_idxB[imm]] > 0.2) or (ev.DGM_relPFiso[ev.DMDM_idxA[imm]] > 0.2 and ev.DGM_relPFiso[ev.DMDM_idxB[imm]] < 0.2)')
      self.MM_iso2l = self.brackets('ev.DGM_relPFiso[ev.DMDM_idxA[imm]] < 0.2 and ev.DGM_relPFiso[ev.DMDM_idxB[imm]] < 0.2')
      self.MM_dPhibackward = self.brackets('ev.DMDM_dPhi[imm] > 3.14/2.0')
      self.MM_dPhiforward = self.brackets('ev.DMDM_dPhi[imm] < 3.14/2.0')
      self.MM_Ixy6high = self.brackets('ev.DMDM_trackIxy_PV[imm] > 6')
      self.MM_Ixy6prompt = self.brackets('ev.DMDM_trackIxy_PV[imm] < 6')
      self.MM_mass15 = self.brackets('ev.DMDM_mass[imm] > 15')
      self.MM_normChi2_10 = self.brackets('ev.DMDM_normalizedChi2[imm] < 10')
      self.MM_cosAlpha0p8 = self.brackets('ev.DMDM_cosAlpha[imm] > -0.80')
      self.MM_dR0p2 = self.brackets('ev.DMDM_dR[imm] > 0.2')
      self.MM_ID = self.brackets('abs(ev.DGM_eta[ev.DMDM_idxA[imm]]) < 2.0 and abs(ev.DGM_eta[ev.DMDM_idxB[imm]]) < 2.0 and ev.DGM_pt[ev.DMDM_idxA[imm]] > 31 and ev.DGM_pt[ev.DMDM_idxB[imm]] > 31 and ev.DGM_ptError[ev.DMDM_idxB[imm]]/ev.DGM_pt[ev.DMDM_idxB[imm]] < 0.3 and ev.DGM_ptError[ev.DMDM_idxA[imm]]/ev.DGM_pt[ev.DMDM_idxA[imm]] < 0.3 and ev.DGM_ndof[ev.DMDM_idxA[imm]] > 0.00001 and ev.DGM_ndof[ev.DMDM_idxB[imm]] > 0.00001 and ev.DGM_chi2[ev.DMDM_idxA[imm]]/ev.DGM_ndof[ev.DMDM_idxA[imm]] < 10 and ev.DGM_chi2[ev.DMDM_idxB[imm]]/ev.DGM_ndof[ev.DMDM_idxB[imm]] < 10 and ev.DGM_numberOfValidHits[ev.DMDM_idxA[imm]] > 22 and ev.DGM_numberOfValidHits[ev.DMDM_idxB[imm]] > 22')
      


   ###########################################
   #######  Logical cut combinations   #######
   ###########################################

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


   #################################
   #######  Cut evaluation   #######
   #################################
      
   def passCut(self, event, cut, index = False):

      if type(index) != bool:
          cut = cut.format(str(index))

      formula = r.TTreeFormula("_auxCut", cut, event)
      formula.GetNdata()
      passed = formula.EvalInstance()

      return passed 






   ###########################################
   #######  Logical cut combinations   #######
   ###########################################

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


   #################################
   #######  Cut evaluation   #######
   #################################
      
   def passCut(self, event, cut, index = False):

      if type(index) != bool:
          cut = cut.format(str(index))

      formula = r.TTreeFormula("_auxCut", cut, event)
      formula.GetNdata()
      passed = formula.EvalInstance()

      return passed 




