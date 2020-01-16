import math

class CutManager:
   'This class serves as an on-demand cut server'

   def __init__(self):

      ########################################################################
      ######Basic Lepton Cuts ################################################
      ########################################################################
      self.haveEEBase = self.brackets('nEEBase > 0')
      self.EESR_dPhi = self.brackets('fabs(EEBase_dPhi[EEBase_maxIxy])< 3.14/2.0')
      self.EECR_dPhi = self.brackets('fabs(EEBase_dPhi[EEBase_maxIxy]) > 3.14/2.0')
      self.haveMMBase = self.brackets('nMMBase > 0')

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


      
