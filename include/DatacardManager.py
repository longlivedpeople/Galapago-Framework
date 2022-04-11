import ROOT as r
import os




class Recipe:

    """ Class Recipe """

    def __init__(self, recipe):

        self.dcs = {}

        _file = open(recipe, 'r')
        _lines = _file.readlines()

        for _l in _lines:

            if _l[0] == '#': continue
            elems = _l.split(' ')
            elems = [x for x in elems if x!= '']

            name      = elems[0]
            channel   = elems[1]
            ctype     = elems[2]
            histo     = elems[3]
            limits    = elems[4]
            directory = elems[5] if '\n' not in elems[5] else elems[5][:-1]

            if name not in self.dcs.keys():
                self.dcs[name] = {}

            if channel not in self.dcs[name].keys():
                self.dcs[name][channel] = self.initChannel()

            self.dcs[name][channel][ctype]['histogram'] = histo
            self.dcs[name][channel][ctype]['limits'][0] = limits.split(',')[0]
            self.dcs[name][channel][ctype]['limits'][1] = limits.split(',')[1]
            self.dcs[name][channel][ctype]['dir'] = directory
                    



    def initChannel(self):

        channel = {}

        channel['bkg'] = {'histogram' : '', 'limits' : [0,0], 'dir' : ''}
        channel['sig'] = {'histogram' : '', 'limits' : [0,0], 'dir' : ''}

        return channel


class Datacard:

    """ Class Datacard """

    def __init__(self, name):

        self.channels = []
        self.channelCounter = 0
        self.width = 15 # Normal column width
        self.hwidth = 20 # First column width
        self.signalName = ''
        self.nBackground = 0
        self.nNuisance = 4 # Provisional: TO BE UPDATED
        self.heading = ''
        self.observation = ''
        self.expected = ''
        self.separator = self.hwidth*'-'
        self.name = name

    def addChannel(self, channel):

        if self.channelCounter == 0:
            self.signalName = channel.signal.name
            self.nBackground = channel.backgroundCounter
        else:
            if self.signalName != channel.signal.name: 
                raise ValueError('Invalid datacard: All channels need to have the same signal')
            if self.nBackground != channel.backgroundCounter:
                raise ValueError('Invalid datacard: All channel need to have the same number of backgrounds')

        self.channels.append(channel)
        self.separator += (channel.backgroundCounter + 1)*self.width*'-'
        self.channelCounter+=1


    def toCell(self, string, width = False):

        """ Function that receives a string and returns a modified string
         with the length of the datacard cell. """ 
  
        nSpace = self.width
        if width: nSpace = width

        cell = string + (nSpace - len(string))*' '

        return cell


    def saveDatacard(self, outputDir = False):

        self.separator += '\n' # we suppose we have added all the backgrounds we wanted

        # heading definition
        self.heading = 'imax ' + str(self.channelCounter) + ' number of bins'
        self.heading += '\n'
        self.heading += 'jmax ' + str(self.nBackground) + ' number of processes minus 1'
        self.heading += '\n'
        self.heading += 'kmax ' + str(self.nNuisance) + ' number of nuisance parameters'       
        self.heading += '\n'

        # observation definition
        self.observation = self.toCell('bin', self.hwidth)
        for channel in self.channels:
            self.observation += self.toCell(channel.name)
        self.observation += '\n'

        self.observation = self.toCell('observation', self.hwidth)
        for channel in self.channels:
            self.observation += self.toCell('{:.0f}'.format(channel.data))
        self.observation += '\n'

        # expected definition
        self.expected = self.toCell('bin', self.hwidth)
        for channel in self.channels:
            self.expected += (channel.backgroundCounter + 1)*(self.toCell(channel.name))
        self.expected += '\n'

        self.expected += self.toCell('process', self.hwidth)
        for channel in self.channels:
            self.expected += self.toCell('signal')
            for background in channel.backgrounds:
                self.expected += self.toCell(background.name)
        self.expected += '\n'

        self.expected += self.toCell('process', self.hwidth)
        for channel in self.channels:
            self.expected += self.toCell(channel.signal.number)
            for background in channel.backgrounds:
                self.expected += self.toCell(str(background.number))
        self.expected += '\n'

        self.expected += self.toCell('rate', self.hwidth)
        for channel in self.channels:
            self.expected += self.toCell('{:.4f}'.format(channel.signal.rate))
            for background in channel.backgrounds:
                self.expected += self.toCell('{:.4f}'.format(background.rate))
        self.expected += '\n'

        # nuisances
        # (por el momeno solo podemos poner constantes, TO BE UPDATED)
        self.nuisance = ''
        self.nuisance += self.toCell('lumi', int(self.hwidth/2))
        self.nuisance += self.toCell('lnN', int(self.hwidth/2))
        for i in range(0, self.channelCounter*(self.nBackground + 1)):
            if ((i%2) == 0):
                self.nuisance += self.toCell('1.025')
            else:
                self.nuisance += self.toCell('-')
        self.nuisance += '\n'

        self.nuisance += self.toCell('PU', int(self.hwidth/2))
        self.nuisance += self.toCell('lnN', int(self.hwidth/2))
        for i in range(0, self.channelCounter*(self.nBackground + 1)):
            if ((i%2) == 0):
                self.nuisance += self.toCell('0.975/1.030')
            else:
                self.nuisance += self.toCell('-')
        self.nuisance += '\n'
 
        self.nuisance += self.toCell('recoID', int(self.hwidth/2))
        self.nuisance += self.toCell('lnN', int(self.hwidth/2))
        for i in range(0, self.channelCounter*(self.nBackground + 1)):
            if ((i%2) == 0):
                self.nuisance += self.toCell('1.10')
            else:
                self.nuisance += self.toCell('-')
        self.nuisance += '\n'

        self.nuisance += self.toCell('dphi', int(self.hwidth/2))
        self.nuisance += self.toCell('lnN', int(self.hwidth/2))
        for i in range(0, self.channelCounter*(self.nBackground + 1)):
            if ((i%2) == 0):
                self.nuisance += self.toCell('-')
            else:
                self.nuisance += self.toCell('1.10')

        ### Write file:

        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        outputPath = outputDir if outputDir else ''

        _f = open(outputDir + self.name, 'w')
        _f.write(self.heading)
        _f.write(self.separator)
        _f.write(self.observation)
        _f.write(self.separator)
        _f.write(self.expected)
        _f.write(self.separator)
        _f.write(self.nuisance)
        _f.close()


class Channel:

    """ Class Channel """

    def __init__(self, name, xmin, xmax):

        self.name = name # ee o mm
        self.backgrounds = []
        self.backgroundCounter = 0
        self.signal = 0
        self.data = 0
        self.xmin = xmin
        self.xmax = xmax


    def addBackground(self, sample, h):

        self.backgroundCounter += 1 
  
        name = sample

        nbinmin = h.FindBin(float(self.xmin)) 
        if self.xmax == 'inf':
            nbinmax = h.GetNbinsX()
        else:
            nbinmax = h.FindBin(float(self.xmax))
        value = h.Integral(nbinmin, nbinmax)

        # add the background
        background = EventSet(name = name, number = self.backgroundCounter, rate = value)
        self.backgrounds.append(background)


    def setSignal(self, sample, h):

        name = sample

        nbinmin = h.FindBin(float(self.xmin))
        if self.xmax == 'inf':
            nbinmax = h.GetNbinsX()
        else:
            nbinmax = h.FindBin(float(self.xmax))
        value = h.Integral(nbinmin, nbinmax)

        signal = EventSet(name = name, number = '0', rate = value)
        self.signal = signal



    def addData(self, h):

        nbinmin = h.FindBin(self.xmin)
        nbinmax = h.FindBin(self.xmax)
        value = h.Integral(nbinmin, nbinmax)

        self.data += value


    


class EventSet:

    def __init__(self, name, number, rate):

        self.name = name
        self.number = number
        self.rate = rate





