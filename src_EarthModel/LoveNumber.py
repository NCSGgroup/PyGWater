import numpy as np
import scipy.io as scio

from src_auxiliary.FileTool import FileTool

from src_auxiliary.EnumClasses import LLN_Data, LLN_variable, Frame


class LoveNumber:
    """
    Load Love numbers (LLNs) h_n, l_n, k_n, all arrays start at degree 0 and have length lmax+1.

    Ref1: Wang et al. (2012) Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135,
          and modified models with refined crustal structure from Crust 2.0
    Ref2: Bevis et al. (2016) On computing the geoelastic response to a disk load
    Ref3: Blewitt (2003) Self-consistency in reference frames, geocenter definition, and surface loading of
          the solid Earth (degree-1 frame transformation)

    Usage (unchanged):
        lln = LoveNumber().config(lmax=1000, method=LLN_Data.PREM).get_Love_number()
        lln.convert(target=Frame.CF)

    Notes
    -----
    - Nothing is read from disk in __init__ any more; data are read by get_Love_number(), or lazily the first time
      `.LLN` is accessed. Previously __init__ always read Wang.dat (lmax=1000), so merely creating the object
      required the data file.
    - Degree-1 frame conversion is now computed from the degree-1 numbers of the LOADED Earth model
      (Blewitt, 2003), instead of hard-coded PREM values that were applied to every model and that differed from
      the CE values in the PREM file itself. Conversion works between any pair of CE / CM / CF.
    """

    '''file name for each Earth model in data/LLN: (file, header lines, max degree available)'''
    _complete_files = {
        LLN_Data.PREM: 'PREM-LLNs-complete.dat',
        LLN_Data.iasp91: 'iasp91-LLNs-complete.dat',
        LLN_Data.ak135: 'ak135-LLNs-complete.dat',
        LLN_Data.ak135hard: 'ak135hard-LLNs-complete.dat',
        LLN_Data.iasp91hard: 'iasp91hard-LLNs-complete.dat',
        LLN_Data.PREMhard: 'PREMhard-LLNs-complete.dat',
        LLN_Data.PREMsoft: 'PREMsoft-LLNs-complete.dat',
    }

    '''degree-0 vertical LLN (compressible PREM); the files start at degree 1'''
    h0 = -0.13205

    def __init__(self):
        self.__path = FileTool.get_data_dir('LLN')
        self.lmax = None
        self.method = None
        self.frame = None
        self._LLN = None
        self._deg1_CE = None
        self.config()
        pass

    def config(self, lmax=1000, method=LLN_Data.Wang):
        self.lmax = int(lmax)
        self.method = method
        self._LLN = None  # force re-reading with the new configuration
        return self

    @property
    def LLN(self):
        if self._LLN is None:
            self.get_Love_number()
        return self._LLN

    @LLN.setter
    def LLN(self, value):
        self._LLN = value

    def get_Love_number(self):
        """
        Read h, l, k up to self.lmax for self.method. Returns self.
        """
        if self.method == LLN_Data.Wang:
            self.__Wang()
        elif self.method == LLN_Data.REF:
            self.__REF()
        elif self.method in self._complete_files:
            self.__complete(self._complete_files[self.method])
        else:
            raise ValueError('Unknown Love number data set: %s' % self.method)

        '''remember the degree-1 numbers in the CE frame of this very Earth model, for frame conversion'''
        h1, l1, k1 = (self._LLN[v][1] for v in (LLN_variable.h, LLN_variable.l, LLN_variable.k))
        if self.frame == Frame.CE:
            self._deg1_CE = (h1, l1, k1)
        elif self.frame == Frame.CM:
            self._deg1_CE = (h1 + 1.0, l1 + 1.0, k1 + 1.0)
        else:
            raise NotImplementedError
        return self

    def __set(self, h, l, k, frame):
        self._LLN = {LLN_variable.k: np.asarray(k, dtype=float),
                     LLN_variable.h: np.asarray(h, dtype=float),
                     LLN_variable.l: np.asarray(l, dtype=float)}
        self.frame = frame

    def __Wang(self):
        """
        Wang.dat (no header), starting from degree 1; degree 0 is prepended. CE frame.
        """
        assert self.lmax <= 13000

        data = np.loadtxt(self.__path / 'Wang.dat', usecols=(1, 2, 3), max_rows=self.lmax)
        h, l, k = data[:self.lmax, 0], data[:self.lmax, 1], data[:self.lmax, 2]

        self.__set(h=np.append(self.h0, h), l=np.append(0, l), k=np.append(0, k), frame=Frame.CE)

    def __REF(self):
        """
        REF_6371 Love numbers (Bevis et al., 2016), starting from degree 0. CM frame.
        """
        assert self.lmax <= 40000

        love = scio.loadmat(self.__path / 'REF_6371_loading_love_numbers.mat')

        k = np.squeeze(love['k_love'])[:self.lmax + 1]
        h = np.squeeze(love['h_love'])[:self.lmax + 1]
        l = np.squeeze(love['l_love'])[:self.lmax + 1]

        self.__set(h=h, l=l, k=k, frame=Frame.CM)

    def __complete(self, filename):
        """
        '*-LLNs-complete.dat' of Wang et al. (2012): one header line, starting from degree 1; degree 0 is prepended.
        CE frame.
        """
        assert self.lmax <= 46000

        data = np.loadtxt(self.__path / filename, usecols=(1, 2, 3), max_rows=self.lmax, skiprows=1)
        h, l, k = data[:self.lmax, 0], data[:self.lmax, 1], data[:self.lmax, 2]

        self.__set(h=np.append(self.h0, h), l=np.append(0, l), k=np.append(0, k), frame=Frame.CE)

    def convert(self, target=Frame.CM):
        """
        Transform the degree-1 Love numbers into the target frame (Blewitt, 2003):
            CE -> CM : h1 - 1,  l1 - 1,  k1 - 1
            CE -> CF : h1 - a,  l1 - a,  k1 - a,  with a = (h1_CE + 2 l1_CE) / 3
        """
        LLN = self.LLN  # make sure data are loaded
        if self.frame == target:
            return self

        h1, l1, k1 = self._deg1_CE
        if target == Frame.CE:
            shift = 0.0
        elif target == Frame.CM:
            shift = -1.0
        elif target == Frame.CF:
            shift = -(h1 + 2 * l1) / 3.0
        else:
            raise ValueError(target)

        LLN[LLN_variable.h][1] = h1 + shift
        LLN[LLN_variable.l][1] = l1 + shift
        LLN[LLN_variable.k][1] = k1 + shift
        self.frame = target

        return self


def demo():
    lln = LoveNumber().config(lmax=300, method=LLN_Data.PREM).get_Love_number()
    for fr in (Frame.CM, Frame.CF, Frame.CE):
        lln.convert(fr)
        print(fr, [lln.LLN[v][1] for v in (LLN_variable.h, LLN_variable.l, LLN_variable.k)])
    pass


if __name__ == '__main__':
    demo()
