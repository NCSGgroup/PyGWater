import warnings

import numpy as np

from src_auxiliary.EnumClasses import PhysicalDimensions
from src_auxiliary.MathTool import MathTool


class Harmonic:
    """
    Harmonic analysis and synthesis: Ordinary 2D integration for computing Spherical Harmonic coefficients
    """

    def __init__(self, lat, lon, lmax: int, option=0):
        """

        :param lat: Co-latitude if option=0, unit[rad]; geophysical latitude if option = others, unit[degree]
        :param lon: If option=0, unit[rad]; else unit[degree]
        :param lmax: int, max degree/order
        :param option:
        """
        if option != 0:
            self.lat, self.lon = MathTool.get_colat_lon_rad(lat, lon)
        else:
            self.lat, self.lon = lat, lon

        self.lmax = lmax

        self._prepare()

    def _prepare(self):
        """pre-process some parameters of the 'two-step' method."""
        self.nlat, self.nlon = len(self.lat), len(self.lon)

        self.pilm = MathTool.get_Legendre(self.lat, self.lmax)
        # Associative Legendre polynomials, indexes stand for (co-lat[rad], degree l, order m)
        # 3-d array [theta, l, m] shape: (nlat * (lmax + 1) * (lmax + 1))

        # self.si = np.sin(self.lat) * np.pi / (2 * self.lmax + 1)

        # get Neumann weights
        # x_mat = np.array([np.cos(self.lat) ** i for i in range(self.nlat)])
        # r_vec = np.ones(self.nlat) * 2 / np.arange(1, self.nlat + 1, 1)
        # r_vec[np.arange(1, self.nlat + 1, 2)] = 0
        # self.wi = np.linalg.pinv(x_mat) @ r_vec

        m = np.arange(self.lmax + 1)
        self.g = m[:, None] @ self.lon[None, :]

        self.factor1 = np.ones((self.nlat, self.lmax + 1))
        self.factor1[:, 0] += 1
        '''Nyquist order m = nlon/2 (only reachable when lmax >= nlon/2): on nlon equally spaced longitudes
        sum(cos^2) or sum(sin^2) equals nlon instead of nlon/2, so this order must be halved as well.
        This is the delta_{m,Nmax} term of Eq. 5 in Yang et al. (2025). Without it, on a cell-centred grid the
        S_{n,Nmax} coefficients were doubled (and C_{n,Nmax} are not observable at all).'''
        self.m_nyquist = self.nlon / 2 if self.nlon % 2 == 0 else None
        if self.m_nyquist is not None and self.m_nyquist <= self.lmax:
            self.factor1[:, int(self.m_nyquist)] += 1
        self.factor1 = 1 / (self.factor1 * self.nlon)

        self._pinv = None  # per-order pseudo-inverses for the least-squares analysis (built on demand)

        '''Driscoll-Healy (type I) latitude weights for the cell-centred equi-angular grid, as in SaGEA
        (sagea.harmonics.harmonic.Harmonic, custom regular grid):
            w_j ~ sin(theta_j) * 4/pi * sum_{l=0}^{nlat/2-1} sin((2l+1) theta_j) / (2l+1),   sum_j w_j = 1
        The series runs over nlat/2 terms (standard DH); SaGEA 0.3.3 stops one term earlier (range(nlat/2-1)),
        which loses exactness at the last admissible degree. With these weights the analysis is exact for
        band-limited fields with lmax <= nlat/2 - 1 (e.g. 44 on a 2-degree grid, 89 on a 1-degree grid).'''
        wi = np.array([np.sum(np.sin((2 * np.arange(self.nlat // 2) + 1) * t) / (2 * np.arange(self.nlat // 2) + 1))
                       for t in self.lat]) * 4 / np.pi
        self.w_dh = np.sin(self.lat) * wi
        self.w_dh /= np.sum(self.w_dh)

        self.factor2 = np.ones((self.lmax + 1, self.lmax + 1))
        self.factor2[:, 0] += 1
        self.factor2 *= np.pi / (2 * self.nlat)

        self.factor3 = np.ones((self.lmax + 1, self.lmax + 1))
        self.factor3[:, 0] += 1
        self.factor3 /= 2
        pass

    def _get_pinv(self):
        """pseudo-inverse of the latitude-dependent design matrix, one per order m (least-squares analysis)"""
        if self._pinv is None:
            self._pinv = [np.linalg.pinv(self.pilm[:, m:, m]) for m in range(self.lmax + 1)]
        return self._pinv

    def analysis(self, gqij: np.ndarray, special_type: PhysicalDimensions = None, method='auto'):
        """
        :param gqij: gridded field(s), shape (nlat, nlon) or (ntime, nlat, nlon)
        :param special_type: None only (horizontal displacement analysis is not supported)
        :param method:
            'auto'          : (default) 'dh' if lmax <= nlat/2 - 1 (where it is exact), otherwise 'midpoint'
                              (which keeps the finite-LGF / spherical-harmonic equivalence of Sect. 4.2 of
                              Yang et al., 2025, for block-mean grids analysed up to Nmax = pi/dlambda).
            'dh'            : SaGEA scheme: Fourier sums in longitude, Driscoll-Healy latitude weights.
                              Exact for band-limited fields if lmax <= nlat/2 - 1 (sampling theorem); beyond
                              that (e.g. lmax = 90 on a 2-degree grid, as Eq. 6 of Yang et al. 2025 suggests)
                              no latitude quadrature on this grid is exact and a warning is issued.
            'midpoint'      : (alias 'quadrature') the original PyGWater scheme used for the figures of
                              Yang et al. (2025): two-step quadrature of Sneeuw (1994), Eqs. 3-5.
                              Grid values are treated as block means. The latitude integral uses the midpoint
                              rule, which is NOT exactly orthogonal on an equi-angular grid: for a band-limited
                              field the zonal (m=0) coefficients are biased by ~0.1% at n=10 and ~1% at n=60.
            'least_squares' : longitude part as above, latitude part solved per order m by least squares. This
                              reproduces any band-limited field exactly when lmax < nlat (e.g. lmax <= 89 on a
                              2-degree grid); for lmax = nlat the m=0 system is under-determined and the
                              minimum-norm solution is returned.
        """
        assert len(gqij.shape) in (2, 3)
        if method == 'quadrature':
            method = 'midpoint'
        if method == 'auto':
            method = 'dh' if self.lmax <= self.nlat // 2 - 1 else 'midpoint'
        assert method in ('dh', 'midpoint', 'least_squares')
        if method == 'dh' and self.lmax > self.nlat // 2 - 1:
            warnings.warn('DH analysis: lmax=%d exceeds nlat/2-1=%d; the latitude quadrature is not exact '
                          '(use a finer grid, a lower lmax, or method="least_squares").'
                          % (self.lmax, self.nlat // 2 - 1), stacklevel=2)

        single = (len(gqij.shape) == 2)
        if single:
            gqij = np.array([gqij])

        assert special_type in (
            PhysicalDimensions.HorizontalDisplacementEast, PhysicalDimensions.HorizontalDisplacementNorth, None)

        if special_type in (
                PhysicalDimensions.HorizontalDisplacementEast, PhysicalDimensions.HorizontalDisplacementNorth):
            assert False, "Horizontal Displacement is not supported yet."

        g = self.g
        co = np.cos(g)  # cos(m phi)
        so = np.sin(g)  # sin(m phi)

        am = np.einsum('pij,mj->pim', gqij, co, optimize='greedy') * self.factor1
        bm = np.einsum('pij,mj->pim', gqij, so, optimize='greedy') * self.factor1

        if special_type is None and method == 'least_squares':
            '''Fourier coefficients a_m(theta), b_m(theta) such that g = sum_m a_m cos(m lon) + b_m sin(m lon)'''
            am_f, bm_f = 2 * am, 2 * bm
            cqlm = np.zeros((am.shape[0], self.lmax + 1, self.lmax + 1))
            sqlm = np.zeros((am.shape[0], self.lmax + 1, self.lmax + 1))
            for m, pinv in enumerate(self._get_pinv()):
                cqlm[:, m:, m] = np.einsum('ni,pi->pn', pinv, am_f[:, :, m])
                sqlm[:, m:, m] = np.einsum('ni,pi->pn', pinv, bm_f[:, :, m])
            sqlm[:, :, 0] = 0.

        elif special_type is None and method == 'dh':
            '''1/nlon Fourier sums for every order (factor1 * (1 + delta_m0)), weights summing to one'''
            scale = np.ones(self.lmax + 1)
            scale[0] = 2.
            cqlm = np.einsum('pim,ilm,i->plm', am * scale, self.pilm, self.w_dh, optimize='greedy')
            sqlm = np.einsum('pim,ilm,i->plm', bm * scale, self.pilm, self.w_dh, optimize='greedy')

        elif special_type is None:
            cqlm = np.einsum('pim,ilm,i->plm', am, self.pilm, np.sin(self.lat), optimize='greedy') * self.factor2
            sqlm = np.einsum('pim,ilm,i->plm', bm, self.pilm, np.sin(self.lat), optimize='greedy') * self.factor2

            # cqlm = np.einsum('pim,ilm,i->plm', am, self.pilm, self.wi, optimize='greedy') * self.factor3
            # sqlm = np.einsum('pim,ilm,i->plm', bm, self.pilm, self.wi, optimize='greedy') * self.factor3

            # cqlm[:, :, 0] = (np.einsum('pim,ilm->plm', am * self.wi[:, None], self.pilm,
            #                               optimize='greedy') * self.factor3)[:, :, 0]
            # sqlm[:, :, 0] = (np.einsum('pim,ilm->plm', bm * self.wi[:, None], self.pilm,
            #                            optimize='greedy') * self.factor3)[:, :, 0]

            # cqlm = np.einsum('pim,ilm->plm', am * self.wi[:, None], self.pilm,
            #                  optimize='greedy') * self.factor3
            # sqlm = np.einsum('pim,ilm->plm', bm * self.wi[:, None], self.pilm,
            #                  optimize='greedy') * self.factor3

        elif special_type == PhysicalDimensions.HorizontalDisplacementNorth:
            """do NOT use this code!"""

            pilm_derivative = MathTool.get_Legendre_derivative(self.lat, self.lmax)

            cqlm = np.einsum('pim,ilm,i->plm', -am, pilm_derivative, np.sin(self.lat)) * self.factor2
            sqlm = np.einsum('pim,ilm,i->plm', -bm, pilm_derivative, np.sin(self.lat)) * self.factor2

        elif special_type == PhysicalDimensions.HorizontalDisplacementEast:
            """do NOT use this code!"""

            mrange = np.arange(self.lmax + 1)
            pilm_divide_sin_theta = self.pilm / np.sin(self.lat)[:, None, None]
            m_pilm_divide_sin_theta = np.einsum("m,ilm->ilm", mrange, pilm_divide_sin_theta)

            am_east = np.einsum('pij,mj->pim', gqij, -so, optimize='greedy') * self.factor1
            bm_east = np.einsum('pij,mj->pim', gqij, co, optimize='greedy') * self.factor1

            cqlm = np.einsum('pim,ilm,i->plm', am_east, m_pilm_divide_sin_theta, np.sin(self.lat)) * self.factor2
            sqlm = np.einsum('pim,ilm,i->plm', bm_east, m_pilm_divide_sin_theta, np.sin(self.lat)) * self.factor2

        else:
            assert False

        if single:
            assert cqlm.shape[0] == 1 and sqlm.shape[0] == 1
            return cqlm[0], sqlm[0]

        else:
            return cqlm, sqlm

    def synthesis(self, cqlm: np.ndarray, sqlm: np.ndarray, special_type: PhysicalDimensions = None):
        assert cqlm.shape == sqlm.shape
        assert len(cqlm.shape) in (2, 3)
        single = (len(cqlm.shape) == 2)
        if single:
            cqlm = np.array([cqlm])
            sqlm = np.array([sqlm])

        assert special_type in (
            PhysicalDimensions.HorizontalDisplacementEast, PhysicalDimensions.HorizontalDisplacementNorth, None)

        cqlm = np.array(cqlm)
        sqlm = np.array(sqlm)

        if special_type is None:
            am = np.einsum('ijk,ljk->ilk', cqlm, self.pilm)
            bm = np.einsum('ijk,ljk->ilk', sqlm, self.pilm)

        elif special_type is PhysicalDimensions.HorizontalDisplacementNorth:
            pilm_derivative = MathTool.get_Legendre_derivative(self.lat, self.lmax)

            am = np.einsum('ijk,ljk->ilk', -cqlm, pilm_derivative)
            bm = np.einsum('ijk,ljk->ilk', -sqlm, pilm_derivative)

        elif special_type is PhysicalDimensions.HorizontalDisplacementEast:
            pilm_divide_sin_theta = self.pilm / np.sin(self.lat)[:, None, None]
            mrange = np.arange(self.lmax + 1)

            am = np.einsum('ijk,k,ljk->ilk', sqlm, mrange, pilm_divide_sin_theta)
            bm = np.einsum('ijk,k,ljk->ilk', -cqlm, mrange, pilm_divide_sin_theta)

        else:
            assert False

        co = np.cos(self.g)
        so = np.sin(self.g)

        gqij = am @ co + bm @ so

        if single:
            assert gqij.shape[0] == 1
            return gqij[0]
        else:
            return gqij
