import pandas as pd
from sim import Sim
from __init__ import ureg,Q_
import numpy
from laminate import Laminate


def quadractic_roots(A, B, C):
    tol = 10 ** -14

    def calc(the_input):
        a, b, c = the_input
        if abs(a) < tol and abs(b) < tol:
            return [0.0, 0.0]
        elif abs(a) < tol:
            root = -c / b
            signs = [root > 0.0, root < 0.0]
            return [root * cond for cond in signs]

        delta = b ** 2 - 4 * a * c
        if delta < 0.0:
            message = "No real roots for %r, %r, %r \n" % (a, b, c)
            message += "Delta : %r" % delta
            raise ValueError(message)

        roots = [(-1.0 * b + x * numpy.sqrt(delta)) / (2.0 * a) for x in (1, -1)]
        pos_roots = max(roots);
        neg_roots = min(roots)

        return [pos_roots, neg_roots]

    R = numpy.array(map(calc, zip(A, B, C)))
    return R

def max_R(input_stress, properties):
    on_stress = numpy.vstack(input_stress)
    xt = numpy.vectorize(lambda x: properties['xt'] / x if x > 0 else 0.0)
    xc = numpy.vectorize(lambda x: numpy.abs(properties['xc'] / x) if x < 0 else 0.0)
    yt = numpy.vectorize(lambda x: properties['yt'] / x if x > 0 else 0.0)
    yc = numpy.vectorize(lambda x: numpy.abs(properties['yc'] / x) if x < 0 else 0.0)
    s = numpy.vectorize(lambda x: numpy.abs(properties['sc'] / x) if numpy.abs(x) > 10 ** -10 else 0.0)
    R = numpy.vstack((xt(on_stress[:, 0]),
                      xc(on_stress[:, 0]),
                      yt(on_stress[:, 1]),
                      yc(on_stress[:, 1]),
                      s(on_stress[:, 2])
    )).T
    return R

def quad_R(input_stress, properties):
    # Unpack stress matrix into three stress vectors
    sigma_x, sigma_y, sigma_s = numpy.vstack(input_stress).T[:]

    Fxx = 1.0 / (properties['xt'] * properties['xc'])
    Fx = 1.0 / properties['xt'] - 1.0 / properties['xc']
    Fyy = 1.0 / (properties['yt'] * properties['yc'])
    Fy = 1.0 / properties['yt'] - 1.0 / properties['yc']
    Fss = 1.0 / (properties['sc'] ** 2)
    Fxy_star = -1.0 / 2.0;
    Fxy = Fxy_star * numpy.sqrt(Fxx * Fyy)

    A = Fxx * sigma_x ** 2 + 2 * Fxy * sigma_x * sigma_y + Fyy * sigma_y ** 2 + Fss * sigma_s ** 2
    B = Fx * sigma_x + Fy * sigma_y
    C = -1.0 * numpy.ones_like(A)

    num_rows = len(sigma_x)
    return numpy.array(quadractic_roots(A, B, C)).reshape(num_rows, 2)


def hashin_R(input_stress, properties):
    sigma_x, sigma_y, sigma_s = numpy.vstack(input_stress).T[:]
    xt = numpy.vectorize(lambda x, s: \
            numpy.sqrt(((x / properties['xt']) ** 2
                         + (s / properties['sc']) ** 2) ** -1
    ) if x > 0 else 0.0)
    xc = numpy.vectorize(lambda x: properties['xc'] / numpy.abs(x) \
        if x < 0 else 0.0)
    yt = numpy.vectorize(lambda y, s: \
             numpy.sqrt(((y / properties['yt']) ** 2
                            + (s / properties['sc']) ** 2) ** -1
    ) if y > 0 else 0.0)
    # For yc we need to solve a quadratic.
    # We find roots for all and then multiply by a boolean to eliminate the y > 0
    A = (1.0 / 4.0) * (sigma_y / properties['sc'])**2 \
    					+ (sigma_s/properties['sc'])**2
    B = (1.0/4.0* (properties['yc']/properties['sc'])**2 - 1) * sigma_y \
    				/ properties['yc']
    C = -1.0 * numpy.ones_like(A)
    yc_roots = quadractic_roots(A, B, C)[:, 0]
    condition = numpy.array(sigma_y < 0, dtype=int)
    yc = yc_roots * condition

    R = numpy.vstack((xt(sigma_x, sigma_s),
                      xc(sigma_x),
                      yt(sigma_y, sigma_s),
                      yc)).T

    return R

def find_R_mins(dataframe,columns):
    lowest_Rs = []; plys = []; 
    modes = [x.split('_')[0] for x in columns]
    nonzero_dat = dataframe[dataframe != 0]
    for col in columns:
        the_data = nonzero_dat[col]
        min_index = the_data.idxmin()
        try:
            lowest_Rs.append(the_data.iloc[min_index])
            plys.append(dataframe['Ply'].iloc[min_index])
        except ValueError:
            lowest_Rs.append(numpy.nan)
            plys.append('N/A')
    return pd.DataFrame({'Mode':modes,
            'Lowest R':lowest_Rs,
            'Ply':plys
            })[['Mode','Lowest R','Ply']]


class FailureAnalysis(object):
    """Begin failure analysis on a solved simulation object.
    Strenghts are converted to GPa within initialization
    """

    table_made = False

    def __init__(self, input_sim):
        try:
            self.stress = input_sim.get_stress()
        except AttributeError:
            message = ("Can't obtain stresses from input sim : %r -- type : %s"
                       % (sim, type(sim)))
            raise AttributeError(message)

        except input_sim.WorkflowError:
            raise WorkflowError("Could not get stresses from input sim")

        PROPS = input_sim.laminate.layers[0].PROPS
        required = ['xt', 'xc', 'yc', 'yt', 'sc']
        self.strength = {i: PROPS[i] / 1000.0 for i in required}  # Strength in MPa
        # self.strength = {i:PROPS[i] for i in required}
        self.sim = input_sim

    def compute_max(self):
        self.max_R = max_R(self.stress, self.strength)

    def compute_quad(self):
        self.quad_R = quad_R(self.stress, self.strength)

    def compute_hashin(self):
        self.hashin_R = hashin_R(self.stress, self.strength)

    def get_max(self):
        try:
            return self.max_R
        except NameError:
            self.compute_max()
            return self.max_R

    def get_quad(self):
        try:
            return self.quad_R
        except NameError:
            self.compute_quad()
            return self.quad_R

    def get_hashin(self):
        try:
            return self.hashin_R
        except NameError:
            self.compute_hashin()
            return self.hashin_R

    def compute_all(self):
        self.compute_max()
        self.compute_quad()
        self.compute_hashin()
        self.all_R = numpy.hstack((self.max_R, self.quad_R, self.hashin_R))
        return self.all_R

    def make_table(self):
        columns = 'FT_max,FC_max,MT_max,MC_max,S,(+),(-),FT_hash,FC_hash,MT_hash,MC_hash'.split(',')
        try:
          data = self.all_R
        except AttributeError:
          self.compute_all()
          data = self.all_R
        df = pd.DataFrame(data=data, columns=columns, dtype=float)
        plies = [('%i (%i) - %s' % (layer.index + 1, layer.theta, pos)
                 ) for layer in self.sim.laminate.layers for pos in ['B', 'T']
        ]
        df['Ply'] = plies
        rearanged_cols = ['Ply']+columns
        self.R_data = df[rearanged_cols]

        self.table_made = True
        return self.R_data

    def find_min(self):
        if not self.table_made:
            dump = self.make_table()
            self.table_made = True
        minimum_R = {}
        max_cols  = ('FT_max','FC_max','MT_max','MC_max','S')
        quad_cols = ('(+)',)
        hash_cols = ('FT_hash','FC_hash','MT_hash','MC_hash')
        all_cols = (max_cols,quad_cols,hash_cols)
        names = ('max','quad','hash')
        for cols,name in zip(all_cols,names):
            R_min_data = find_R_mins(self.R_data,cols)
            the_min_idx = R_min_data['Lowest R'].idxmin()
            minimum_R[name] = tuple(
                R_min_data[['Lowest R','Mode']].iloc[the_min_idx])

        return minimum_R

if __name__ == '__main__':
    import numpy

    #Create a Simulation object composed of a laminate.
    sim = Sim(laminate = Laminate('0_2/p25/0_2s',
                               materialID = 5, #5
                               core_thick = 0.01))
    #Define and apply load
    P = Q_(1000,'N'); b = Q_(0.11,'m'); L = Q_(0.51,'m')
    M1 = -P*L/(4*b); N1 = 0.5*P/b;
    M = Q_([M1.magnitude,0,0],M1.units)
    N = Q_([N1.magnitude,0,0],N1.units)
    #Apply load
    sim.apply_M(M)
    sim.apply_N(N)
    sim.solve()
    fail = FailureAnalysis(sim)
    R_data = fail.make_table()
# print M
# print N