from composites.failure import FailureAnalysis
from composites.sim import Sim, ureg, Q_, make_test_sim
from composites.laminate import Laminate
import pandas as pd

def test_init():
	sim = make_test_sim()
	fail = FailureAnalysis(sim)

def test_ass5():
  sim = Sim(laminate = Laminate('0_2/p25/0_2s',
                               materialID = 5, #5
                               core_thick = 0.01))
  P = Q_(1000,'N'); b = Q_(0.11,'m'); L = Q_(0.51,'m')
  M1 = -P*L/(4*b); N1 = 0.5*P/b;
  M = Q_([M1.magnitude,0,0],M1.units)
  N = Q_([N1.magnitude,0,0],N1.units)
  sim.apply_M(M)
  sim.apply_N(N)
  sim.solve()
  df = sim.return_results()
  fail = FailureAnalysis(sim)
  R = fail.compute_all()
  R_data = fail.make_table()
  def print_R_analyz(data):
    the_min_idx = data['Lowest R'].idxmin()
    lowest_R = data['Lowest R'].iloc[the_min_idx]
    ply = data['Ply'].iloc[the_min_idx].split(' ')[0]
    print ("Lowest R is %.1f and occurs at top of ply number %s."
               % (lowest_R,ply))
    print "The load vectors R(M) and R(N) which cause failure are:"
    print "R(M) [N] : "
    print [round(M1.magnitude*lowest_R,2),0,0]
    print "R(N) [N/m] : "
    print [round(N1.magnitude*lowest_R,2),0,0]

  def find_min(dataframe,columns):
    lowest_Rs = []; plys = []; 
    modes = [x.split('_')[0] for x in columns]
    nonzero_dat = dataframe[dataframe != 0]
    for col in columns:
        the_data = nonzero_dat[col]
        min_index = the_data.idxmin()
        lowest_Rs.append(the_data.iloc[min_index])
        plys.append(dataframe['Ply'].iloc[min_index])
    return pd.DataFrame({'Mode':modes,
            'Lowest R':lowest_Rs,
            'Ply':plys
            })[['Mode','Lowest R','Ply']]

  max_data = find_min(R_data,'FT_max,FC_max,MT_max,MC_max,S'.split(','))
  print_R_analyz(max_data)
  quad_data = find_min(R_data,['(+)'])
  print_R_analyz(quad_data)
  print R_data
  hashin_data = find_min(R_data,'FT_hash,FC_hash,MT_hash,MC_hash'.split(','))
  hashin_data
  print_R_analyz(hashin_data)

  fail.find_min




	