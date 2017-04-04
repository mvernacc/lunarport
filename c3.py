''' Plot payload vs C3 for Eexploration Upper Stage (EUS) refueled at Lunarport.'''
import numpy as np
from matplotlib import pyplot as plt
from dv_r import dv_circ_to_c3, c3_to_dv_circ



# Acceleration due to gravity at Earth's surface [units: meter second**-2].
g = 9.81

def sls_delta_v(m_pl):
    ''' Compute the delta-v which SLS can impart to a given payload.'''
    # specific impulses [second].
    Isp_0 = 259
    Isp_cv = 452
    Isp_2 = 460

    # Masses in megagrams
    m_bbo = 204
    m_bp = 1256.8
    m_cbo = 113.9
    m_cp = 960
    m_cbbo = 881
    m_2bo = 13.6
    m_2p = 129

    dv_0 = Isp_0 * g * np.log((m_bbo + m_cbo + m_cp + m_bp + m_2bo + m_2p + m_pl) /
        (m_bbo + m_cbbo + m_2bo + m_2p + m_pl))
    dv_1 = Isp_cv * g * np.log((m_cbbo + m_2bo + m_2p + m_pl) /
        (m_cbo + m_2bo + m_2p + m_pl))
    dv_2 = Isp_cv * g * np.log((m_2bo + m_2p + m_pl) /
        (m_2bo + m_pl))

    return (dv_0 + dv_1 + dv_2)


# EUS burnout mass [units: megagram]
m_eus_bo = 13.6

# EUS propellant mass [units: megagram]
m_eus_prop = 129.

# Range of c3 values [m**2 s**-2]
c3 = np.linspace(0, 170e6, 100)

# Convert c3 to delta v applied at L1 [units: meter ssecond**-1]
# Rough guess at conversion, dv from an orbit matching the moon's + some margin.
r_moon = 364e6
dv_L1 = c3_to_dv_circ(r_moon, c3) + 500
# EUS engine ISP (assuming RL10) [units: second]
Isp = 460. 
# Use Tsiolkovsky equation to find final mass [units: megagram]
m_f = (m_eus_prop + m_eus_bo) * np.exp(-dv_L1 / (Isp * g))
# Payload mass [units: megagram]
m_pl = m_f - m_eus_bo
# Knock out all payload masses > than the payload capacity of SLS to L1.
m_pl[m_pl > 39] = 39
# Plot results
plt.plot(c3*1e-6, m_pl, color='blue', linestyle='--', label='EUS refueled at L1, direct escape')


# Compare to EUS refueled at L1, escape via Earth and moon flyby. Use data from SpaceWorks study.
# http://www.ulalaunch.com/uploads/docs/Published_Papers/Exploration/Cryogenic_Propulsive_Stage_Mission_Sensitivity_Studies_-_Earth-Moon_L1_Departure.pdf
dv_L1_flyby = 520 + 1770. / 27e6 * c3
# Use Tsiolkovsky equation to find final mass [units: megagram]
m_f_flyby = (m_eus_prop + m_eus_bo) * np.exp(-dv_L1_flyby / (Isp * g))
# Payload mass [units: megagram]
m_pl_flyby = m_f_flyby - m_eus_bo
# Knock out all payload masses > than the payload capacity of SLS to L1.
m_pl_flyby_tug = np.copy(m_pl_flyby)
m_pl_flyby[m_pl_flyby > 39] = 39
# Plot results
plt.plot(c3*1e-6, m_pl_flyby, color='blue', linestyle='-', label='EUS refueled at L1, Earth + Moon flyby')
plt.plot(c3*1e-6, m_pl_flyby_tug, color='purple', linestyle='--',
    label='EUS tug to L1, refuel at L1, Earth + Moon flyby')

# Now compare to a full EUS in LEO.
# Earth gravitation parameter [meter**3 second**-2].
mu = 3.986e14
# LEO radius [meter].
r_leo = 6371e3 + 300e3
# Delta-v [meter second**-1].
dv_leo = c3_to_dv_circ(r_leo, c3)
# Use Tsiolkovsky equation to find final mass [units: megagram]
m_f_leo = (m_eus_prop + m_eus_bo) * np.exp(-dv_leo / (Isp * g))

# Payload mass [units: megagram]
m_pl_leo = m_f_leo - m_eus_bo

# Plot results
plt.plot(c3*1e-6, m_pl_leo, color='red', label='EUS refueled at 300km LEO')

# Plot my calculated SLS results for comparison
m_pl_sls_me = np.linspace(0, 50, 100)
dv_sls_me = sls_delta_v(m_pl_sls_me)
dv_leo = 9500 + 500
c3_sls_me = dv_circ_to_c3(r_leo, dv_sls_me - dv_leo)
plt.plot(c3_sls_me[c3_sls_me>100e6]*1e-6, m_pl_sls_me[c3_sls_me>100e6],
    color='green', linestyle=':')

# Load comparison data
d = np.genfromtxt('SLS_LUS.csv', delimiter=',')
c3_sls = d[:,0]
m_pl_sls = d[:,1]
d = np.genfromtxt('Delta_IV_Heavy.csv', delimiter=',')
c3_d4h = d[:,0]
m_pl_d4h = d[:,1]

# Plot comparison data
plt.plot(c3_sls, m_pl_sls, color='green', label='EUS launched on SLS')
# plt.plot(c3_d4h, m_pl_d4h, color='grey', label='Delta IV Heavy')

# Add c3 references
bbox = bbox=dict(facecolor='white', edgecolor='white', alpha=0.5)
plt.axvline(x=11, color='black')
plt.text(12, 7, '< Mars', bbox=bbox)
plt.axvline(x=40, color='black', ymax=1)
plt.text(42, 7, '< Mars\nfree return')
plt.axvline(x=91, color='black', ymax=0.7)
plt.text(93, 9, '< Europa', bbox=bbox)
plt.axvline(x=106, color='black', ymax=0.7)
plt.text(108, 20, '< Enceladus', bbox=bbox)
plt.axvline(x=128, color='black', ymax=0.7)
plt.text(130, 13, '< Uranus', bbox=bbox)
plt.text(145, 20, 'Sun\nescape >', bbox=bbox)

# Dress up the plot.
fig = plt.gcf()
fig.set_size_inches(8, 8)
plt.legend(loc='best', framealpha=0.5, fontsize=12)
plt.title('Payload capacity of SLS EUS after refueling at Lunarport')
plt.xlabel('C3 injection energy [km^2 s^-2]')
plt.ylabel('Payload capacity [Mg]')
plt.xlim([0, 170])
plt.ylim([0, 100])
plt.grid(True)

# Plot percent improvement in payload.
# Covert basis of SLS data
m_pl_sls = np.interp(c3, np.flip(c3_sls_me, 0), np.flip(m_pl_sls_me, 0))
# Payload mass ratios
mr = m_pl / m_pl_sls
mr_flyby = m_pl_flyby / m_pl_sls
mr_flyby_tug = m_pl_flyby_tug / m_pl_sls
mr_leo = m_pl_leo / m_pl_sls

# import pdb; pdb.set_trace()

plt.figure(figsize=(8, 6))
plt.plot(c3*1e-6, mr, color='blue', linestyle='--', label='EUS refueled at L1, direct escape')
plt.plot(c3*1e-6, mr_flyby, color='blue', linestyle='-', label='EUS refueled at L1, Earth + Moon flyby')
plt.plot(c3*1e-6, mr_flyby_tug, color='purple', linestyle='--',
    label='EUS tug to L1, refuel at L1, Earth + Moon flyby')
plt.plot(c3*1e-6, mr_leo, color='red', label='EUS refueled at 300km LEO')

plt.axvline(x=11, color='black', ymax=0.75)
plt.text(12, 0.6, '< Mars', bbox=bbox)
plt.axvline(x=40, color='black', ymax=0.75)
plt.text(42, 0.6, '< Mars\nfree return')
plt.axvline(x=91, color='black', ymax=0.75)
plt.text(93, 1, '< Europa', bbox=bbox)
plt.axvline(x=106, color='black', ymax=0.75)
plt.text(108, 0.6, '< Enceladus', bbox=bbox)
plt.axvline(x=128, color='black', ymax=1)
plt.text(130, 1, '< Uranus', bbox=bbox)
plt.text(145, 0.6, 'Sun\nescape >', bbox=bbox)

plt.legend(loc='best', framealpha=0.5, fontsize=12)
plt.title('Improvement of SLS EUS payload capacity with Lunarport')
plt.xlabel('C3 injection energy [km^2 s^-2]')
plt.ylabel('Payload multiplication due to refueling [-]')
plt.xlim([0, 170])
plt.ylim([0, 5])
plt.grid(True)

plt.show()

