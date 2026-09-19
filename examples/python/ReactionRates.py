# Only for dev purposes
import os

import matplotlib.pyplot as plt
from pySMOKEPostProcessor.postprocessor import PostProcessor

kineticFolder = os.path.join("..", "data", "ROPA", "kinetics")
resultsFolder = os.path.join("..", "data", "ROPA", "Output")

pp = PostProcessor(kineticFolder, resultsFolder)
rate = pp.GetReactionRates(reaction_name=['O2+H=O+OH'])[0]

_, temperature = pp.getTemperatureProfile()

fig_1, ax_1 = plt.subplots(nrows=1)
ax_1.plot(temperature, rate, 'r')
ax_1.set_xlabel('Temperature [K]', fontsize=15)
ax_1.set_ylabel(r'Reaction Rate $\left[ \dfrac{kmol}{m^{3}s} \right]$', fontsize=15)
ax_1.grid()

fig_2, ax_2 = plt.subplots(nrows=1)

names = ['O2+H=O+OH', 'H2+O=H+OH', 'O2+H(+M)=HO2(+M)']

rates = pp.GetReactionRates(reaction_name=names)
ax_2.plot(temperature, rates[0], 'r', label=names[0])
ax_2.plot(temperature, rates[1], 'b', label=names[1])
ax_2.plot(temperature, rates[2], 'g', label=names[2])
ax_2.set_xlabel('Temperature [K]', fontsize=15)
ax_2.set_ylabel(r'Reaction Rate $\left[ \dfrac{kmol}{m^{3}s} \right]$', fontsize=15)
ax_2.grid()
ax_2.legend()

plt.show()
