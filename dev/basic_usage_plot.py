import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

def f(x):
    return np.sin(x) - 0.1*x**2 + 1

rts = [[3.14959, 3.1496], [-4.42654, -4.42653], [-1.08205, -1.08204], [-3.10682, -3.10681]]

fig, axmain = plt.subplots()

xx = np.linspace(-6, 6, 10000)

xzoom = [-1.08206, -1.08203]
yzoom = [-0.00001, 0.00001]

axins = inset_axes(axmain, width="35%", height="45%", loc=4)
axins.set_xticks([])
axins.set_yticks([])

axins.set_xlim(xzoom)
axins.set_ylim(yzoom)

mark_inset(axmain, axins, loc1=1, loc2=3, alpha=0.5)

for ax in (axmain, axins):
    ax.plot(xx, f(xx))
    ax.axhline(0, color="k")

    for x, y in rts:
        ax.axvspan(x, y, alpha=0.7, color="C2")

fig.savefig("docs/basic_usage.png")
