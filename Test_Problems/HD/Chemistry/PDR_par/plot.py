from __future__ import print_function
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from prizmo_user_utils import Utils_prizmo
import matplotlib

font = {'family': 'normal',
        'weight': 'normal',
        'size': 14}

# matplotlib.rc('font', **font)


# ***************************
# ***************************
# ***************************
class Report:

    # ***************************
    def __init__(self):
        self.output_folder = "plots/"
        self.fname_status = "slab_final.out"
        self.fname_cooling = "cool.out"
        self.cadence = 30  # number of time-steps to skip for plotting

        self.format = "png"

        lu = Utils_prizmo()
        self.species_names = lu.prizmo_names
        self.reaction_names = lu.prizmo_reaction_names
        self.cmap = matplotlib.cm.get_cmap('viridis')

        if os.path.exists(self.output_folder):
            for fname in glob.glob(self.output_folder + "/*"):
                os.remove(fname)
        else:
            os.makedirs(self.output_folder)

    # ***************************
    def plot_all(self):
        self.plot_cooling()
        self.plot_variables()

    # *************************
    @staticmethod
    def compact_name(arg, latex=True):
        if arg.count("+") < 2:
            return arg
        else:
            if latex:
                return arg.replace("+", "") + ("$^{+%d}$" % arg.count("+"))
            else:
                return arg.replace("+", "") + ("+%d" % arg.count("+"))

    # **************************
    def plot_variables(self):

        labs = ["t", "z", "Tgas", "Tdust", "Ntot", "N(H2)", "N(CO)", "Av", "crflux", "k_H2", "k_CO", "k_C"] \
               + ["abundance %s" % x for x in self.species_names]
               # + ["flux %s" % x for x in self.reaction_names]

        print("loading slab results from " + self.fname_status)
        data = np.loadtxt(self.fname_status).T

        tidx = 0
        xidx = 7

        utime = sorted(np.unique(data[tidx]))

        for i, lab in enumerate(tqdm(labs)):

            plt.clf()
            for j, t in enumerate(utime):
                if j % self.cadence != 0 and t != utime[-1]:
                    continue
                idxs = data[tidx, :] == t
                jnorm = j / (len(utime) - 1)
                label = "%.2e yr" % t
                plt.loglog(data[xidx, idxs], data[i, idxs], label=label, color=self.cmap(jnorm))

            plt.legend(loc='best', fontsize=10)
            plt.title(lab)
            plt.xlabel("$A_v$")
            plt.ylabel(labs[i])
            plt.tight_layout()
            plt.savefig(self.output_folder +
                        "/plot_%s.%s" % (lab.replace(")", "").replace("(", "__"), self.format))


    # **************************
    def plot_cooling(self):
        lab_vars = ["t", "z", "Av"]
        lab_cool = ["H2", "CO", "atomic", "c_chem", "dust", "H2O", "BS"]
        lab_heat = ["CR", "phe", "ph", "PAH", "h_chem"]
        lab_skip = lab_vars + ["BS"]
        labs = lab_vars + lab_cool  + lab_heat
        if not os.path.isfile(self.fname_cooling):
            print("WARNING: file %s not found!" % self.fname_cooling)
            return
        data = np.loadtxt(self.fname_cooling).T

        tidx = 0
        xidx = 2
        utime = sorted(np.unique(data[tidx]))

        for j, t in enumerate(tqdm(utime)):
            if j % self.cadence != 0 and t != utime[-1]:
                continue
            plt.clf()
            idxs = data[tidx, :] == t
            tot_cool = np.zeros_like(data[xidx, idxs])
            tot_heat = np.zeros_like(tot_cool)

            for i, lab in enumerate(labs):
                if lab in lab_skip:
                    continue
                if lab in lab_cool:
                    tot_cool += data[i, idxs]
                    ls = "-"
                else:
                    tot_heat += data[i, idxs]
                    ls = "--"
                plt.loglog(data[xidx, idxs], data[i, idxs], label=lab, ls=ls)

            plt.loglog(data[xidx, idxs], tot_cool, label="tot_cool", lw=3, color="k", alpha=0.3)
            plt.loglog(data[xidx, idxs], tot_heat, label="tot_heat", lw=3, ls="--", color="k", alpha=0.3)

            jnorm = j / (len(utime) - 1)

            plt.legend(loc='lower left', fontsize=8, ncol=2, labelspacing=.2)
            plt.xlabel("$A_v$")
            plt.ylabel("cooling, heating [erg/cm$^{-3}$/s]")
            plt.ylim(1e-25, 3e-20)
            plt.title("t=%.2e yr" % t) #, color=self.cmap(jnorm))
            ax = plt.gca()
            ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10.0, numticks=10))
            locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1)
            ax.yaxis.set_minor_locator(locmin)
            ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            plt.tight_layout()
            plt.savefig(self.output_folder +
                        "/cool_%05d.%s" % (j, self.format))


if __name__ == "__main__":
    r = Report()
    r.plot_all()
