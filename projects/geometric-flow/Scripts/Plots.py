import matplotlib.pyplot as plt 
import numpy as np 
import matplotlib as mpl
def load_pretty_figure_setup():
    # if mpl seems to work with the wrong latex executable, uncomment and adapt path according to your installation
    # os.environ["PATH"] = os.path.expanduser("~/texlive/bin/x86_64-linux") + ":" + os.environ["PATH"]
    
    _preamble_shared = R"""
        \usepackage{graphicx}
        \DeclareMathOperator{\arcsinh}{arcsinh}
        \DeclareMathOperator{\km}{k_\mathrm{m}}
        \DeclareMathOperator{\fbi}{f_\mathrm{bi}}
        \DeclareMathOperator{\eps}{\epsilon_{\mathrm{mc}}}
        \DeclareMathOperator{\epscrit}{\epsilon_{\mathrm{mc}}^*}
        \DeclareMathOperator{\uf}{u_{\mathrm{f}}}
        \DeclareMathOperator{\kBT}{k_\mathrm{B}T}
        """[
        1:
    ]

    def mpl_rcParams_avenir():
        rcParams = {}
        rcParams["font.family"] = "sans-serif"
        rcParams["font.cursive"] = ["Optima"]
        rcParams["text.usetex"] = True
        # rcParams['text.latex.unicode']= True
        rcParams["pgf.texsystem"] = "lualatex"
        rcParams["pgf.rcfonts"] = False
        rcParams["pgf.preamble"] = (
            R"""
        \usepackage[utf8x]{inputenc}
        \usepackage[T1]{fontenc}
        \usepackage{fontspec}
        \usepackage{amsmath}
        \setmainfont{Avenir}[Scale=.9]
        \renewcommand{\setmainfont}{}
        \renewcommand{\sffamily}{}
        """[
                1:
            ]
            + "\n"
            + _preamble_shared
        )
        return rcParams
    
    def rc_params_setup():
        mpl.rcParams["font.family"] = "serif"
        mpl.rcParams["text.usetex"] = True
        mpl.rcParams["figure.constrained_layout.use"] = True
        mpl.rcParams.update(mpl_rcParams_avenir())
        # mpl.rcParams["pgf.texsystem"] = "lualatex"
        # mpl.rcParams["text.latex.preamble"] = mpl.rcParams['pgf.preamble'] #R"\usepackage{amsmath}\usepackage{lmodern}"
        mpl.rcParams["text.latex.preamble"] = (
            R"""
        \usepackage{lmodern}
        \usepackage{amsmath}
        """
            + "\n"
            + _preamble_shared
        )

    rc_params_setup()
    print("Pretty figure set-up loaded.")






load_pretty_figure_setup()



# I want to do A piechart plot

# labels = 'Remeshing', 'Saving mesh', 'Integrating'
# sizes = [28643,3231, 11327126]

# fig, ax = plt.subplots()
# ax.pie(sizes, labels=labels)
# plt.show()


def plot1():



    from matplotlib.ticker import LinearLocator

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    # Make data.
    rc = 1.0
    X = np.arange(0, 2, 0.25)
    Y = np.arange(-5, 10, 0.25)
    X, Y = np.meshgrid(X, Y)

    # So now i need to write the function
    Z = -1*(2*X**3/(rc**3)- 3*X**2/(rc**2) +1 )*Y*(X<rc)
    # R = np.sqrt(X**2 + Y**2)
    # Z = np.sin(R)

    # Plot the surface.
    surf = ax.plot_surface(X, Y, Z, cmap="coolwarm",
                        linewidth=0, antialiased=False)

    # Customize the z axis.
    ax.set_zlim(-10.01, 10.01)
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel('$z$')
    ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter('{x:.02f}')

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()



def plot2():
    # Here we reload
    dir ="../Results/Wrapping_July_rescale/"
    file = dir+"EnergyScale.txt"
    data = np.loadtxt(file,usecols=(1))
    bins_trial= np.linspace(-0.011,0.011,200)
    plt.yscale('log')
    plt.hist(data,bins=bins_trial)
    plt.show()
    print(max(data))
    print(min(data))


# OK 
def spacing_contour():
    xmin = 0.015625
    xmax = 1.5625

    xs = 10**np.linspace(np.log10(xmin), np.log10(xmax), 10)

    spacing = xs[1]/xs[0]

    xs2 = []
    xs2.append(xs[0]/(spacing*spacing))

    xs2.append(xs[0]/(spacing))
    for i in xs:
        xs2.append(i)

    xs2.append(xs[-1]*spacing)
    # xs2.append(xs[-1]*spacing*spacing)

    # THis is useful cause i can add things 


    xs2 = np.array(xs2)
    print(xs2)

    ys = np.linspace(0,5,11)
    lamda = 1/np.sqrt(xs2)
    X,Y = np.meshgrid(lamda,ys)
    Z = np.ones_like(X)
    plt.scatter(X,Y)
    plt.xscale('log')
    plt.axvline(10)
    plt.show()


def contour():
    folder = "../Results/WrappingPhaseSpaceNewFast/"
    filepath = folder + "Coverage_data.txt"

    Data = np.loadtxt(filepath,delimiter = ' ', skiprows = 1, usecols = (1,2,3,4,5,6))
    # 
    plt.scatter(Data[:,2],Data[:,5])
    # plt.show()
    plt.clf()
    
    lamda = np.sqrt(Data[:,1]/Data[:,0])
    a = Data[:,3]
    print(a)
    X = a/lamda 
    Y = Data[:,2]/Data[:,0]

    # Create a pcolormesh by binning the scattered values (NumPy-only)
    nx, ny = 24, 16
    # use log-spaced x edges because X is plotted on a log scale
    x_min, x_max = np.nanmin(X), np.nanmax(X)
    if x_min <= 0:
        x_min = np.nextafter(np.nanmin(X[X>0]), 0)
    x_edges = np.logspace(np.log10(x_min), np.log10(x_max), nx + 1)
    y_edges = np.linspace(np.nanmin(Y), np.nanmax(Y), ny + 1)

    sum_grid, _, _ = np.histogram2d(X, Y, bins=[x_edges, y_edges], weights=Data[:,5])
    cnt_grid, _, _ = np.histogram2d(X, Y, bins=[x_edges, y_edges])
    avg_grid = np.full_like(sum_grid, np.nan, dtype=float)
    mask = cnt_grid > 0
    avg_grid[mask] = sum_grid[mask] / cnt_grid[mask]

    Xe, Ye = np.meshgrid(x_edges, y_edges)
    plt.pcolormesh(Xe, Ye, avg_grid.T, shading='auto', cmap='viridis')
    plt.xscale('log')
    plt.ylabel(r"$w/\sigma$")
    plt.xlabel(r"$a\,\sqrt{ \sigma / \kappa_\text{B}}$")
    # plt.axvline(x=4.4, ls='dashed', color='black')
    # plt.axhline(y=1.37, ls='dashed', color='black')
    cbar = plt.colorbar(ticks=[0.1*i for i in range(11)])
    cbar.set_label('Coverage', rotation=90)
    plt.savefig(folder + "PhaseSpaceVesicle.png",bbox_inches = 'tight')
    plt.show()
    

    print("THe values of X are {}".format(X))
    print("The values of Y are {}".format(Y))
    
    return

# contour()

def coverage_lowest_energy(batch="../Results/Phase_one_vesicle", legs=("in", "out")):
    """All the points of the Coverage_data_*.csv of a folder (one file per chain, Coverage_data_chain.py), one row per point: where a chain has two numbers
    for the same scanned value, the way in and the way back, the one with the lower Total_E is kept (the same constants, so the lower energy is the more
    stable state; a tie keeps the way in). The cold start is left out. Returns the rows (dicts of strings, as the csv has them, plus "file").
    batch: the folder, or a glob of csv's. Needs the Total_E column: the csv's written before 2026-10-17 do not have it, run Collect_chain_results.sh again."""
    import csv
    import glob
    import os

    files = sorted(glob.glob(os.path.join(batch, "Coverage_data_*.csv") if os.path.isdir(batch) else batch))
    if not files:
        raise FileNotFoundError("no Coverage_data_*.csv in " + str(batch))
    best, seen = {}, {}
    for fn in files:
        with open(fn) as f:
            for r in csv.DictReader(f):
                if r["leg"] not in legs:
                    continue
                if "Total_E" not in r:
                    raise KeyError("no Total_E column in {}: run Collect_chain_results.sh again".format(fn))
                r["file"] = os.path.basename(fn)
                key = (fn, r["scan_var"], float(r["scan_value"]))
                seen[key] = seen.get(key, 0) + 1
                if key not in best or float(r["Total_E"]) < float(best[key]["Total_E"]):
                    best[key] = r
    pairs = [k for k, n in seen.items() if n > 1]
    back = sum(best[k]["leg"] == "out" for k in pairs)
    print("{} points from {} files; {} have both legs, the way back has the lower energy in {}".format(len(best), len(files), len(pairs), back))
    return list(best.values())


def phase_one_vesicle_map(batch="../Results/Phase_one_vesicle", legs="lowest", column="CoveredArea", kappa_over_KB=1.0, w_over_KI=1.0,
                          nx=24, ny=16, save="Phase_one_vesicle_coverage.png", show=True):
    """pcolormesh of the chains of Phase_one_vesicle in the fashion of contour(): the points are binned (nx bins log-spaced in x, ny linear in y) and the
    colour of a bin is the mean of the column of its points; bins without points are left empty.
        x = a sqrt(KA / (kappa_over_KB KB))  (log axis; kappa_over_KB = 1 is a sqrt(sigma/KB) as in contour(), 0.5 the a/lambda of deserno_fig5),   y = w_over_KI KI / KA = w/sigma
    batch: the folder with the Coverage_data_<Tag>.csv of every chain (Coverage_data_chain.py, run by Collect_chain_results.sh on the cluster), or a glob of csv's.
    legs: "lowest" (default) takes every point once, the way in or the way back, whichever has the lower Total_E (coverage_lowest_energy); ("in",) the way
    along the KB list, ("out",) the way back, ("in", "out") both (a bin then takes the mean of the two numbers of its points).
    column: CoveredArea (as contour()), CoverageUnion, cap_cfg or f_ad.
    save: file name of the figure, a bare name goes to the folder of the csv's (None: not saved).
    """
    import csv
    import glob
    import os

    files = sorted(glob.glob(os.path.join(batch, "Coverage_data_*.csv") if os.path.isdir(batch) else batch))
    if not files:
        raise FileNotFoundError("no Coverage_data_*.csv in " + str(batch))
    if legs == "lowest":
        rows = [r for r in coverage_lowest_energy(batch) if float(r["KA"]) > 0]                  # w/sigma = KI/KA is undefined at zero tension
    else:
        rows = []
        for fn in files:
            with open(fn) as f:
                rows += [r for r in csv.DictReader(f) if r["leg"] in legs and float(r["KA"]) > 0]
    KA, KB, KI, r0 = (np.array([float(r[k]) for r in rows]) for k in ("KA", "KB", "KI", "BeadRadius"))
    C = np.array([float(r[column]) for r in rows])
    X = r0 * np.sqrt(KA / (kappa_over_KB * KB))
    Y = w_over_KI * KI / KA
    good = np.isfinite(X) & np.isfinite(Y) & np.isfinite(C) & (X > 0)
    X, Y, C = X[good], Y[good], C[good]
    print("{} points from {} files".format(len(X), len(files)))

    x_min, x_max = X.min(), X.max()
    y_min, y_max = (Y.min(), Y.max()) if Y.max() > Y.min() else (Y.min() - 0.5, Y.max() + 0.5)
    x_edges = np.logspace(np.log10(x_min), np.log10(x_max), nx + 1)
    y_edges = np.linspace(y_min, y_max, ny + 1)
    sum_grid, _, _ = np.histogram2d(X, Y, bins=[x_edges, y_edges], weights=C)
    cnt_grid, _, _ = np.histogram2d(X, Y, bins=[x_edges, y_edges])
    avg_grid = np.full_like(sum_grid, np.nan, dtype=float)
    mask = cnt_grid > 0
    avg_grid[mask] = sum_grid[mask] / cnt_grid[mask]

    fig, ax = plt.subplots()
    ax.pcolormesh(x_edges, y_edges, np.ma.masked_invalid(avg_grid.T), shading="flat", cmap="viridis", vmin=0.0 if column != "f_ad" else None, vmax=1.0 if column != "f_ad" else None)
    ax.set_xscale("log")
    ax.set_ylabel(r"$w/\sigma$")
    ax.set_xlabel(r"$a\,\sqrt{ \sigma / \kappa_\text{B}}$" if kappa_over_KB == 1.0 else r"$a/\lambda$")
    cbar = fig.colorbar(ax.collections[0], ticks=[0.1 * i for i in range(11)] if column != "f_ad" else None)
    cbar.set_label("Coverage", rotation=90)
    if save is not None:
        path = save if os.path.dirname(save) else os.path.join(batch if os.path.isdir(batch) else os.path.dirname(files[0]), save)
        fig.savefig(path, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax


def EdgePlot():
    dir = "../Results/Mem_shape_PR/41/"
    filepath = dir + "Edge_data_step_9850.txt"

    E_data = np.loadtxt(filepath)   

    # plt.scatter(E_data[:,1],np.abs(E_data[:,3]))
    # plt.show()
    plt.xlabel("Edge length")
    plt.ylabel("Count")
    plt.hist(E_data[:,1],bins='auto')
    plt.show()

    plt.hist(np.abs(E_data[:,2]),bins='auto')
    plt.xlabel("Dihedral angle")
    plt.ylabel("Count")
    plt.show()

    plt.hist(E_data[:,3],bins='auto')
    plt.xlabel(r"$\bar{H}$")
    plt.ylabel("Count")
    plt.show()

    plt.scatter(E_data[:,1],E_data[:,3])
    plt.xlabel("Edge length")
    plt.ylabel(r"$$\bar{H}$$")
    plt.show()

# EdgePlot()

def bending_phase(filepath="../Results/TwoBeadsCov/Bending_data_phase.txt", column="bending", save=None,
                  contour_levels=10):
    """pcolormesh of the final energy over the (bead distance, target coverage) phase space.

    File columns: run, distance, target coverage, bending energy, coverage energy, coverage strength.
    column: "bending" or "coverage", selects the energy that is coloured.
    contour_levels: number of (or list of values for) the contour lines drawn on top, None to skip them.
    """
    columns = {
        "bending": (3, r"$E_\mathrm{bend}$"),
        "coverage": (4, r"$E_\mathrm{cov}$"),
        "total": (6, r"$E_\mathrm{total}$"),
    }
    col, label = columns[column]

    # Lines starting with '#' are skipped; the runs come unordered
    data = np.loadtxt(filepath, comments="#")
    distance, coverage, energy = data[:, 1], data[:, 2], data[:, col]

    # Sorted unique grid values, the pivot fills Z[coverage, distance]
    xs = np.unique(distance)
    ys = np.unique(coverage)
    Z = np.full((len(ys), len(xs)), np.nan)
    ix = np.searchsorted(xs, distance)
    iy = np.searchsorted(ys, coverage)
    Z[iy, ix] = energy

    # Cell edges halfway between the centres, so each cell is centred on its run
    def edges(c):
        mid = 0.5 * (c[1:] + c[:-1])
        return np.concatenate(([c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]))

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(edges(xs), edges(ys), Z, cmap="viridis", shading="flat")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(label, rotation=90)

    # Contour lines through the cell centres, on top of the pcolormesh (NaN cells are left out)
    if contour_levels is not None:
        cs = ax.contour(xs, ys, np.ma.masked_invalid(Z), levels=contour_levels,
                        colors="black", linewidths=0.8)
        ax.clabel(cs, fontsize="small", fmt="%.2f")

    ax.set_xlabel(r"Bead distance $d$")
    ax.set_ylabel(r"Target coverage")
    ax.set_xticks(xs)
    ax.set_yticks(ys)
    ax.set_aspect("auto")

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    plt.show()
    return fig, ax

# bending_phase(filepath="../Results/TwoBeadsCov/Bending_data_phase_Bendi.txt", column="bending", save="../Results/TwoBeadsCov/Bending_phase_Bendi_bendingE.pdf")


def obtained_coverage(filepath="../Results/TwoBeadsCov/Bending_data_phase.txt", sign=-1,
                      deviation=False, save=None):
    """Obtained coverage vs target coverage, one line per bead distance.

    The coverage energy is E = sum_beads K (Omega_i - target)^2 (E_Handler::E_Coverage, no 1/2).
    Both beads have the same target and, by symmetry, the same Omega, so E = 2 K (Omega - target)^2
    and |Omega - target| = sqrt(E / (2 K)). The energy loses the sign of the deviation, so it is
    set by `sign` (-1: under-covered, +1: over-covered).
    deviation: plot Omega - target instead of Omega, to resolve the (small) offset from y = x.
    """
    # File columns: run, distance, target coverage, bending energy, coverage energy, K
    data = np.loadtxt(filepath, comments="#")
    distance, target, E_cov, K = data[:, 1], data[:, 2], data[:, 4], data[:, 5]

    n_beads = 2
    delta = sign * np.sqrt(E_cov / (n_beads * K))
    obtained = target + delta
    y = delta if deviation else obtained

    distances = np.unique(distance)
    colors = plt.get_cmap("viridis")(np.linspace(0, 0.9, len(distances)))

    fig, ax = plt.subplots()
    for d, c in zip(distances, colors):
        m = distance == d
        order = np.argsort(target[m])
        ax.plot(target[m][order], y[m][order], "o-", color=c, ms=4, lw=1.2,
                label=r"$d=%.1f$" % d)

    if deviation:
        ax.axhline(0, color="black", ls="dashed", lw=0.8)
        ax.set_ylabel(r"Obtained $-$ target coverage")
    else:
        ax.plot([0, 1], [0, 1], color="black", ls="dashed", lw=0.8)
        ax.set_ylabel(r"Obtained coverage")
        ax.set_ylim(0, 1)
    ax.set_xlabel(r"Target coverage")
    ax.set_xlim(0, 1)
    ax.legend(title="Bead distance", fontsize="small", ncol=2)

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    plt.show()
    return fig, ax

# obtained_coverage()


def bending_vs_distance(filepath="../Results/TwoBeadsFast/Final_energies.txt", energy="Bending_tan", save=None):
    """Final energy vs distance between the beads, from the Final_energies file of Two_bead.py.

    File columns: directory, distance, <energy terms...>, Total_E (the header line names them).
    energy: name of the column to plot, by default the bending energy.
    """
    with open(filepath) as f:
        names = f.readline().lstrip("#").split()
    col = names.index(energy)

    # The first column is the run directory (a string), so only read distance and the energy
    data = np.loadtxt(filepath, comments="#", usecols=(1, col), ndmin=2)
    order = np.argsort(data[:, 0])
    distance, E = data[order, 0], data[order, 1]

    fig, ax = plt.subplots()
    ax.plot(distance, E, "o-", ms=4, lw=1.2)
    ax.set_xlabel(r"Bead distance $d$")
    ax.set_ylabel(r"$E_\mathrm{bend}$" if energy == "Bending_tan" else energy.replace("_", " "))

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    plt.show()
    return fig, ax

# bending_vs_distance(energy="Surface_tension")


def coverage_phase(filepath="../Results/WrappingPhaseVaryingKA/Coverage_data.txt", column="union", save=None):
    """pcolormesh of the coverage over the (KA, interaction strength KI) phase space.

    File columns: DIR KA KB KI BeadRadius rc CoveredArea CoverageUnion MultilayerFrac.
    column: "union" (CoverageUnion) or "area" (CoveredArea), selects the coverage that is coloured.
    """
    columns = {
        "union": (7, "Coverage"),
        "area": (6, "Covered area"),
    }
    col, label = columns[column]

    # Header starts with '#'; usecols skips the DIR string. The runs come unordered
    data = np.loadtxt(filepath, comments="#", usecols=(1, 3, col))
    KA, KI, cov = data[:, 0], data[:, 1], data[:, 2]

    # Sorted unique grid values, the pivot fills Z[strength, KA]
    xs = np.unique(KA)
    ys = np.unique(KI)
    Z = np.full((len(ys), len(xs)), np.nan)
    Z[np.searchsorted(ys, KI), np.searchsorted(xs, KA)] = cov

    # Cell edges halfway between the centres, so each cell is centred on its run
    def edges(c, log=False):
        if log:
            return np.exp(edges(np.log(c)))
        mid = 0.5 * (c[1:] + c[:-1])
        return np.concatenate(([c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]))

    # KA is sampled roughly geometrically, so it goes on a log axis
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(edges(xs, log=True), edges(ys), Z, cmap="viridis", shading="flat",
                         vmin=0, vmax=1)
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(label, rotation=90)

    ax.set_xscale("log")
    ax.set_xlabel(r"$K_A$")
    ax.set_ylabel(r"Interaction strength")
    ax.set_yticks(ys)

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    plt.show()
    return fig, ax

# coverage_phase()


def planar_phase(filepath="../Results/Wrapping_planar/Coverage_data.txt", column="CoverageUnion",
                 contour_levels=None, xlim=(0.75, 3), ylim=(0.0, 2.0), save=None):
    """pcolormesh of the planar membrane phase space, from the Coverage_data.txt of PostProcessing.

    x = KI r^2 / KB (adhesion strength), y = KA r^2 / KB (surface tension), with r the bead radius.
    File columns: DIR KA KB KI BeadRadius rc CoveredArea CoverageUnion MultilayerFrac BeadX Area <energies> Total_E,
    the header line names them.
    column: name of the column that is coloured, e.g. "CoverageUnion", "CoveredArea", "Bending_tan", "Total_E".
    contour_levels: number of (or list of values for) the contour lines drawn on top, None for no lines.
    xlim, ylim: axis limits.
    """
    with open(filepath) as f:
        names = f.readline().lstrip("#").split()
    col = {n: i + 1 for i, n in enumerate(names[1:])}  # the DIR column is a string: usecols skips it

    # The first header line can be repeated when PostProcessing appends to an old file
    data = np.loadtxt(filepath, comments="#", usecols=range(1, len(names)), ndmin=2)
    KA, KB, KI, r = (data[:, col[n] - 1] for n in ("KA", "KB", "KI", "BeadRadius"))
    value = data[:, col[column] - 1]
    # KB = KB/2
    x = np.round( KI * r**2 / (KB), 6)
    y = np.round(KA * r**2 / (KB), 6)
    
    # Sorted unique grid values, the pivot fills Z[y, x]
    xs = np.unique(x)
    ys = np.unique(y)
    Z = np.full((len(ys), len(xs)), np.nan)
    Z[np.searchsorted(ys, y), np.searchsorted(xs, x)] = value

    # Cell edges halfway between the centres, so each cell is centred on its run
    def edges(c):
        if len(c) == 1:
            return np.array([c[0] - 0.5, c[0] + 0.5])
        mid = 0.5 * (c[1:] + c[:-1])
        return np.concatenate(([c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]))

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(edges(xs), edges(ys), Z, cmap="viridis", shading="flat")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(column.replace("_", " "), rotation=90)

    # Contour lines through the cell centres, on top of the pcolormesh (NaN cells are left out)
    if contour_levels is not None:
        cs = ax.contour(xs, ys, np.ma.masked_invalid(Z), levels=contour_levels,
                        colors="white", linewidths=0.8)
        ax.clabel(cs, fontsize="small", fmt="%.2f")
    

    ax.set_xlabel(r"$K_I r^2 / K_B$")
    ax.set_ylabel(r"$K_A r^2 / K_B$")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("auto")
    ax.axvline(x=1.0, ls="dashed", color="black", lw=0.8)

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    plt.show()
    return fig, ax

def planar_phase_comp(filepath_flat="../Results/Wrapping_planar_excess/Coverage_data.txt",
                      filepath_full="../Results/Wrapping_planar_excess/Coverage_data_full.txt",
                      column="CoverageUnion", contour_levels=None, xlim=(0.75, 3), ylim=(0.0, 2.0), save=None):
    """planar_phase of the lower energy state of each point, flat start vs full (wrapped) start.

    The two files are the Coverage_data.txt of the batch that starts from the flat disk and of the one
    that starts from Planar_full.obj (same axes as planar_phase). At every (x, y) the run with the lower
    Total_E is kept and its `column` is coloured; a point missing in one file takes the other one.
    """
    def load(filepath):
        with open(filepath) as f:
            names = f.readline().lstrip("#").split()
        col = {n: i for i, n in enumerate(names[1:])}  # the DIR column is a string: usecols skips it
        data = np.loadtxt(filepath, comments="#", usecols=range(1, len(names)), ndmin=2)
        KA, KB, KI, r = (data[:, col[n]] for n in ("KA", "KB", "KI", "BeadRadius"))
        KB = KB/2  # same axes as planar_phase
        x = np.round(KI * r**2 / (2*KB), 6)
        y = np.round(KA * r**2 / KB, 6)
        return x, y, data[:, col[column]], data[:, col["Total_E"]]

    x_a, y_a, v_a, E_a = load(filepath_flat)
    x_b, y_b, v_b, E_b = load(filepath_full)

    # Sorted unique grid values over both files, the pivot fills [y, x]
    xs = np.unique(np.concatenate((x_a, x_b)))
    ys = np.unique(np.concatenate((y_a, y_b)))

    def grid(x, y, values):
        G = np.full((len(ys), len(xs)), np.nan)
        G[np.searchsorted(ys, y), np.searchsorted(xs, x)] = values
        return G

    value_flat, value_full = grid(x_a, y_a, v_a), grid(x_b, y_b, v_b)
    # NaN energy means no run: that start can not win
    E_flat = np.where(np.isnan(grid(x_a, y_a, E_a)), np.inf, grid(x_a, y_a, E_a))
    E_full = np.where(np.isnan(grid(x_b, y_b, E_b)), np.inf, grid(x_b, y_b, E_b))
    full_wins = E_full < E_flat
    Z = np.where(full_wins, value_full, value_flat)

    both = np.isfinite(E_flat) & np.isfinite(E_full)
    print("The full start has the lower energy in {} of {} points".format((full_wins & both).sum(), both.sum()))

    # Cell edges halfway between the centres, so each cell is centred on its run
    def edges(c):
        if len(c) == 1:
            return np.array([c[0] - 0.5, c[0] + 0.5])
        mid = 0.5 * (c[1:] + c[:-1])
        return np.concatenate(([c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]))

    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(edges(xs), edges(ys), Z, cmap="viridis", shading="flat")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(column.replace("_", " "), rotation=90)

    # Contour lines through the cell centres, on top of the pcolormesh (NaN cells are left out)
    if contour_levels is not None:
        cs = ax.contour(xs, ys, np.ma.masked_invalid(Z), levels=contour_levels,
                        colors="white", linewidths=0.8)
        ax.clabel(cs, fontsize="small", fmt="%.2f")

    ax.set_xlabel(r"$K_I r^2 / K_B$")
    ax.set_ylabel(r"$K_A r^2 / K_B$")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("auto")
    ax.axvline(x=1.0, ls="dashed", color="black", lw=0.8)

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    plt.show()
    return fig, ax

# planar_phase_comp(contour_levels=[0.25, 0.5, 0.75])

# Axes of the structural wrapping phase diagram of Deserno (arXiv:cond-mat/0303656, Fig. 2), in terms of the
# constants of the simulation (verified against the C++ energies): kappa = KB/2 (the code has KB H^2, the
# paper kappa/2 (c1+c2)^2 with H = (c1+c2)/2), w = KI (adhesion energy per area), sigma = KA, a = BeadRadius, so
#   w~ = 2 w a^2 / kappa = 4 KI r^2 / KB   and   sigma~ = sigma a^2 / kappa = 2 KA r^2 / KB,
# the degree of wrapping is z = 2 * (covered fraction of the 4 pi solid angle) = 2 * CoverageUnion,
# a/lambda = sqrt(sigma~) and w/sigma = w~ / (2 sigma~) = KI / KA (axes of Fig. 5).
EULER_GAMMA = 0.5772156649015329


def deserno_lines_small_tension(n=200):
    """Exact asymptotic E boundary (partially wrapped <-> fully enveloped), Eq. 26 of the paper, valid for
    sigma~ < 0.433 (w~ up to 4 + 8 exp(-1-2 gamma) = 4.928). Returns w~ and sigma~."""
    from scipy.special import lambertw
    w = np.linspace(4.0, 4.0 + 8 * np.exp(-1 - 2 * EULER_GAMMA), n)[1:]
    W = lambertw(-(w - 4) / 8 * np.exp(2 * EULER_GAMMA), -1).real
    sigma = (w - 4) / 4 * (1 + np.sqrt(1 + 1 / (2 * W) + 1 / (2 * W) ** 2))
    return w, sigma


def _deserno_axes_data(filepath, columns):
    """x = w~, y = sigma~ and the requested columns of one Coverage_data.txt of PostProcessing (planar)."""
    with open(filepath) as f:
        names = f.readline().lstrip("#").split()
    col = {n: i for i, n in enumerate(names[1:])}  # the DIR column is a string: usecols skips it
    data = np.loadtxt(filepath, comments="#", usecols=range(1, len(names)), ndmin=2)
    KA, KB, KI, r = (data[:, col[n]] for n in ("KA", "KB", "KI", "BeadRadius"))
    w_t = np.round(4 * KI * r**2 / KB, 6)
    sigma_t = np.round(2 * KA * r**2 / KB, 6)
    return w_t, sigma_t, {c: data[:, col[c]] for c in columns}


def planar_phase_deserno(filepath="../Results/Wrapping_planar_excess/Coverage_data.txt",
                         filepath_full="../Results/Wrapping_planar_excess/Coverage_data_full.txt",
                         column="CoverageUnion", as_z=True, regions=False, region_thresholds=(0.3, 0.95),
                         contour_levels=None, theory=True, xlim=(2.0, 10.0), ylim=(0.0, 2.0), save=None):
    """Planar membrane phase diagram in the axes of Fig. 2 of Deserno: w~ = 4 KI r^2/KB (x), sigma~ = 2 KA r^2/KB (y).

    filepath, filepath_full: the Coverage_data.txt of PostProcessing for the flat start and for the full (wrapped)
    start. When filepath_full is not None, every point takes the run with the lower Total_E of the two (None: only
    filepath). Compare Total_E only between runs with the same tension energy.
    column: coloured column, by default CoverageUnion; as_z multiplies the coverage columns by 2 so that the colour
    is the degree of wrapping z in [0, 2] of the paper (CoverageUnion and CoveredArea are fractions of 4 pi).
    regions: instead of the colour scale, three greys like Fig. 2: free / partially wrapped / fully enveloped,
    with the CoverageUnion thresholds region_thresholds (arbitrary, the soft adhesion shell gives the free
    state a coverage of about 0.2).
    theory: dashed vertical line W (w~ = 4), dotted line w~ = 4 + 2 sigma~ (zero energy of the enveloped state),
    the E line (solid) and the spinodals S1, S2 (short dashed) solved numerically from the shape equations
    (deserno_theory.py and deserno_table.npz next to this file; without them only the exact small tension
    E line of the paper, Eq. 26, sigma~ < 0.433, is drawn).
    contour_levels: number of (or list of values for) contour lines of the coloured column, None for no lines.
    """
    w_a, s_a, d_a = _deserno_axes_data(filepath, [column, "Total_E", "CoverageUnion"])
    if filepath_full is not None:
        w_b, s_b, d_b = _deserno_axes_data(filepath_full, [column, "Total_E", "CoverageUnion"])
    else:
        w_b, s_b, d_b = w_a, s_a, d_a

    # Sorted unique grid values over both files, the pivot fills [sigma~, w~]
    ws = np.unique(np.concatenate((w_a, w_b)))
    ss = np.unique(np.concatenate((s_a, s_b)))

    def grid(w, s_, values):
        G = np.full((len(ss), len(ws)), np.nan)
        G[np.searchsorted(ss, s_), np.searchsorted(ws, w)] = values
        return G

    E_a = np.where(np.isnan(grid(w_a, s_a, d_a["Total_E"])), np.inf, grid(w_a, s_a, d_a["Total_E"]))
    E_b = np.where(np.isnan(grid(w_b, s_b, d_b["Total_E"])), np.inf, grid(w_b, s_b, d_b["Total_E"]))
    full_wins = (E_b < E_a) if filepath_full is not None else np.zeros_like(E_a, dtype=bool)
    if filepath_full is not None:
        both = np.isfinite(E_a) & np.isfinite(E_b)
        print("The full start has the lower energy in {} of {} points".format((full_wins & both).sum(), both.sum()))

    def pick(name):
        return np.where(full_wins, grid(w_b, s_b, d_b[name]), grid(w_a, s_a, d_a[name]))

    Z = pick(column)
    cov_columns = ("CoverageUnion", "CoveredArea")
    label = column.replace("_", " ")
    if as_z and column in cov_columns:
        Z = 2 * Z
        label = r"$z$ (degree of wrapping)"

    # Cell edges halfway between the centres, so each cell is centred on its run
    def edges(c):
        if len(c) == 1:
            return np.array([c[0] - 0.5, c[0] + 0.5])
        mid = 0.5 * (c[1:] + c[:-1])
        return np.concatenate(([c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]))

    fig, ax = plt.subplots()
    if regions:
        from matplotlib.colors import ListedColormap
        cov = pick("CoverageUnion")
        R = np.where(np.isnan(cov), np.nan, np.where(cov < region_thresholds[0], 0.0, np.where(cov < region_thresholds[1], 1.0, 2.0)))
        ax.pcolormesh(edges(ws), edges(ss), R, cmap=ListedColormap(["#ffffff", "#e3e3e3", "#c9c9c9"]), vmin=0, vmax=2, shading="flat")
    else:
        mesh = ax.pcolormesh(edges(ws), edges(ss), Z, cmap="viridis", shading="flat")
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label(label, rotation=90)
        cbar.set_label(r" Coverage ", rotation=90,fontsize=25)

    # Contour lines through the cell centres, on top of the pcolormesh (NaN cells are left out)
    if contour_levels is not None:
        cs = ax.contour(ws, ss, np.ma.masked_invalid(Z), levels=contour_levels,
                        colors="black" if regions else "white", linewidths=0.8)
        ax.clabel(cs, fontsize="small", fmt="%.2f")

    if theory:
        ax.axvline(x=4.0, ls="dashed", color="black", lw=1.2)                       # W
        sig = np.linspace(0.0, ylim[1], 50)
        ax.plot(4.0 + 2.0 * sig, sig, ls="dotted", color="black", lw=0.9)           # enveloped state has zero energy
        try:
            import deserno_theory
            deserno_theory.plot_fig2_lines(ax, ylim=ylim)                           # E, S1, S2
        except ImportError:
            w_E, s_E = deserno_lines_small_tension()
            ax.plot(w_E, s_E, color="black", lw=1.6)                                # E, exact for sigma~ < 0.433

    ax.set_xlabel(r"$\tilde w = 2 w a^2/\kappa $",fontsize=25) # = 4 K_I r^2/K_B
    ax.set_ylabel(r"$\tilde\sigma = \sigma a^2/\kappa $",fontsize=25) #= 2 K_A r^2/K_B
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("auto")

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    plt.show()
    return fig, ax

def deserno_fig5(filepath="../Results/WrappingPhaseSpaceNewFast/CoverageData.txt", filepath_full=None,
                 points=None, column="CoverageUnion", as_z=True, kappa_over_KB=0.5, w_over_KI=1.0,
                 xlim=(0.8, 1.0e3), ylim=(0.0, 5.0), shade=True, save=None):
    """Fig. 5 of Deserno: influence of the particle radius on wrapping, w/sigma against a/lambda (log axis).

    Solid curve = E boundary (partially wrapped <-> fully enveloped), dashed = W boundary (free <-> partially
    wrapped, w/sigma = 2 (lambda/a)^2), both from deserno_theory (numerical shape equations, table next to it).
    Simulation points (default: the vesicle runs in WrappingPhaseSpaceNewFast/CoverageData.txt of PostProcessing,
    KA = 1 fixed and KB varied), coloured by z = 2 CoverageUnion when as_z:
        a/lambda = r sqrt(sigma/kappa) = r sqrt(KA / (kappa_over_KB KB)),   w/sigma = w_over_KI KI / KA
    kappa_over_KB = 0.5 because the code has KB H^2 with H = (c1+c2)/2 (the paper: kappa/2 (c1+c2)^2); w_over_KI = 1
    is exact for the Adhesion interaction (planar runs), but the vesicle runs of that file use the Frenkel
    potential of every vertex, whose adhesion per area is not KI: w_over_KI is then an unknown calibration.
    filepath_full: the full-start file of a planar batch, to keep the lower Total_E of the two starts (needs a
    Total_E column, the vesicle file has none).
    points: instead of a file, (a_over_lambda, w_over_sigma, color) arrays to draw, color can be None.
    """
    import deserno_theory as dt

    fig, ax = plt.subplots()
    a = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 400)
    w_E, w_W = dt.fig5_curves(a)
    if shade:
        ax.fill_between(a, w_E, ylim[1], color="0.72", lw=0)         # fully enveloped
        ax.fill_between(a, w_W, w_E, color="0.86", lw=0)             # partially wrapped
        ax.fill_between(a, 0, w_W, color="0.97", lw=0)               # free
    ax.plot(a, w_E, color="black", lw=2.0)
    ax.plot(a, w_W, color="black", lw=2.0, ls="dashed")
    ax.plot([40, xlim[1]], [2, 2], color="black", lw=1.0, ls="dashed")   # E tends to w/sigma = 2

    if points is None and filepath is not None:
        def load(path):
            with open(path) as f:
                names = f.readline().lstrip("#").split()
            col = {n: i for i, n in enumerate(names[1:])}
            data = np.loadtxt(path, comments="#", usecols=range(1, len(names)), ndmin=2)
            return {n: data[:, col[n]] for n in names[1:]}
        d = load(filepath)
        if filepath_full is not None:
            d2 = load(filepath_full)
            # same (KA, KI) in both files: keep the lower energy
            first = {(round(k, 6), round(i, 6)): n for n, (k, i) in enumerate(zip(d["KA"], d["KI"]))}
            for n, (k, i) in enumerate(zip(d2["KA"], d2["KI"])):
                m = first.get((round(k, 6), round(i, 6)))
                if m is None:
                    for key in d:
                        d[key] = np.append(d[key], d2[key][n])
                elif d2["Total_E"][n] < d["Total_E"][m]:
                    for key in d:
                        d[key][m] = d2[key][n]
        keep = d["KA"] > 0                                   # w/sigma = KI/KA is undefined at zero tension
        z = d[column][keep] * (2 if as_z and column in ("CoverageUnion", "CoveredArea") else 1)
        points = (d["BeadRadius"][keep] * np.sqrt(d["KA"][keep] / (kappa_over_KB * d["KB"][keep])),
                  w_over_KI * d["KI"][keep] / d["KA"][keep], z)
    if points is not None:
        x, y, c = points
        if c is None:
            ax.scatter(x, y, s=14, color="black", zorder=3)
        else:
            sc = ax.scatter(x, y, c=c, s=22, cmap="viridis", edgecolor="black", linewidth=0.3, zorder=3)
            cbar = fig.colorbar(sc, ax=ax)
            cbar.set_label(r"$z$ (degree of wrapping)" if as_z else column.replace("_", " "), rotation=90)

    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$a/\lambda$")
    ax.set_ylabel(r"$w/\sigma$")

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    plt.show()
    return fig, ax


def phase_one_vesicle(batch="../Results/Phase_one_vesicle", legs=("in", "out"), column="CoverageUnion", as_z=True,
                      kappa_over_KB=0.5, w_over_KI=1.0, xlim=(0.8, 1.0e3), ylim=(0.0, 5.0), shade=True, save=None, show=True):
    """The chains of Phase_one_vesicle (Mem_one_vesicle_up.sh: KB scanned at fixed KI, KA = 1) in the axes and scaling of deserno_fig5:
    w/sigma against a/lambda (log axis), the E and W curves of Deserno, coloured by z = 2 * CoverageUnion.
        a/lambda = r sqrt(KA / (kappa_over_KB KB)),   w/sigma = w_over_KI KI / KA      (as deserno_fig5; points with KA = 0 are left out)
    batch: the folder with the Coverage_data_<Tag>.csv of every chain (Coverage_data_chain.py, run by Collect_chain_results.sh on the cluster), or a glob of csv's.
    legs: which legs to draw. "in" (the way along the KB list) as circles, "out" (the way back) as smaller squares on top, so that both numbers of a point
    show where they differ; "cold" (the cold start, the same KB as the first link) is left out by default.
    column: the coloured column of the csv (CoverageUnion; CoveredArea, cap_cfg or f_ad are the others); as_z multiplies the coverage columns by 2 (z in [0, 2]).
    """
    import csv
    import glob
    import os
    import deserno_theory as dt

    files = sorted(glob.glob(os.path.join(batch, "Coverage_data_*.csv") if os.path.isdir(batch) else batch))
    if not files:
        raise FileNotFoundError("no Coverage_data_*.csv in " + str(batch))
    rows = []
    for fn in files:
        with open(fn) as f:
            rows += [r for r in csv.DictReader(f)]
    print("{} points from {} files".format(len(rows), len(files)))

    fig, ax = plt.subplots()
    a = np.logspace(np.log10(xlim[0]), np.log10(xlim[1]), 400)
    w_E, w_W = dt.fig5_curves(a)
    if shade:
        ax.fill_between(a, w_E, ylim[1], color="0.72", lw=0)         # fully enveloped
        ax.fill_between(a, w_W, w_E, color="0.86", lw=0)             # partially wrapped
        ax.fill_between(a, 0, w_W, color="0.97", lw=0)               # free
    ax.plot(a, w_E, color="black", lw=2.0)
    ax.plot(a, w_W, color="black", lw=2.0, ls="dashed")
    ax.plot([40, xlim[1]], [2, 2], color="black", lw=1.0, ls="dashed")   # E tends to w/sigma = 2

    scale = 2.0 if as_z and column in ("CoverageUnion", "CoveredArea") else 1.0
    sc = None
    for leg, marker, size in (("in", "o", 34), ("cold", "D", 20), ("out", "s", 14)):       # the way back is drawn last, smaller, on top
        if leg not in legs:
            continue
        sel = [r for r in rows if r["leg"] == leg and float(r["KA"]) > 0]
        if not sel:
            continue
        KA, KB, KI, r0 = (np.array([float(r[k]) for r in sel]) for k in ("KA", "KB", "KI", "BeadRadius"))
        x = r0 * np.sqrt(KA / (kappa_over_KB * KB))
        y = w_over_KI * KI / KA
        c = scale * np.array([float(r[column]) for r in sel])
        sc = ax.scatter(x, y, c=c, s=size, marker=marker, cmap="viridis", vmin=0.0 if scale == 2.0 else None, vmax=2.0 if scale == 2.0 else None,
                        edgecolor="black", linewidth=0.3, zorder=3, label="way " + leg if leg != "cold" else "cold start")
    if sc is not None:
        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label(r"$z$ (degree of wrapping)" if scale == 2.0 else column.replace("_", " "), rotation=90)
        if len(legs) > 1:
            ax.legend(loc="upper left", fontsize="small", framealpha=0.9)

    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$a/\lambda$")
    ax.set_ylabel(r"$w/\sigma$")

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    if show:
        plt.show()
    return fig, ax


def Shell_test_comparison(filepath="../Results/CoverageShellSum/CoverageSum.txt", save=None, deserno=True,
                          quantities=("coverage", "bending", "total")):
    """Bead coverage (CoverageUnion), bending energy and total energy vs adhesion strength w~ = 4 KI r^2/KB for the shell test.

    Columns follow the CoverageSum header: P_order (0), shell_width (1),
    KI (5), CoverageUnion (9), Bending_tan (14), Total_E (16). Rows are grouped by (p, s), sorted by KI and
    duplicates (same p, s, KI) are collapsed, so the lines are properly ordered. The x axis is w~ = 4 KI (KB = r = 1).
    One figure per fixed p (curves for each s) and one per fixed s (curves per p), for each quantity in `quantities`
    ("coverage", "bending", "total"; a single name as a string is fine); the figures are saved as <save>_<p1|s0.25>.png (coverage),
    <save>_bending_<...>.png and <save>_total_<...>.png.
    deserno=True adds Deserno's zero-tension result (the shell test is KA = 0, KB = r = 1) as a black line: nothing wrapped below
    w~ = 4 (KI = 1), fully wrapped above (coverage z/2 = 1, bending 8 pi kappa = 4 pi KB, total E = -pi (w~ - 4), a ramp).
    """
    from collections import defaultdict

    if isinstance(quantities, str):
        quantities = (quantities,)
    spec = {"coverage": dict(col=9, ylabel="Bead coverage", full=lambda w: 1.0, tag=""),
            "bending": dict(col=14, ylabel=r"Bending energy $E_{\rm bend}$", full=lambda w: 4 * np.pi, tag="bending_"),
            "total": dict(col=16, ylabel=r"Total energy $E$", full=lambda w: -np.pi * (w - 4), tag="total_")}
    curves = {q: defaultdict(dict) for q in quantities}          # quantity -> (p, s) -> {KI: value}
    for line in open(filepath):
        if line.startswith("#") or not line.strip():
            continue
        d = line.split()
        for q in quantities:
            curves[q][(int(d[0]), float(d[1]))][float(d[5])] = float(d[spec[q]["col"]])

    keys = curves[quantities[0]]
    ps = sorted({p for p, _ in keys})
    ss = sorted({s for _, s in keys}, reverse=True)

    def panel(q, varying, make_key, make_label, title, fname):
        fig, ax = plt.subplots(figsize=(5.5, 4))
        cmap = plt.get_cmap("inferno")
        markers = "osD^v"
        for i, v in enumerate(varying):
            key = make_key(v)
            if key in curves[q]:
                c = curves[q][key]
                ki = sorted(c)
                ax.plot([4 * k for k in ki], [c[k] for k in ki], "-", marker=markers[i % len(markers)], ms=4, lw=1.6,
                        color=cmap(i / max(len(varying) - 1, 1) * 0.85), label=make_label(v))
        if deserno:
            xl = ax.get_xlim()
            full = spec[q]["full"]
            ax.plot([xl[0], 4.0, 4.0, xl[1]], [0, 0, full(4.0), full(xl[1])], "k-", lw=1.3, label=r"Deserno, $\tilde\sigma = 0$")
            ax.set_xlim(xl)
        ax.set_xlabel(r"Adhesion strength $\tilde w = 4 K_I r^2/K_B$")
        ax.set_ylabel(spec[q]["ylabel"])
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend(frameon=False)
        fig.tight_layout()
        if save is not None:
            fig.savefig(f"{save}_{spec[q]['tag']}{fname}.png", dpi=200, bbox_inches="tight")
        plt.show()

    for q in quantities:
        for p in ps:
            panel(q, ss, lambda s: (p, s), lambda s: f"s = {s:g}", f"p = {p}", f"p{p}")
        for s in ss:
            panel(q, ps, lambda p: (p, s), lambda p: f"p = {p}", f"s = {s:g}", f"s{s:g}")


def Potentials_comparison(filepath="../Results/Potentials_sum/CoverageSum.txt", save=None,
                          conds=("cos", "frenkel", "lj25", "lj20", "morse2", "morse10"), deserno=True):
    """Planar wrapping at zero tension with other radial profiles of the adhesion (Process_potentials.py output).

    Columns are found by name in the header (Cond Init Profile shell_width shell_alpha shell_power DIR KA KB KI ...).
    One figure with four panels, against the adhesion strength KI (top axis: w~ = 4 KI r^2/KB):
    coverage (CoverageUnion), Total_E, Total_E(flat) - Total_E(full) (the lower one is the equilibrium state; 0 means
    no hysteresis) and the bending energy Bending_tan (the dashed line is 4 pi K_B, the energy of a fully wrapped sphere of
    radius 1 with H = 1). Lines are the flat start, hollow markers the full start. BeadX is not plotted: the free sheet drifts
    along x, flat and full starts with the same energy end at different BeadX.
    The dotted vertical line is Deserno's wrapping transition at zero tension (w~ = 4, KI = 1 for r = KB = 1).
    deserno=True adds Deserno's zero-tension result as a black step line: nothing is wrapped for w~ < 4 (coverage 0, E = 0), the
    bead is fully wrapped above (coverage 1 = z/2 with z = 2, E~ = -2 (w~ - 4), i.e. E = -4 pi (KI - 1) for r = KB = 1, bending
    8 pi kappa = 4 pi KB); checked against deserno_theory.py at sigma~ = 1e-4.
    Prints the KI where the coverage of each condition crosses 0.5 and 0.9.
    """
    from collections import defaultdict

    names, data = None, defaultdict(lambda: defaultdict(dict))   # (cond, init) -> column -> {KI: value}
    labels = {}
    for line in open(filepath):
        if not line.strip():
            continue
        if line.startswith("#"):
            names = line.lstrip("#").split()
            continue
        d = dict(zip(names, line.split()))
        key = (d["Cond"], d["Init"])
        ki = float(d["KI"])
        for col in ("CoverageUnion", "Total_E", "Bending_tan"):
            data[key][col][ki] = float(d[col])
        prof, sw, al = d["Profile"], float(d["shell_width"]), d["shell_alpha"]
        labels[d["Cond"]] = {"cosine": f"cosine, s = {sw:g}", "frenkel": f"Frenkel, rc = {1 + sw:.2f} r",
                             "lj": f"LJ, rc = {1 + sw:.2f} r", "morse": f"Morse a = {al}, rc = {1 + sw:.2f} r"}[prof]

    conds = [c for c in conds if (c, "flat") in data]
    # Cove categorical hues, fixed order
    colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
    markers = "osD^vP"

    def curve(cond, init, col):
        c = data[(cond, init)][col]
        ki = np.array(sorted(c))
        return ki, np.array([c[k] for k in ki])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    panels = [("CoverageUnion", "Coverage"), ("Total_E", r"Total energy $E$"),
              ("dE", r"$E_{\rm flat} - E_{\rm full}$"), ("Bending_tan", r"Bending energy $E_{\rm bend}$")]
    for ax, (col, ylabel) in zip(axes.ravel(), panels):
        for i, cond in enumerate(conds):
            color, marker = colors[i % len(colors)], markers[i % len(markers)]
            if col == "dE":
                ki, flat = curve(cond, "flat", "Total_E")
                if (cond, "full") in data:
                    ki_f, full = curve(cond, "full", "Total_E")
                    common = np.intersect1d(ki, ki_f)
                    ax.plot(common, [data[(cond, "flat")]["Total_E"][k] - data[(cond, "full")]["Total_E"][k] for k in common], "-",
                            marker=marker, ms=4, lw=1.4, color=color, label=labels.get(cond, cond))
                continue
            ki, y = curve(cond, "flat", col)
            ax.plot(ki, y, "-", marker=marker, ms=4, lw=1.6, color=color, label=labels.get(cond, cond))
            if (cond, "full") in data:
                ki_f, y_f = curve(cond, "full", col)
                ax.plot(ki_f, y_f, "none", marker=marker, ms=7, mfc="none", mec=color, mew=1.0)
        ax.axvline(1.0, color="gray", ls=":", lw=1.0)
        if deserno and col != "dE":
            kd = np.array([0.5, 1.0, 1.0, 3.0])
            step = {"CoverageUnion": [0, 0, 1, 1], "Total_E": [0, 0, 0, -8 * np.pi],
                    "Bending_tan": [0, 0, 4 * np.pi, 4 * np.pi]}[col]
            ax.plot(kd, step, color="black", lw=1.3, label=r"Deserno, $\tilde\sigma = 0$")
        if col == "dE":
            ax.axhline(0.0, color="black", lw=0.6)
        if col == "Bending_tan":
            ax.axhline(4 * np.pi, color="black", lw=0.8, ls="--")
            ax.text(0.98, 4 * np.pi, r"$4\pi K_B$", ha="right", va="bottom", fontsize=8, transform=ax.get_yaxis_transform())
        ax.set_xlabel(r"Adhesion strength $K_I$")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        top = ax.secondary_xaxis("top", functions=(lambda x: 4 * x, lambda x: x / 4))
        top.set_xlabel(r"$\tilde w = 4 K_I r^2/K_B$")
    axes[0, 0].legend(frameon=False, fontsize=8, loc="lower right")
    axes[0, 0].set_ylim(-0.02, 1.02)
    for ax in axes.ravel():
        ax.set_xlim(0.72, 2.48)
    fig.suptitle(r"Planar membrane, $K_A = 0$: lines = flat start, hollow markers = full start; black = Deserno ($\tilde w = 4$)", y=1.0, fontsize=10)
    fig.tight_layout()
    if save is not None:
        fig.savefig(f"{save}.png", dpi=200, bbox_inches="tight")
    plt.show()

    print("KI where the coverage (flat start) crosses 0.5 and 0.9:")
    for cond in conds:
        ki, y = curve(cond, "flat", "CoverageUnion")
        order = np.argsort(y)
        cross = [float(np.interp(level, y[order], ki[order])) for level in (0.5, 0.9)]
        print(f"  {cond:8s} {labels.get(cond, cond):26s} KI(0.5) = {cross[0]:.3f}   KI(0.9) = {cross[1]:.3f}")


def Thin_shell_vs_potentials(filepath_shell="../Results/CoverageShellSum/CoverageSum.txt",
                             filepath="../Results/Potentials_sum/CoverageSum.txt", p=1, s=0.05, save=None,
                             conds=("cos", "frenkel", "lj25", "lj20", "morse2", "morse10"), deserno=True):
    """The thin cosine shell of the shell test (power p, half width s; default p = 1, s = 0.05) against all the potentials of the
    planar potentials scan, all at zero tension (KA = 0, KB = r = 1): total energy, bending energy and coverage against KI.

    filepath_shell: CoverageSum.txt of the shell test (columns KI 5, CoverageUnion 9, Bending_tan 14, Total_E 16, a flat start);
    filepath: the output of Process_potentials.py, flat start (lines). The thin shell is the thick black line, Deserno's zero-tension
    step is gray. Total_E and the bending energy are comparable between all curves. CoverageUnion is NOT quite: it counts the faces
    inside the shell of the run (half width 0.05 for the thin shell), the potentials are measured with the plain 0.25 shell.
    """
    names, rows = None, {}                                    # cond -> {column: {KI: value}}
    labels = {}
    for line in open(filepath):
        if not line.strip():
            continue
        if line.startswith("#"):
            names = line.lstrip("#").split()
            continue
        d = dict(zip(names, line.split()))
        if d["Init"] != "flat" or d["Cond"] not in conds:
            continue
        ki = float(d["KI"])
        r = rows.setdefault(d["Cond"], {"CoverageUnion": {}, "Total_E": {}, "Bending_tan": {}})
        for col in r:
            r[col][ki] = float(d[col])
        prof, sw, al = d["Profile"], float(d["shell_width"]), d["shell_alpha"]
        labels[d["Cond"]] = {"cosine": f"cosine, s = {sw:g}", "frenkel": f"Frenkel, rc = {1 + sw:.2f} r",
                             "lj": f"LJ, rc = {1 + sw:.2f} r", "morse": f"Morse a = {al}, rc = {1 + sw:.2f} r"}[prof]
    thin = {"CoverageUnion": {}, "Total_E": {}, "Bending_tan": {}}
    for line in open(filepath_shell):
        if line.startswith("#") or not line.strip():
            continue
        d = line.split()
        if int(d[0]) == p and abs(float(d[1]) - s) < 1e-9:
            ki = float(d[5])
            thin["CoverageUnion"][ki], thin["Bending_tan"][ki], thin["Total_E"][ki] = float(d[9]), float(d[14]), float(d[16])
    if not thin["Total_E"]:
        raise ValueError(f"no rows with p = {p}, s = {s} in {filepath_shell}")

    colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
    markers = "osD^vP"
    panels = [("CoverageUnion", "Coverage (CoverageUnion)", 1.0), ("Total_E", r"Total energy $E$", None),
              ("Bending_tan", r"Bending energy $E_{\rm bend}$", 4 * np.pi)]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (col, ylabel, full) in zip(axes, panels):
        if deserno:
            kd = np.array([0.4, 1.0, 1.0, 2.6])
            y = {"CoverageUnion": [0, 0, 1, 1], "Total_E": [0, 0, 0, -4 * np.pi * 1.6], "Bending_tan": [0, 0, 4 * np.pi, 4 * np.pi]}[col]
            ax.plot(kd, y, color="gray", lw=1.5, label=r"Deserno, $\tilde\sigma = 0$")
        for i, cond in enumerate([c for c in conds if c in rows]):
            ki = np.array(sorted(rows[cond][col]))
            ax.plot(ki, [rows[cond][col][k] for k in ki], "-", marker=markers[i % len(markers)], ms=4, lw=1.3,
                    color=colors[i % len(colors)], label=labels[cond])
        ki = np.array(sorted(thin[col]))
        ax.plot(ki, [thin[col][k] for k in ki], "-", marker="o", ms=5, lw=2.4, color="black", label=f"cosine, p = {p:g}, s = {s:g} (shell test)")
        ax.axvline(1.0, color="gray", ls=":", lw=1.0)
        ax.set_xlim(0.45, 2.45)
        ax.set_xlabel(r"Adhesion strength $K_I$")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        top = ax.secondary_xaxis("top", functions=(lambda x: 4 * x, lambda x: x / 4))
        top.set_xlabel(r"$\tilde w = 4 K_I r^2/K_B$")
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    axes[1].set_ylim(-19.5, 1.0)
    fig.suptitle(rf"Planar membrane, $K_A = 0$: thin cosine shell (p = {p:g}, s = {s:g}) against the other potentials (flat start)", y=1.02, fontsize=10)
    fig.tight_layout()
    if save is not None:
        fig.savefig(f"{save}.png", dpi=200, bbox_inches="tight")
    plt.show()


def Softshell_prediction(filepath="../Results/Potentials_sum/CoverageSum.txt",
                         prediction="../Studies/deserno_shell/data/softshell_sigma0.csv", cond="cos", save=None):
    """Continuum soft-shell prediction against the simulation at zero tension, cosine shell (the 'cos' condition, s = 0.25).

    prediction: ../Studies/deserno_shell/theory/softshell_sigma0.py writes it (the constrained minimum e(G) of softshell_direct.py
    for sigma~ = 0, s = 0.25, p = 1, R = 10: columns G wt e Ebend Aex, a = 1, kappa = KB/2, energies in pi kappa). At w~ = 4 KI
    the state with wt(G) = w~ has coverage G, E~ = e - w~ G and bending Ebend; Total_E = (pi KB/2) E~ and E_bend = (pi KB/2) Ebend
    (KB = r = 1). The simulation's G is the adhesion-weighted covered solid angle, G = -Bead/(2 pi KI) (not CoverageUnion, which counts
    the directions hidden by any selected face). Simulation: flat start filled markers, full start hollow; prediction red, Deserno's
    zero-tension step black. The fourth panel is Total_E(simulation) - Total_E(prediction).
    """
    names, sim = None, {"flat": [], "full": []}
    for line in open(filepath):
        if not line.strip():
            continue
        if line.startswith("#"):
            names = line.lstrip("#").split()
            continue
        d = dict(zip(names, line.split()))
        if d["Cond"] == cond:
            ki = float(d["KI"])
            sim[d["Init"]].append((ki, -float(d["Bead"]) / (2 * np.pi * ki), float(d["Total_E"]), float(d["Bending_tan"])))
    sim = {k: np.array(sorted(v)).T for k, v in sim.items() if v}
    G, wt, e, Eb = np.loadtxt(prediction, comments="#", usecols=(0, 1, 2, 3), unpack=True)

    def predicted(ki):
        """coverage G, total energy and bending energy of the stationary state at w~ = 4 KI (wt(G) is monotonic at zero tension)"""
        w = 4 * np.asarray(ki, float)
        Gp = np.interp(w, wt, G)
        return Gp, 0.5 * np.pi * (np.interp(w, wt, e) - w * Gp), 0.5 * np.pi * np.interp(w, wt, Eb)

    ki_line = np.linspace(max(0.72, wt.min() / 4), min(2.48, wt.max() / 4), 300)
    Gp, Ep, Ebp = predicted(ki_line)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax_G, ax_E, ax_B, ax_R = axes.ravel()
    kd = np.array([0.5, 1.0, 1.0, 3.0])
    for ax, y_pred, y_des, ylabel, col in ((ax_G, Gp, [0, 0, 2, 2], r"Adhesion-weighted coverage $G = -E_{\rm adh}/(2\pi K_I)$", 1),
                                           (ax_E, Ep, [0, 0, 0, -8 * np.pi], r"Total energy $E$", 2),
                                           (ax_B, Ebp, [0, 0, 4 * np.pi, 4 * np.pi], r"Bending energy $E_{\rm bend}$", 3)):
        ax.plot(kd, y_des, color="black", lw=1.3, label=r"Deserno, $\tilde\sigma = 0$")
        ax.plot(ki_line, y_pred, color="#d03b3b", lw=2.0, label="continuum soft shell, $s = 0.25$")
        for init, ms, mfc in (("flat", 5, "#2a78d6"), ("full", 8, "none")):
            if init in sim:
                ax.plot(sim[init][0], sim[init][col], "none", marker="o", ms=ms, mfc=mfc, mec="#2a78d6", label=f"simulation, {init} start")
        ax.set_ylabel(ylabel)
    for init, ms, mfc in (("flat", 5, "#2a78d6"), ("full", 8, "none")):
        if init in sim:
            ki = sim[init][0]
            ax_R.plot(ki, sim[init][2] - predicted(ki)[1], "o", ms=ms, mfc=mfc, mec="#2a78d6", label=f"{init} start")
    ax_R.axhline(0, color="black", lw=0.6)
    ax_R.set_ylabel(r"$E_{\rm sim} - E_{\rm soft\,shell}$")
    ax_R.legend(frameon=False, fontsize=8)
    for ax in axes.ravel():
        ax.axvline(1.0, color="gray", ls=":", lw=1.0)
        ax.set_xlabel(r"Adhesion strength $K_I$")
        ax.set_xlim(0.72, 2.48)
        ax.grid(alpha=0.3)
        top = ax.secondary_xaxis("top", functions=(lambda x: 4 * x, lambda x: x / 4))
        top.set_xlabel(r"$\tilde w = 4 K_I r^2/K_B$")
    ax_G.set_ylim(-0.05, 2.1)
    ax_G.legend(frameon=False, fontsize=8, loc="lower right")
    fig.suptitle(r"Planar membrane, $K_A = 0$, cosine shell $s = 0.25$: simulation vs continuum prediction", y=1.0, fontsize=10)
    fig.tight_layout()
    if save is not None:
        fig.savefig(f"{save}.png", dpi=200, bbox_inches="tight")
    plt.show()

def _profile_spec(name):
    """Parameters of an adhesion profile from its name. Named variants of the runs:
        cos (cosine p = 1, s = 0.25), frenkel (s = 0.25), lj25 / lj20 (rc = 2.5 / 2 sigma_LJ), morse2 / morse10 (a = 2 / 10, rc = 2.5 sigma_LJ)
    or any `kind_tokens`, kind in cos, frenkel, lj, morse and tokens p<exponent>, s<half width or cutoff width>, a<Morse alpha>,
    e.g. cos_p4_s0.05, frenkel_s0.5, lj_s1, morse_a3_s1.2 (cos: p = 1, s = 0.25; others s = 0.25; morse a = 6 when not given).
    s is the shell_width of the input file: the cosine shell reaches out to r = 1 +- s, the potentials are cut at rc = 1 + s (radius = 1)."""
    s_lj25, s_lj20 = 2.5 * 2 ** (-1 / 6) - 1, 2.0 * 2 ** (-1 / 6) - 1
    named = {"cos": dict(kind="cos", p=1.0, s=0.25), "frenkel": dict(kind="frenkel", s=0.25),
             "lj25": dict(kind="lj", s=s_lj25), "lj20": dict(kind="lj", s=s_lj20),
             "morse2": dict(kind="morse", a=2.0, s=s_lj25), "morse10": dict(kind="morse", a=10.0, s=s_lj25)}
    if name in named:
        return dict(dict(p=1.0, s=0.25, a=6.0, named=True), **named[name])
    kind, *tokens = name.split("_")
    if kind not in ("cos", "frenkel", "lj", "morse"):
        raise ValueError(f"unknown potential '{name}': use {sorted(named)}, a group (cosine, scan, all) or kind_p<p>_s<s>_a<a> with kind in cos, frenkel, lj, morse")
    spec = dict(kind=kind, p=1.0, s=0.25, a=6.0, named=False)
    for t in tokens:
        if t[:1] not in ("p", "s", "a") or not t[1:]:
            raise ValueError(f"cannot read '{t}' in '{name}' (expected p<exponent>, s<width> or a<alpha>)")
        spec[{"p": "p", "s": "s", "a": "a"}[t[0]]] = float(t[1:])
    return spec


def _profile_weight(r, spec):
    """w(r) of the Adhesion weight (radius = 1, w = 1 at the minimum), the formulas of bead_geometry::shellProfileWeight."""
    r = np.asarray(r, float)
    kind, sw = spec["kind"], spec["s"]
    if kind == "cos":
        x = (r - 1.0) / sw
        return np.where(np.abs(x) < 1.0, (0.5 * (1.0 + np.cos(np.pi * np.clip(x, -1, 1)))) ** spec["p"], 0.0)
    rc = 1.0 + sw
    if kind == "frenkel":
        s2 = 2.0 * rc ** 2 / (3.0 * rc ** 2 - 1.0)                  # sigma_F^2, so that the minimum is at r = 1
        alpha = 2.0 * (rc ** 2 / s2) * (3.0 / (2.0 * (rc ** 2 / s2 - 1.0))) ** 3
        w = -alpha * (s2 / r ** 2 - 1.0) * (rc ** 2 / r ** 2 - 1.0) ** 2
    else:
        if kind == "lj":
            V = lambda x: x ** -12.0 - 2.0 * x ** -6.0
        else:
            V = lambda x: (1.0 - np.exp(-spec["a"] * (x - 1.0))) ** 2 - 1.0
        w = (V(rc) - V(r)) / (1.0 + V(rc))                          # shifted to zero at rc, depth renormalised to 1
    return np.where(r < rc, w, 0.0)


def Potential_data(tag, x=None, n=2000):
    """The radial profile of the adhesion between the bead and the membrane, as numpy arrays [x, U(x)].

    tag: `kind_p<p>_s<s>_a<a>` as in the input files, kind in cos, frenkel, lj, morse, for example cos_p1_s0.25, cos_p4_s0.05,
         frenkel_s0.25, lj_s1.2272468, morse_a2_s1.2272468; or one of the named variants of the runs: cos, frenkel, lj25, lj20,
         morse2, morse10 (see _profile_spec). p is the exponent of the cosine shell, s the shell_width (the cosine shell reaches
         r = 1 +- s, the potentials are cut at rc = 1 + s, LJ and Morse also shifted to zero there), a the Morse range alpha.
    x: where to evaluate, in units of the bead radius (default: n points from 0.3, or from just below the cosine shell, to 1.08 s
       beyond r = 1: the whole support).
    Returns x (bead radius = 1) and U(x) = -w(x) in units of W = inter_str: -1 at x = 1 (the minimum of the well), 0 beyond the cutoff,
    positive in the repulsive core of the Frenkel, LJ and Morse (the cosine shell has none). The adhesion energy density is W U(x).
    """
    spec = _profile_spec(tag)
    if x is None:
        lo = max(0.0, 1.0 - 1.08 * spec["s"]) if spec["kind"] == "cos" else 0.3
        x = np.linspace(max(lo, 1e-3), 1.0 + 1.08 * spec["s"], n)
    else:
        x = np.asarray(x, float)
    return x, 0.0 - _profile_weight(x, spec)


def potential_curve(tag, init="flat", filepath="../Results/Potentials_sum/CoverageSum.txt"):
    """x, Cov(x), Total_E(x) of one potential of the planar potentials scan (Process_potentials.py output).

    tag: the Cond column of the file ("cos", "frenkel", "lj25", "lj20", "morse2", "morse10").
    init: "flat" or "full", the start of the run (the Init column).
    x is the adhesion strength KI, Cov is CoverageUnion; the three arrays are sorted by KI, so
    x, cov, E = potential_curve("cos") can be plotted directly.
    """
    names, rows, tags = None, [], set()
    for line in open(filepath):
        if not line.strip():
            continue
        if line.startswith("#"):
            names = line.lstrip("#").split()
            continue
        d = dict(zip(names, line.split()))
        tags.add(d["Cond"])
        if d["Cond"] == tag and d["Init"] == init:
            rows.append((float(d["KI"]), float(d["CoverageUnion"]), float(d["Total_E"])))
    if not rows:
        raise ValueError(f"no rows with Cond = {tag!r}, Init = {init!r} in {filepath} (tags: {sorted(tags)})")
    x, cov, E = np.array(sorted(rows)).T
    return x, cov, E



def CosPlots():
    """Profile U(r)/K_I of the cosine adhesion shell, drawn like Shell_test_comparison (5.5 x 4 figure, grid, legend without
    frame). First figure: p = 1 with the three widths s, second: s = 0.1 with the three exponents p; both share the x and y
    window. No markers: the profile is analytic, not data points."""
    cmap = plt.get_cmap("inferno")
    xlim, ylim = (0.7, 1.3), (-1.1, 0.1)

    def panel(tags, labels, title, save):
        fig, ax = plt.subplots(figsize=(5.5, 4))
        ax.axhline(0,ls='solid',color='black')
        ax.axhline(-1,ls='dashed',color='grey')
        for tag, label, shade in zip(tags, labels, (0.75, 0.5, 0.25)):
            x,U = Potential_data(tag)
            ax.plot(x,U,color=cmap(shade),lw=1.6,label=label)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel(r"$r/a$")
        ax.set_ylabel(r"$U(r)/K_I$")
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(save, dpi=200, bbox_inches="tight")
        plt.show()

    panel(("cos_p1_s0.25", "cos_p1_s0.1", "cos_p1_s0.05"), ("s = 0.25", "s = 0.1", "s = 0.05"),
          "p = 1", "../Results/SoftHour/CosPot1.png")
    panel(("cos_p1_s0.1", "cos_p2_s0.1", "cos_p4_s0.1"), ("p = 1", "p = 2", "p = 4"),
          "s = 0.1", "../Results/SoftHour/CosPot2.png")


# CosPlots()

def SoftHourPlots():
    # OK lets see
    # plt.scatter(1.0,-1.0,marker='*',color='black')
    x, U = Potential_data("lj20")
    plt.xlim(0.5,2.25)
    plt.axhline(0,ls='solid',color='black')
    plt.axhline(-1,ls='dashed',color='grey')
    plt.ylim(-2,2)
    # plt.plot(x,U,label=r'LJ 2.0',color='purple')
    x, U = Potential_data("lj25")
    plt.plot(x,U,label=r'LJ 2.5',color='purple')
    
    # plt.xlabel(r"$r/a$",fontsize=20)
    # plt.legend(fontsize=15)
    # plt.ylabel(r"$U(r)/K_I$",fontsize=20)
    # plt.savefig("../Results/SoftHour/LJs.png",bbox_inches='tight')
    # plt.show()

    x, U = Potential_data("morse2")
    # plt.xlim(0.5,2.25)
    plt.axhline(0,ls='solid',color='black')
    plt.axhline(-1,ls='dashed',color='grey')
    # plt.ylim(-2,2)
    # plt.plot(x,U,label=r'Morse 2',color='purple')
    x, U = Potential_data("morse10")
    plt.plot(x,U,label=r'Morse 10',color='magenta')
    
    plt.xlabel(r"$r/a$",fontsize=20)
    plt.legend(fontsize=15)
    plt.ylabel(r"$U(r)/K_I$",fontsize=20)
    # plt.savefig("../Results/SoftHour/Morses.png",bbox_inches='tight')
    # plt.show()

    # x, U = Potential_data("morse2")
    # plt.xlim(0.5,2.25)
    # plt.axhline(0,ls='solid',color='black')
    # plt.axhline(-1,ls='dashed',color='grey')
    # plt.ylim(-2,2)
    # plt.plot(x,U,label=r'Frenkel',color='purple')
    x, U = Potential_data("frenkel")
    plt.plot(x,U,label=r'Frenkel',color='pink')
    
    plt.xlabel(r"$r/a$",fontsize=20)
    plt.legend(fontsize=15)
    plt.ylabel(r"$U(r)/K_I$",fontsize=20)
    # plt.scatter(1.0,-1.0,color='black')
    plt.savefig("../Results/SoftHour/Potentials.png",bbox_inches='tight')
    
    plt.show()


def SoftHourCov():
    """Zero-tension theory and the Morse / Frenkel / LJ potentials against w~ = 4 KI, drawn like Shell_test_comparison
    (5.5 x 4 figure, grid, markers, same axis labels, legend without frame). Solid lines: full start, dashed lines: flat start.
    All four figures share the x window; the two coverage figures share their y window and so do the two energy figures."""
    kd = np.array([0.0, 4.0, 4.0, 2.6*4])
    yCov = [0, 0, 1, 1]
    yEnergy = [0, 0, 0, -4 * np.pi * 1.6]
    xlabel = r"Adhesion strength $\tilde w = 4 K_I r^2/K_B$"
    theory = r"Deserno, $\tilde\sigma = 0$"
    xlim = (2.0, 7.2)                                       # KI = 0.5 ... 1.8 of the shell test, in w~ = 4 KI
    ylimCov = (-0.05, 1.05)
    ylimE = (-11.6, 0.6)                                    # the lowest energy in the window is -11.1 (LJ)
    style = dict(lw=1.6, ms=4)                              # as the curves of Shell_test_comparison

    def finish(fig, ax, ylabel, ylim, save):
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(save, dpi=200, bbox_inches="tight")
        plt.show()

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(kd, yCov, color="gray", lw=1.5, label=theory)
    finish(fig, ax, "Bead coverage", ylimCov, "../Results/SoftHour/TheoryCov.png")

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(kd, yEnergy, color="gray", lw=1.5, label=theory)
    finish(fig, ax, r"Total energy $E$", ylimE, "../Results/SoftHour/TheoryE.png")

    uMorse , covMorse, eMorse = potential_curve("morse10",init='full')
    uFrenkel, covFrenkel,eFrenkel = potential_curve("frenkel",init='full')
    uLJ, covLJ, eLJ = potential_curve("lj25",init='full')

    uMorse2 , covMorse2, eMorse2 = potential_curve("morse10",init='flat')
    uFrenkel2, covFrenkel2,eFrenkel2 = potential_curve("frenkel",init='flat')
    uLJ2, covLJ2, eLJ2 = potential_curve("lj25",init='flat')


    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(kd, yCov, color="gray", lw=1.5, label=theory)
    ax.plot(uMorse*4,covMorse,color="magenta",marker='o',label="Morse",**style)
    ax.plot(uFrenkel*4,covFrenkel,color='pink',marker='s',label="Frenkel",**style)
    ax.plot(uLJ*4,covLJ,color='purple',marker='D',label="LJ",**style)

    ax.plot(uMorse2*4,covMorse2,color="magenta",marker='o',ls='dashed',**style)
    ax.plot(uFrenkel2*4,covFrenkel2,color='pink',marker='s',ls='dashed',**style)
    ax.plot(uLJ2*4,covLJ2,color='purple',marker='D',ls='dashed',**style)
    finish(fig, ax, "Bead coverage", ylimCov, "../Results/SoftHour/TheoryCovPotentials.png")


    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(uMorse*4,eMorse,color="magenta",marker='o',label="Morse",**style)
    ax.plot(uFrenkel*4,eFrenkel,color='pink',marker='s',label="Frenkel",**style)
    ax.plot(uLJ*4,eLJ,color='purple',marker='D',label="LJ",**style)

    ax.plot(uMorse2*4,eMorse2,color="magenta",marker='o',ls='dashed',**style)
    ax.plot(uFrenkel2*4,eFrenkel2,color='pink',marker='s',ls='dashed',**style)
    ax.plot(uLJ2*4,eLJ2,color='purple',marker='D',ls='dashed',**style)

    ax.plot(kd, yEnergy, color="gray", lw=1.5, label=theory)
    finish(fig, ax, r"Total energy $E$", ylimE, "../Results/SoftHour/TheoryEPotentials.png")


    # Good now we add the other potentials.



        
def shell_curve(p, s, filepath="../Results/CoverageShellSum/CoverageSum.txt"):
    """x, Cov(x), Total_E(x) of one cosine shell (power p, half width s) of the shell test, like potential_curve (flat start only).

    x is the adhesion strength KI, Cov is CoverageUnion; the arrays are sorted by KI."""
    names, rows = None, []
    for line in open(filepath):
        if not line.strip():
            continue
        if line.startswith("#"):
            names = line.lstrip("#").split()
            continue
        d = dict(zip(names, line.split()))
        if int(d["P_order"]) == p and abs(float(d["shell_width"]) - s) < 1e-9:
            rows.append((float(d["KI"]), float(d["CoverageUnion"]), float(d["Total_E"])))
    if not rows:
        raise ValueError(f"no rows with p = {p}, s = {s} in {filepath}")
    x, cov, E = np.array(sorted(rows)).T
    return x, cov, E


def CosPotentialsPlots():
    """Coverage and total energy against w~ = 4 KI of two cosine shells (p = 1, s = 0.25 and p = 2, s = 0.05, from the shell test)
    and of the Morse / Frenkel / LJ potentials, drawn like SoftHourCov (same colors, limits and style). All curves are the flat
    start, the only one the shell test has."""
    cmap = plt.get_cmap("inferno")
    kd = np.array([0.0, 4.0, 4.0, 2.6*4])
    yCov = [0, 0, 1, 1]
    yEnergy = [0, 0, 0, -4 * np.pi * 1.6]
    xlabel = r"Adhesion strength $\tilde w = 4 K_I r^2/K_B$"
    theory = r"Deserno, $\tilde\sigma = 0$"
    xlim = (2.0, 7.2)                                       # same windows as SoftHourCov
    style = dict(lw=1.6, ms=4)                              # as the curves of Shell_test_comparison

    curves = [("Cosine, p = 1, s = 0.25", shell_curve(1, 0.25), cmap(0.75), '^'),
              ("Cosine, p = 2, s = 0.05", shell_curve(2, 0.05), cmap(0.1), 'v'),
              ("Morse", potential_curve("morse10"), "magenta", 'o'),
              ("Frenkel", potential_curve("frenkel"), "pink", 's'),
              ("LJ", potential_curve("lj25"), "purple", 'D')]       # (x, cov, E) each

    for ylabel, col, y_theory, ylim, save in (("Bead coverage", 1, yCov, (-0.05, 1.05), "../Results/SoftHour/CosPotentialsCov.png"),
                                              (r"Total energy $E$", 2, yEnergy, (-11.6, 0.6), "../Results/SoftHour/CosPotentialsE.png")):
        fig, ax = plt.subplots(figsize=(5.5, 4))
        ax.plot(kd, y_theory, color="gray", lw=1.5, label=theory)
        for label, curve, color, marker in curves:
            ax.plot(curve[0]*4, curve[col], color=color, marker=marker, label=label, **style)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(saave, dpi=200, bbox_inches="tight")
        plt.show()


def Plot_data_foce(directory):
    # Ok so he datafile would be
    Data = np.loadtxt(directory+"Force_data.csv")
    Distance_x= []
    Force_x = []
    Distance_x_2 = []
    Force_x_2 = []
    for i in range(len(Data)):
        # Ok so thi is gonna go over whats in data
        # print(Data[i,0])
        if(Data[i,2]==0.3 ):
            # print(Data[i,3])
            if(Data[i,1] == 16):
                Distance_x.append(Data[i,3])
                Force_x.append(-1*Data[i,4])
            else:
                Distance_x_2.append(Data[i,3])
                Force_x_2.append(-1*Data[i,4])
                

    plt.axhline(2*np.pi*np.sqrt(16*4),ls='dashed',color='purple')
    plt.axhline(2*np.pi*np.sqrt(32*4),ls='dashed',color='magenta')
    plt.scatter(np.array(Distance_x)+0.3,Force_x,color='purple',label=r'$\kappa=8$')
    plt.xlabel("Membrane extension")
    plt.ylabel("Membrane response to pulling")
    # plt.axvline(2.1)
    # plt.axvline(2.3)
    plt.scatter(np.array(Distance_x_2)+0.3,Force_x_2,color='magenta',label=r'$\kappa=16$')
    plt.ylim(0,150)
    # print(Distance_x)
    # print(Force_x)
    plt.title(r'$\sigma = 4$')
    plt.legend( )
    plt.savefig(directory + "Scatter_plot.png",bbox_inches='tight')
    plt.show()

# load_pretty_figure_setup()
# SoftHourPlots()
# SoftHourCov()
# CosPotentialsPlots()
# Plot_data_foce("../Results/Pulling_and_relaxing_2/")
# phase_one_vesicle(batch="../Results/Phase_one_vesicle/",legs="out",as_z=False)
phase_one_vesicle_map(batch="../Results/Phase_one_vesicle/")

# Shell_test_comparison(filepath="../Results/CoverageShellSum/CoverageSum.txt",save="../Results/CoverageShellSum/CoverageSum2",quantities="total")

# Potentials_comparison(save="../Results/Potentials_sum/Potentials_comparison")                                      # Frenkel / LJ / Morse adhesion, KA = 0
# Thin_shell_vs_potentials(save="../Results/Potentials_sum/Thin_shell_vs_potentials")                                 # cosine p = 1, s = 0.05 (shell test) against all the potentials
# Softshell_prediction(save="../Results/Potentials_sum/Softshell_prediction")                                           # cosine shell vs continuum model, needs Studies/deserno_shell/theory/softshell_sigma0.py
# x, U = Potential_data("cos_p1_s0.25")                                                                         # radial profile [x, U(x)] of an adhesion variant, x in bead radii, U in units of W
# deserno_fig5(filepath = "../Results/WrappingPhaseSpaceNewFast/Coverage_data.txt",as_z = False)                         # vesicle data, kappa = KB/2
# deserno_fig5(filepath="../Results/Wrapping_planar_excess/Coverage_data.txt", filepath_full="../Results/Wrapping_planar_excess/Coverage_data_full.txt")
# planar_phase_deserno(filepath="../Results/Wrapping_planar_excess/Coverage_data_thin.txt",filepath_full="../Results/Wrapping_planar_excess/Coverage_data_thin.txt",as_z = False,save="../Results/Wrapping_planar_excess/Planar_phasepace_thin.png",theory=True)
# planar_phase_deserno(as_z = False,save="../Results/Wrapping_planar_excess/Planar_phasepace_theory.png",theory=True)
# planar_phase_deserno(regions=True)
# planar_phase_deserno(filepath_full="../Results/Wrapping_planar_excess/Coverage_data_thin.txt", filepath="../Results/Wrapping_planar_excess/Coverage_data_thin.txt",as_z=False)
# planar_phase(contour_levels=[ 0.5 ,0.75],filepath="../Results/Wrapping_planar_excess/Coverage_data.txt")
