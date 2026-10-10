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
    cbar.set_label('Wrapping fraction', rotation=90)
    plt.savefig(folder + "PhaseSpaceVesicle.pdf",bbox_inches = 'tight')
    plt.show()
    

    print("THe values of X are {}".format(X))
    print("The values of Y are {}".format(Y))
    
    return

# contour()

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
                         contour_levels=None, theory=True, xlim=(3.0, 7.0), ylim=(0.0, 1.0), save=None):
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

    ax.set_xlabel(r"$\tilde w = 2 w a^2/\kappa = 4 K_I r^2/K_B$")
    ax.set_ylabel(r"$\tilde\sigma = \sigma a^2/\kappa = 2 K_A r^2/K_B$")
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

# deserno_fig5(filepath = "../Results/WrappingPhaseSpaceNewFast/Coverage_data.txt",as_z = False)                         # vesicle data, kappa = KB/2
# deserno_fig5(filepath="../Results/Wrapping_planar_excess/Coverage_data.txt", filepath_full="../Results/Wrapping_planar_excess/Coverage_data_full.txt")
# planar_phase_deserno(contour_levels=[0.5, 1.0, 1.5])
# planar_phase_deserno(regions=True)
planar_phase_deserno(filepath="../Results/Wrapping_planar_excess/Coverage_data_thin.txt",as_z=False)
# planar_phase(contour_levels=[ 0.5 ,0.75],filepath="../Results/Wrapping_planar_excess/Coverage_data.txt")
