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

def bending_phase(filepath="../Results/TwoBeadsCov/Bending_data_phase.txt", column="bending", save=None):
    """pcolormesh of the final energy over the (bead distance, target coverage) phase space.

    File columns: run, distance, target coverage, bending energy, coverage energy, coverage strength.
    column: "bending" or "coverage", selects the energy that is coloured.
    """
    columns = {
        "bending": (3, r"$E_\mathrm{bend}$"),
        "coverage": (4, r"$E_\mathrm{cov}$"),
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

    ax.set_xlabel(r"Bead distance $d$")
    ax.set_ylabel(r"Target coverage")
    ax.set_xticks(xs)
    ax.set_yticks(ys)
    ax.set_aspect("auto")

    if save is not None:
        fig.savefig(save, bbox_inches="tight")
    plt.show()
    return fig, ax

bending_phase(filepath="../Results/TwoBeadsCov/Bending_data_phase_Bendi.txt", column="bending", save="../Results/TwoBeadsCov/Bending_phase_Bendi.pdf")


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
