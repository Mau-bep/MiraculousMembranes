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
    filepath = "../Results/Wrapping_PhaseLogscale/Coverage_data.txt"

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

    plt.scatter(X,Y,c=Data[:,5])
    plt.xscale('log')

    # plt.xlim(0.)
    plt.axvline(x=4.4,ls='dashed',color='black')
    plt.axhline(y=1.37,ls='dashed',color='black')
    plt.colorbar()
    plt.show()
    

    print("THe values of X are {}".format(X))
    print("The values of Y are {}".format(Y))
    
    return

contour()