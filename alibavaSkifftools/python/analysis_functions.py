#!/usr/bin/env python
"""Pool of functions for the sensor analysis
"""
__author__ = "Jordi Duarte-Campderros"
__credits__ = ["Jordi Duarte-Campderros"]
__version__ = "0.1"
__maintainer__ = "Jordi Duarte-Campderros"
__email__ = "jorge.duarte.campderros@cern.ch"
__status__ = "Development"

## TODO: 
## - remove the OptStat
## - a new canvas with the time fit extraction, time window used
##   and statistics

def spectrum_searcher_landau_parameters(h,xmin,xmax):
    """Search the peak for the Landau using the TSpectrum

    Parameters
    ----------
    h: ROOT.TH1
        The histogram
    xmin: float
        The minimum value the peak should be
    xmax: float
        The maximum value the peak should be

    Return
    ------
    peak: float
    """
    import ROOT 
    # Uset TSpectrum to find the available peaks
    batch=ROOT.gROOT.IsBatch()
    ROOT.gROOT.SetBatch(1)
    s=ROOT.TSpectrum()
    nfound = s.Search(h,1,'nobackground new')
    # Get the first peak comming from the left
    for i in range(nfound-1,-1,-1):
        peak = s.GetPositionX()[i]
        if abs(xmin) > abs(peak) or abs(peak) > abs(xmax):
            continue
        break
    ROOT.gROOT.SetBatch(batch)

    return peak

def landau_gaus(x,par):
    """Definition of a Landau convulated with gaus 

    Based on 
    https://root.cern.ch/root/html/tutorials/fit/langaus.C.html

    Parameters
    ----------
    x: array
        The variable
    par: array
        The parameters [MPV,Sigma landau, Normalization, Sigma Gauss]
    """
    from math import pi,sqrt
    import ROOT 
    
    # Landau maximum location: in the Landau distribution
    # (represented by the CERNLIB approximation), the maximum
    # is located at -0.22278298 with the location parameter = 0
    # This shift is corrected within this function, so that the actual 
    # maximum is identical to the MP parameter
    mpshift  = -0.22278298
    
    # MP shift correction
    mpc = par[0]-mpshift*par[1]
    # -- Fit parameters:
    #   par[0] = Most Probable (MP, location) parameter of Landau density
    #   par[1] = Width (scale) parameter of Landau density
    #   par[2] = Total area (integral -inf to inf, normalization constant)
    #   par[3] = Gaussian smearing
    
    # -- Control constants
    #  Number of convolution steps
    np = 100.0
    #  Convolution extends to +- sc gaussian sigmas
    sc = 5.0
    # -- Range of convolution integral
    xlow = max(0.0,x[0]-sc*par[3])
    xupp = x[0]+sc*par[3]
    step =(xupp-xlow)/np
    
    # -- Convolution integral of Landau and Gaussian by sum
    xlow_i = lambda i: xlow+(i-0.5)*step
    xup_i  = lambda i: xupp-(i-0.5)*step
    
    total = 0.0
    for i in xrange(1,int(np/2.0)+1):
        xi = xlow_i(i)
        xxi = xup_i(i)
        total += sum(map(lambda ix: ROOT.TMath.Landau(ix,mpc,par[1])/par[1]*\
                ROOT.TMath.Gaus(x[0],ix,par[3]),[xlow_i(i),xup_i(i)]))
    
    return par[2]*1./sqrt(pi)*step*total/par[3]


def fit_langaus(h,xmin=0,xmax=60,peak_min=0,peak_max=60,force_range=False):
    """Fit a Landau * Gauss to the input histogram, returning
    the main parameters of the fitted function, as well as if
    the fit was succesful. The function gets attached to the
    histogram and can be retrieved with the h.GetFunction method

    Parameters
    ----------
    h: ROOT.TH1
        The histogram to fit
    xmin: float, optional
        The range minimum to fit
    xmax: float, optional
        The range maximum to fit
    peak_min: float, optional
        The minimum allowed value where the peak should be
    peak_max: float, optional
        The maximum allowed value where the peak should be
    force_range: bool, optional
        Whether or not use tunned range for the fit, 
        If yes, it will apply a range defined around
        the found peak and using an estimated sigma (peak/7.0):
          [peak - 2.5 *Sigma , peak + 6.0 *Sigma ]

    Return
    ------
    (MPV,Sigma gaus,TResultsStatus,TF1 fitted)
    """
    import ROOT

    # Define the TF1 
    name = h.GetName()+"_landaugaus_tf1"
    lg = ROOT.TF1(name,landau_gaus,xmin,xmax,4)
    # Start values for the parameters. Some guese
    # -- the Most probable value
    #peak = h.GetBinCenter(h.GetMaximumBin())
    peak = spectrum_searcher_landau_parameters(h,peak_min,peak_max)
    # -- Sigma
    sm = peak/7.0
    if force_range:
        xmin = peak-2.5*sm
        xmax = peak+6.0*sm
        if xmin < 0.5:
            xmin = 0.5
    # -- Gaussian noise
    lg.SetParameter(0,peak)
    #lg.SetParLimits(0,1,50)
    # -- Scale of the Landau
    lg.SetParameter(1,sm)
    # Looks like things below 0.1-0.15 
    # Only make sense if ToT
    if peak < 50:
        lg.SetParLimits(1,0.10,10)
    # -- Normalization constant
    norm=h.Integral()
    lg.SetParameter(2,norm)
    #lg.SetParLimits(2,norm*0.001,norm*3.0)
    # -- Gauss width (initially as the Landau width)
    lg.SetParameter(3,sm)
    #lg.SetParLimits(3,1.e-3,10)

    lg.SetLineColor(46)
    lg.SetLineWidth(2)
    results = h.Fit(lg,"RSQ+","",xmin,xmax)
    return lg.GetParameter(0),lg.GetParameter(3),results.Status(),lg

def get_time_window(tree,cut,\
        win_halfsize=2.0,\
        fog_brname="cluster_charge",time_brname="eventTime",
        fog_min=0,fog_max=500):
    """Get the time window by obtaining the maximum of the 
    time profile 
    
    Parameters
    ----------
    tree: ROOT.TTree
    cut: str
        The cut string to be applied to the selection
    win_halfsize: float, default: 2.0 
        The (half) size of the window around the maximum
    fog_brname: str, optional
        The name of the branch will act as figure of merit
        to extract their maximum (values are averaged vs. x)
    time_brname: str, optional
        The name of the branch will act as x-coordinate 
    fog_min: float
        The min value for the fog in the fit
    fog_max: float
        The max value for the fog in the fit
    """
    import ROOT
    ROOT.gROOT.SetBatch()

    c = ROOT.TCanvas()
    # Check the branches exist
    for brname in [ fog_brname, time_brname ]:
        if brname not in map(lambda x: x.GetName(),tree.GetListOfBranches()):
            raise AttributeError("get_time_window: no branch '{0}' "\
                    "found in the input tree".format(brname))
    # First get the eta-distribution
    tree.Draw("{0}:{1}>>_test(30,-0.5,29.5,100,{2},{3})".format(\
            fog_brname,time_brname,fog_min,fog_max),"{0}".format(cut))
    _h2 = ROOT.gDirectory.Get("_test")
    if type(_h2) == ROOT.TObject:
        raise AttributeError("Histogram generation failed. Probably"\
                " due to some unexisting branch name introduced in the"\
                " cut definition: '{0}' or in the Figure of merit branch: '{1}'"\
                " or in the Time branch: '{2} or not enough Events' ".format(cut,fog_brname,time_brname))
    _h2.SetDirectory(0)
    prof_window = _h2.ProfileX()
    # Fit a polynomial
    t0=3.0
    t1=30.0
    pol_func = ROOT.TF1("pol_func","pol4",t0,t1)
    status = prof_window.Fit(pol_func,"QS")
    if status.Status() != 0:
        print "[analysis_functions.get_time_window] WARNING : Fit failed"
        ROOT.gROOT.SetBatch(0)
        c.Close()
        return (-1,-1)
    # Find the maximum and define time window as +- 2 ns.
    # To find the maximum df/dt ~= 0
    N = 200
    Dx = (t1-t0)/float(N)
    tmax=0.0
    found_max=False
    last_derivate = pol_func.Derivative(t0)
    for i in xrange(1,N):
        t = t0+i*Dx
        # Change of sign?
        current_derivate =  pol_func.Derivative(t)
        if(current_derivate*last_derivate < 0.0):
            # Get the medium 
            tmax = t-Dx/2.0
            found_max = True
            break
        last_derivate = current_derivate
    # If did not find the maximum, means is monotically
    # creasing or decreasing, take the absolute maximum
    if not found_max:
        # Creasing case
        if pol_func(t0) < pol_func(t1):
            tmax = t1-win_halfsize
        # Decreasing case
        else:
            tmax = t0
    ROOT.gROOT.SetBatch(0)
    c.Close()
    return (tmax-win_halfsize,tmax+win_halfsize)

def get_eta_distribution(tree,eta_brname,cut):
    """Obtain the eta distribution from a tree

    Parameters
    ----------
    tree: ROOT.TTree
    eta_brname: str
        The branch of the eta variable
    cut: str
        The cut condition

    Returns
    -------
    eta: ROOT.TH1F 
        the eta distribution with the applied cut 
        conditions
    """
    import ROOT
    ROOT.gROOT.SetBatch()
    
    # Check the branches exist (not trivial to extract branches
    # names from the cut
    for brname in [eta_brname]:
        if brname not in map(lambda x: x.GetName(),tree.GetListOfBranches()):
            raise AttributeError("get_time_window: no branch '{0}' "\
                    "found in the input tree".format(brname))
    # First get the eta-distribution
    c = ROOT.TCanvas()
    tree.Draw("{0}>>eta_2(100,-0.1,1.1)".format(eta_brname),"{0}".format(cut))
    eta = ROOT.gDirectory.Get("eta_2")
    if type(eta) == ROOT.TObject:
        raise AttributeError("Histogram generation failed. Probably"\
                " due to some unexisting branch name introduced in the"\
                " cut definition: '{0}'".format(cut))
    eta.SetDirectory(0)
    eta.SetTitle("#eta-distribution;#eta;Clusters")

    c.Close()
    ROOT.gROOT.SetBatch(0)
    return eta

class memoize_eta(object):
    """Memoize the eta distribution and converting it
    to the cumulative
    XXX: DOC
    """
    def __init__(self,f):
        """
        Parameters
        ----------
        f: functor
        """
        self.histo_eta_c = None
        self.norm  = None
        self.eta_c = None
        self.histo_eta_density = None
    def __call__(self,eta,eta_d=None):
        """Obtain the wrapper which actually 
        returns the f(eta)=Int(0,eta)/Norm
        XXX DOC
        """
        import ROOT
        from math import sqrt

        if not self.histo_eta_c and not self.norm:
            self.histo_eta_c = eta_d.GetCumulative()
            self.histo_eta_c.SetDirectory(0)
            # Get the interpolation
            self.eta_c = self.histo_eta_c.Interpolate
            # See also: Beta distribution and Kumaraswamy distribution
            # the logit-normal: https://en.wikipedia.org/wiki/Logit-normal_distribution
            self.norm  = 1./float(eta_d.Integral(1,eta_d.GetNbinsX()))
            self.histo_eta_density = eta_d
            self.histo_eta_density.SetDirectory(0)
            # Update as well the histograms
        return self.eta_c(eta)*self.norm

@memoize_eta
def frac_position(eta,eta_distribution):
    """Get the fractionary position gived an eta.
    See the func wrapper.

    Parameters
    ----------
    eta: float
        The eta value
    eta_distribution: ROOT.TH1F
        The eta distribution 
    """
    pass
