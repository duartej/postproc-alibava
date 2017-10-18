#!/usr/bin/env python
"""Post-processor module to be applied to the output of the EUTelescope 
Marlin jobs. This module acts as liason between the Marlin EUTelTreeCreator
processor (the creator of the ntuple) and the post-processing python
"""
__author__ = "Jordi Duarte-Campderros"
__credits__ = ["Jordi Duarte-Campderros"]
__version__ = "0.1"
__maintainer__ = "Jordi Duarte-Campderros"
__email__ = "jorge.duarte.campderros@cern.ch"
__status__ = "Development"
# XXX:
# 2. improve doc
# 3. Get all the metadata from SPS2017TB class
# 4. Be sure plots are in the LOCAL coordinates of the plane as well, i.e. channel (see the formula to convert in biblio)
# 5. Be sure the binning is taking into account the strip resolution (at least bin_width = 1/4*pitch or so)

# Resolution of the DUT plus the error in the fit (?) 0.04 mm 
# Incluir covariance en los hits fitted
# TODO OBTAIN FROM SPS
REFPLANE=7
DUTPLANE=6
RESOLUTION= { DUTPLANE: 0.2, REFPLANE : 0.1 } # in mm
SENSOR_NAME = { REFPLANE: "REF", DUTPLANE: "DUT" }

# UNITS
MM = 1.0
UM = 1e3

class metadata_container(object):
    """Container of some globals which are going to be
    used along the module; these globals are run/file
    dependent.
    """
    def __init__(self,tree,dut_name):
        """Use the tree structure to define the DUT 
        planes ids, their resolutions, etc

        Parameters
        ----------
        tree: ROOT.Tree
            The ntuple containing tracks and DUT, REF sensors
            hits
        dut_name: str
            The 

        Raises
        ------
        RuntimeError
            When the tree doesn't contain the Branch hit_X_#
        """
        from .SPS2017TB_metadata import sensor_name_spec_map
        import math

        # Names (following conventions at .SPS2017TB_metadata
        self.ref_name = "REF_0_b1"
        self.dut_name = dut_name

        # The REF plane is explicitely set to 7 ALWAYS
        self.ref_plane = 7
        # Get the DUT by using the hit branches
        try:
            self.dut_plane = map(lambda br: int(br.GetName().replace("hit_X_","")),
                filter(lambda x: x.GetName().find("hit_X_") == 0, tree.GetListOfBranches()))[0]
        except IndexError:
            raise RuntimeError("Unexpected tree structure")
        # Resolutions: REF and DUT (binary so far)
        _p = 1.0/math.sqrt(12.0)
        self.resolution = { self.ref_plane: sensor_name_spec_map[self.ref_name].pitchX*_p/MM , 
                self.dut_plane: sensor_name_spec_map[self.dut_name].pitchX*_p/MM }

class dummy_list(list):
    """A dummy list which returns always
    the value in the first element
    """
    def __getitem__(self,i):
        return list.__getitem__(self,0)

class hits_plane_accessor(object):
    """Access to the TBranches contained in a 
    TTree created by the EUTelTreeCreator class from
    the https://github.com/duartej/eutelescope/ package.
    
    Attributes
    ----------
    id: int
        The telescope plane ID
    x: ROOT.std.vector(float)()
        The vector of x-position of the hits 
    y: ROOT.std.vector(float)()
        The vector of y-position of the hits
    z: dummy_list instance
        The vector of z-position of the hits, which is 
        the same for all the hits belonging to the same 
        plane. Therefore, it is extracted in the first
        event from the branch, and afterwards used in
        all the events
    x_local: ROOT.std.vector(float)()
        The vector of x local position of the hits 
    y_local: ROOT.std.vector(float)()
        The vector of y local position of the hits
    charge: ROOT.std.vector(int)()
        The hit total charge in ADCs counts 
    """
    def __init__(self,tree,planeid):
        """Accessors for the branches related with Alibava
        sensors hits

        Parameters
        ----------
        tree: ROOT.TTree
            The ntuple 
        planeid: int
            The sensor plane ID
        
        """
        self.id = planeid
        # Check if there is any problem
        if len(filter(lambda br: br.GetName() == "hit_X_{0}".format(self.id) ,tree.GetListOfBranches())) == 0:
            raise RuntimeError("The INDEX '{0}' does not identify any Alibava"\
                    " sensor plane".format(self.id))
        # Activate the tree if wasn't
        dummy = tree.GetEntry(0)
        # Fill the elements
        self.x = getattr(tree,"hit_X_{0}".format(self.id))
        self.y = getattr(tree,"hit_Y_{0}".format(self.id))
        self.x_local = getattr(tree,"hit_XLocal_{0}".format(self.id))
        self.y_local = getattr(tree,"hit_YLocal_{0}".format(self.id))
        self.charge = getattr(tree,"hit_total_charge_{0}".format(self.id))
        
        # Let's assume that z doesn't change between events, as 
        # it is a fixed constrain, so find the first entry with
        # a valid value and then obtain
        branch_z = getattr(tree,"hit_Z_{0}".format(self.id))
        self.z = dummy_list()
        for i in xrange(tree.GetEntries()):
            dum = tree.GetEntry(i)
            if branch_z.size() > 0:
                self.z.append(branch_z[0])
                break
        if len(self.z) == 0:
            raise RuntimeError("No event was found with a valid value "\
                    "at Z for the {0}-plane. Something in the input file"\
                    " is corrupted")
    
    @property
    def n(self):
        """The number of hit elements

        Return
        ------
        int
        """
        return self.x.size()
    
    def equal(self,i,other,j):
        """Allow comparation beetween hits, we assume equal
        if (x,z) are equal

        Parameters
        ----------
        i: int
            The hit index of this (self) instance to be 
            compare
        other: hit_plane_accesor instance (or any other instance
            with a method 'z.__getitem__')
            The other instance to compare with
        j: int
            The hit index of the 'other' instance

        Return
        ------
        bool
        """
        return abs(self.x[i]-other.x[j]) < 1e-9 \
                and abs(self.z[i]-other.z[j]) < 1e-9


    def matched_hits(self,track_acc,histo):
        """Search for tracks matching (within a resolution) the list 
        of hits in the current event

        Parameters
        ----------
        track_acc: tracks_accessor instance
            The accessor to the telescope tracks branches
        histo: ROOT.TH2
            The histogram to store correlation between the 
            hit in the sensor and the hit in the track

        Return
        ------
        hitlist: dict(int,list(int))
            The list of indices of the hits matched with a track 
        """
        matched_hits = []
        # Given the hits in that sensor,
        for ihit in xrange(self.n):
            closest = {}
            # Evaluate the distance from the hit to any tracks
            for itrk in xrange(track_acc.n):
                xpred,ypred=track_acc.get_point(itrk,self.z[ihit])
                closest[abs(xpred-self.x[ihit])] = (itrk,xpred,self.x[ihit])
            if len(closest) == 0:
                continue
            # if any track matched, get the closest one 
            distance,(trk_el,xpredicted,x_in_hit) =sorted(closest.iteritems())[0]

            # Fill correlation histograms
            histo.Fill(x_in_hit,xpredicted)
            # define resolution (an input maybe?)
            # and olny add the track matched within the 
            # resolution
            res = RESOLUTION[self.id]
            if res > distance:
                matched_hits.append( (ihit,trk_el) )

        return matched_hits


class track_hits_plane_accessor(object):
    """Access to the TBranches contained in a 
    TTree created by the EUTelTreeCreator class from
    the https://github.com/duartej/eutelescope/ package.

    Attributes
    ----------
    id: int
        The telescope plane ID
    x: ROOT.std.vector(float)()
        The vector of x-position of the hits 
    y: ROOT.std.vector(float)()
        The vector of y-position of the hits
    z: ROOT.std.vector(float)()
        The vector of z-position of the hits. In the case
        of measured hits, the plane is obviously restricting
        the z to be the same for all hits
    x_local: ROOT.std.vector(float)()
        The vector of x local position of the hits 
    y_local: ROOT.std.vector(float)()
        The vector of y local position of the hits
    track_index: ROOT.std.vector(int)()
        The index (a number identifying univocaly a 
        track in the event) of the track which the 
        hit is associated
    """
    def __init__(self,tree,planeid,ismeas,ref_plane_id,dut_plane_id):
        """Accessors for the branches related with the
        telescope hits. The hits are associated to 
        reconstructed tracks of the telescope, and therefore
        can refer to the measured or the fitted one.
        

        Parameters
        ----------
        tree: ROOT.TTree
            The ntuple 
        planeid: int
            The telescope plane ID
        ismeas: bool
            Whether this instance should deal with associated
            measured hits (true) or fitted hits (false)
        ref_plane_id: int
            The plane ID of the Alibava-REF sensor
        dut_plane_id: int
            The plane ID of the Alibava-DUT sensor

        Raises
        ------
        AssertionError
            When the REF and DUT plane ID introduced are exactly
            the same
        RuntimeError
            The plane index introduced does not belong to any
            (checking the 'trk_hit_meas(|fit)_X_<planeID>' brach)
        """
        if ismeas:
            prefix = "trk_hit_meas"
        else:
            prefix = "trk_hit_fit"

        self.id = planeid
        self.is_dut = (self.id == dut_plane_id)
        self.is_ref = (self.id == ref_plane_id)
        # A consistency check
        if self.is_ref and self.is_dut:
            raise AssertionError("The plane index for the DUT and the "\
                    "REF Alibava sensor planes are exactly the same: {0}".format(dut_plane_id))
        # Check if there is any problem
        if len(filter(lambda br: br.GetName() == "{0}_X_{1}".format(prefix,self.id) ,tree.GetListOfBranches())) == 0:
            raise RuntimeError("The INDEX '{0}' does not identify any Telescope"\
                    " sensor plane".format(self.id))
        # Activate the tree if wasn't
        dummy = tree.GetEntry(0)
        # Fill the elements
        self.x = getattr(tree,"{0}_X_{1}".format(prefix,self.id))
        self.y = getattr(tree,"{0}_Y_{1}".format(prefix,self.id))
        self.x_local = getattr(tree,"{0}_XLocal_{1}".format(prefix,self.id))
        self.y_local = getattr(tree,"{0}_YLocal_{1}".format(prefix,self.id))
        self.z     = getattr(tree,"{0}_Z_{1}".format(prefix,self.id))
        self.track_index = getattr(tree,"{0}_index_{1}".format(prefix,self.id))
    
    @property
    def n(self):
        """The number of hit elements

        Return
        ------
        int
        """
        return self.x.size()
    
    def equal(self,i,other,j):
        """Allow comparation beetween hits, we assume equal
        if (x,z) are equal

        Parameters
        ----------
        i: int
            The hit index of this (self) instance to be 
            compare
        other: track_hits_plane_accesor instance (or any other instance
            with a method 'z.__getitem__')
            The other instance to compare with
        j: int
            The hit index of the 'other' instance

        Return
        ------
        bool
        """
        return abs(self.x[i]-other.x[j]) < 1e-9 \
                and abs(self.z[i]-other.z[j]) < 1e-9

    def get_hits_with_same_track_indices(self,itrk):
        """Return the list of indices for those hits belonging 
        to the itrk-track

        Parameters
        ----------
        itrk: int
            The index of the track

        Return
        ------
        list(int): The list of hits (indices) of the itrk-track
        """
        return filter(lambda i: self.track_index[i] == itrk,xrange(self.n))

    def get_hit_with_same_track_index(self,i,other):
        """A track has a list of associated hits (either measured
        or fitted). This function returns the index of the hit
        of a different plane which belongs to the same track than
        the i-hit of this instance.

        Parameters
        ----------
        i: int
            The index of the hit to check
        other: track_hits_plane_accessor
            The other plane to check

        Return
        ------
        int: The corresponding hit-index of the other instance associated
            to the same track

        Raises
        ------
        RuntimeError
            Whenever there is more than one element matched in
            the other plane. 
        """
        matched = filter(lambda other_i: other.track_index[other_i] == self.track_index[i],\
                xrange(other.n))
        if len(matched) == 0:
            return None
        elif len(matched) == 1:
            return matched[0]
        else:
            raise RuntimeError("Unexpected behaviour: more than elements match")


    def is_dut(self):
        """Whether the plane corresponds to the Alibava DUT
        """
        return self.is_dut
    
    def is_ref(self):
        """Whether the plane corresponds to the Alibava REF
        """
        return self.is_ref
    

class tracks_accessor(object):
    """Access to the TBranches contained in a 
    TTree created by the EUTelTreeCreator class from
    the https://github.com/duartej/eutelescope/ package.
    
    A track is defined as a straight line using 
    a reference point (r0) and a vector director (v).
    Any point of the track can be obtained with
     x(t)=r0+v(t)
    using the parameter t
    Each track has a list of associated hits as well.

    
    Attributes
    ----------
    refid: int
        The Alibava REF ID
    dutid: int
        The Alibava DUT ID
    tel_planes: list(int)
        The list of telescope planes
    matched_hits: functor
        See _matched_hits method
    x0: float
        The x-position of the reference point
    y0: float
        The y-position of the reference point
    z0: float
        The z-position of the reference point
    z0: float
        The z-position of the reference point
    dxdz: float
        The vector director (x component)
    dydz: float
        The vector director (y component)
    telescope_measured_hits: list(track_hits_plane_accesor)
        The telescope (measured) hits organized by planes
    telescope_fitted_hits: list(track_hits_plane_accesor)
        The telescope (fitted) hits organized by planes
    dut_measured_hits: track_hits_plane_accesor
        The DUT (measured hits), note that this data member
        is only usable if the tracks were fitted in the 
        Marlin processor using the DUT sensor as telescope
        plane
    dut_fitted_hits: track_hits_plane_accesor
        The DUT (fitted hits), note that this data member
        is only usable if the tracks were fitted in the 
        Marlin processor using the DUT sensor as telescope
        plane
    ref_measured_hits: track_hits_plane_accesor
        The REF (measured hits), note that this data member
        is only usable if the tracks were fitted in the 
        Marlin processor using the DUT sensor as telescope
        plane
    ref_fitted_hits: track_hits_plane_accesor
        The REF (fitted hits), note that this data member
        is only usable if the tracks were fitted in the 
        Marlin processor using the DUT sensor as telescope
        plane
    """
    def __init__(self,tree,tel_planes,ref_plane,dut_plane):
        """Accessors for the branches related with the
        telescope tracks.       
        
        Parameters
        ----------
        tree: ROOT.TTree
            The ntuple 
        tel_planes: list(int)
            The list of telescope planes ID
        ref_plane_id: int
            The plane ID of the Alibava-REF sensor
        dut_plane_id: int
            The plane ID of the Alibava-DUT sensor

        Raises
        ------
        AssertionError
            When the REF and DUT plane ID introduced are exactly
            the same
        RuntimeError
            The plane index introduced does not belong to any
            (checking the 'trk_hit_meas(|fit)_X_<planeID>' brach)
        """
        self.refid = ref_plane
        self.dutid = dut_plane
        self.tel_planes = tel_planes
        # Check if the tracks were fitted with the DUT and REF 
        # planes, in order to activate or not the matched_hits 
        # function
        dutref_present=True
        for _i in [self.dutid+self.refid]:
            dutref_present = dutref_present and (len(filter(lambda br: br.GetName() == "trk_hit_fit_X_{0}".format(_i),\
                tree.GetListOfBranches())) == 0)
        if dutref_present:
            self.matched_hits = self._matched_hits
        else:
            print "\033[1;33mWARNING\033[1;m No DUT or/and REF hits were used "\
                    "to fit the tracks. No sense to use 'matched_hits' function,"\
                    " use the one from the 'hits_plane_accessor' class"
            self.matched_hits = lambda _d,_d1,_d3: {}

        # Track specific
        # Reference point
        self.x0 = tree.trk_refPoint_X
        self.y0 = tree.trk_refPoint_Y
        self.z0 = tree.trk_refPoint_Z
        # Director vector
        self.dxdz = tree.trk_dxdz
        self.dydz = tree.trk_dydz
        # -- the associated hits
        # --- telescope
        self.telescope_measured_hits = map(lambda i: track_hits_plane_accessor(tree,i,True,self.refid,self.dutid),self.tel_planes)
        self.telescope_fitted_hits = map(lambda i: track_hits_plane_accessor(tree,i,False,self.refid,self.dutid),self.tel_planes)
        # --- dut
        self.dut_measured_hits = track_hits_plane_accessor(tree,self.dutid,True,self.refid,self.dutid)
        self.dut_fitted_hits = track_hits_plane_accessor(tree,self.dutid,False,self.refid,self.dutid)
        # --- ref
        self.ref_measured_hits = track_hits_plane_accessor(tree,self.refid,True,self.refid,self.dutid)
        self.ref_fitted_hits = track_hits_plane_accessor(tree,self.refid,False,self.refid,self.dutid)
        
    @property
    def n(self):
        """The number of tracks

        Return
        ------
        int
        """
        return self.x0.size()

    def _matched_hits(self,duthits,refhits,h):
        """Function only usable if the tracks were fitted including
        both the DUT and REF sensor hits. 
        Extract those measured hits in the DUT and REF which were 
        used to create a track, meaning that there is a match with
        the track and those hits. 

        Note
        ----
        In order to have an unbiased set of matched hits, the input
        tracks must have been created using a special treatment for
        the DUT/REF hits in the Y-resolution, this is not the case
        so far. Therefore: DO NOT USE THIS FUNCTION until opposite
        established
        
        Parameters
        ----------
        duthits: hit_planes_accessor
            The accessor to the DUT hits
        refhits: hit_planes_accessor
            The accessor to the REF hits

        Return
        ------
        hitlist: dict(int,list(int))
            The list of indices of the hits used to fit a track per
            detector (identified by its plane ID). Note that the indices
            are referring to the hit_planes_accessor collections used
            as the input arguments, i.e. duthits and refhits.
        """
        # Just defined a helper dict to be able to perform the same algorithm 
        # to both sensors:
        # -- The hit_planes accessor instances
        original_measured_hits = { self.dutid: duthits, self.refid:refhits}
        # -- And the measured and fitted accessor, data members of these tracks
        automatic_list = [ (self.dutid,self.dut_measured_hits,self.dut_fitted_hits),\
                (self.refid,self.ref_measured_hits,self.ref_fitted_hits)]

        # The output dict: will contain the element index of those DUT/REF 
        # hits which were used to fit a track (note that the index is referring
        # to the hit_planes_accessor instances, introduced as arguments)
        hitlist ={ self.dutid: [], self.refid: []}
        for itrk in xrange(self.n):
            for sID,measured_hits,fitted_hits in automatic_list:
                # -- Get those hits in the DUT and REFthat were actually used
                #    in the fit
                for i in measured_hits.get_hits_with_same_track_indices(itrk):
                    ihit_fit = measured_hits.get_hit_with_same_track_index(i,fitted_hits)
                    if not ihit_fit:
                        continue
                    # The itrk-track was using in the fit the measured i-hit: it's a match
                    # Get the index of the hit_plane_accessor instance
                    i_in_hit_plane = filter(lambda k: measured_hits.equal(i,original_measured_hits[sID],k),\
                            xrange(original_measured_hits[sID].n))[0]
                    #print self.get_point(itrk,original_measured_hits[sID].z[i_in_hit_plane]),fitted_hits.x[ihit_fit],fitted_hits.y[ihit_fit]
                    h[sID].Fill(original_measured_hits[sID].x[i_in_hit_plane],self.get_point(itrk,original_measured_hits[sID].z[i_in_hit_plane])[0])
                    hitlist[sID].append( (i_in_hit_plane,itrk) )
        return hitlist
    
    def get_point(self,i,z):
        """Give the predicted position (x,y) for the i-track

        Parameters
        ----------
        i: int
            The index of the track
        z: float
            The z-position where to predict 
        
        Return
        ------
        (float,float): Predicted position at the z-plane
        """
        t = z-self.z0[i]
        return (self.x0[i]+self.dxdz[i]*t,self.y0[i]+self.dydz[i]*t)

class processor(object):
    """Results extractor. Histograms are defined and filled here
    XXX MISSING DOC
    """
    def __init__(self,minst):
        """ XXX MISSING DOC
        Parameters
        ----------
        minst: metadata_container instance
        """
        import ROOT

        # some info numbers
        self.total_events    = 0
        self.events_with_dut = 0
        self.dut_hits        = 0
        self.events_with_ref = 0
        self.ref_hits        = 0
        self.events_with_dut_no_tracks = 0
        self.events_with_ref_no_tracks = 0
        self.events_with_ref_dut_tracks = 0
        # Some statistic histograms
        self.residual_projection = { minst.dut_plane: ROOT.TH1F("res_projection_dut"," ; x_{DUT}-x_{trk}^{pred} [mm]; Entries",200,-0.5*MM,0.5*MM),
                minst.ref_plane: ROOT.TH1F("res_projection_ref"," ; _x_{REF}-x_{trk}^{pred} [mm]; Entries",200,-0.5*MM,0.5*MM) }
        # -- histograms DUTS/REF
        self.hmap = { minst.dut_plane: ROOT.TH2F("local_map_dut" ,";x [mm]; y [mm]; Entries", 130,-10.0*MM,10.0*MM,100,-10.0*MM,10.0*MM),
                minst.ref_plane: ROOT.TH2F("local_map_ref",";x [mm]; y [mm]; Entries", 130,-10.0*MM,10.0*MM,100,-10.0*MM,10.0*MM) }
        self.hcharge = { minst.dut_plane: ROOT.TH2F("precharge_map_dut",";x [mm]; y [mm]; charge cluster [ADC]",250,-10.0*MM,10.0*MM,250,-10.0*MM,10.0*MM),
                minst.ref_plane: ROOT.TH2F("precharge_map_ref",";x [mm]; y [mm]; charge cluster [ADC]",250,-10.0*MM,10.0*MM,250,-10.0*MM,10.0*MM) }
        # Keep the number of entries used in each bin
        self.haux_charge = { minst.dut_plane: ROOT.TH2F("aux_charge_dut",";x [mm]; y [mm]; Entries",250,-10.0*MM,10.0*MM,250,-10.0*MM,10.0*MM),
                minst.ref_plane: ROOT.TH2F("aux_charge_ref",";x [mm]; y [mm]; Entries",250,-10.0*MM,10.0*MM,250,-10.0*MM,10.0*MM) }
        # Correlations
        self.hcorr_trkX = { minst.dut_plane: ROOT.TH2F("corr_trkX_dut",";x_{DUT} [mm]; x_{pred}^{trk} [mm]; Entries",200,-10.0,10.0,200,-10.0,10.0),
                minst.ref_plane: ROOT.TH2F("corr_trkX_ref",";x_{REF} [mm]; x_{pred}^{trk} [mm]; Entries",200,-10.0,10.0,200,-10.0,10.0)}
        self.hcorrX = ROOT.TH2F("corrX_dut_ref",";x_{DUT} [mm]; x_{REF} [mm]; Entries",100,-10.0,10.0,100,-10.0,10.0)
        self.hcorrY = ROOT.TH2F("corrY_dut_ref",";y_{DUT} [mm]; y_{REF} [mm]; Entries",100,-10.0,10.0,100,-10.0,10.0)
        # efficiendy
        self.hpass = ROOT.TH2F("pass_map",";x [mm]; y [mm]; #varepsilon",130,-10.0,10.0,100,-10.0,10.0)
        self.htotal = ROOT.TH2F("total_map",";x [mm]; y [mm]; #varepsilon",130,-10.0,10.0,100,-10.0,10.0)

        self._allhistograms = self.residual_projection.values()+self.hmap.values()+\
                self.hcharge.values()+self.haux_charge.values()+self.hcorr_trkX.values()+\
                [self.hcorrX,self.hcorrY,self.hpass,self.htotal]
        dummy=map(lambda h: h.SetDirectory(0),self._allhistograms)

    #def draw_all_in_canvas():
    #    """
    #    """
    #    import ROOT
    def create_efficiency(self):
        """Create the hit efficiency map using the total and pass
        histograms
        """
        import ROOT
        self.eff = ROOT.TEfficiency(self.hpass,self.htotal)
        self.eff.SetName("sensor_efficiency")
        self.eff.SetDirectory(0)
        self._allhistograms.append(self.eff)

    def create_charge_map(self):
        """Create the charge map using the two auxiliary histograms,
        one containing the entries for a x,y position and the other
        the total charge of all the entries
        """
        import ROOT
        self.hcharge_average = dict([(s,h.Clone(h.GetName().replace("pre",""))) \
                for s,h in self.hcharge.iteritems()])
        for s,h in self.hcharge_average.iteritems():
            h.GetZaxis().SetTitle("<charge cluster> [ADC]")
            h.Divide(self.haux_charge[s])
            self._allhistograms.append(h)

    def store_histos(self,filename):
        """
        """
        import ROOT
        # Create the efficiency before recording them
        self.create_efficiency()
        # Create the map of average charges
        self.create_charge_map()
        # Store all the histograms
        f=ROOT.TFile.Open(filename,"RECREATE")
        for h in self._allhistograms:
            h.Write()
        f.Close()

    def fill_statistics(self,ndut,nref,tracks):
        """Fill some statistic related with the number of dut, 
        ref hits and tracks found.

        Parameters
        ----------
        ndut: int
            Number of DUT hits
        nref: int
            Number of REF hits
        tracks: int
            Number of reconstructed tracks
        """
        self.total_events += 1
        if ndut != 0:
            self.events_with_dut += 1
            self.dut_hits += ndut
            if tracks == 0:
                self.events_with_dut_no_tracks += 1 
        if nref != 0:
            self.events_with_ref += 1
            self.ref_hits += nref
            if tracks == 0:
                self.events_with_ref_no_tracks += 1 
        if ndut != 0 and nref != 0 and tracks != 0:
            self.events_with_ref_dut_tracks += 1
        

    def process(self,t,trks,refhits,duthits):
        """Fill the relevant counters
        """
        import math
        
        # Get the number of hits in the DUT and REF plane
        self.fill_statistics(duthits.n,refhits.n,t.ntracks)

        # Fill histograms for the DUT/REF
        # -- Charge map: using a matched track to the sensor, in order to
        # find the missing Y
        # -- Efficiency map: matched track to the reference, and extrapolate to
        # the DUT plane (total histogram), match to the DUT  (pass histogram)
        
        # - a dict to get the sensors by its plane ID
        sensor_hits = { refhits.id: refhits, duthits.id: duthits }

        # Get the matched hits { sensorID: [ (hit_index,track_index,. ..] , }
        #matched_hits = trks.matched_hits(duthits,refhits,self.hcorr_trkX)
        matched_hits = { duthits.id: duthits.matched_hits(trks,self.hcorr_trkX[duthits.id]),
                refhits.id: refhits.matched_hits(trks,self.hcorr_trkX[refhits.id]) }
        predictions = {}
        # Filling histograms
        for sensorID,pointlist in matched_hits.iteritems():
            hits = sensor_hits[sensorID]
            for hit_el,trk_index in pointlist:
                xpred,ypred = trks.get_point(trk_index,hits.z[hit_el])
                # Residuals
                self.residual_projection[sensorID].Fill(xpred-hits.x[hit_el])
                # Event and charge maps (using track predictions)
                self.hmap[sensorID].Fill(xpred,ypred)
                self.hcharge[sensorID].Fill(xpred,ypred,hits.charge[hit_el])
                self.haux_charge[sensorID].Fill(xpred,ypred)
        # Matching beetwen DUT and REF: hits with the same track ID
        for ihit_ref,itrk in matched_hits[refhits.id]:
            # Fill the total histogram (using the prediction in DUT)
            x,y = trks.get_point(itrk,duthits.z[0])
            self.htotal.Fill(x,y)
            dut_match = filter(lambda (ihit_dut,itrk_dut): itrk_dut == itrk,matched_hits[duthits.id])
            if len(dut_match) > 0:
                self.hpass.Fill(x,y)
                # And some correlation plots of the predicted values
                x_at_ref,y_at_ref = trks.get_point(itrk,refhits.z[ihit_ref])
                self.hcorrX.Fill(x,x_at_ref)
                self.hcorrY.Fill(y,y_at_ref)


    def __str__(self):
        """Summarize the information
        """
        m = "|-------------------------------------------|\n"
        m+= " Events with DUT: {0} (eff: {1:.2f}%)\n".format(self.events_with_dut,float(self.events_with_dut)/float(self.total_events)*100.)
        m+= " Events with REF: {0} (eff: {1:.2f}%)\n".format(self.events_with_ref,float(self.events_with_ref)/float(self.total_events)*100.)
        m+= "\n Events with DUT but no tracks: {0}\n".format(self.events_with_dut_no_tracks)
        m+= " Events with REF but no tracks: {0}\n".format(self.events_with_ref_no_tracks)
        m+= " Events with REF and DUT and tracks: {0}\n".format(self.events_with_ref_dut_tracks)
        m+= "\n Total processed events: {0}\n".format(self.total_events)
        m+= "|-------------------------------------------|"

        return m

def sensor_map_production(fname):
    """Produce the charge maps and the efficiency maps (hit 
    reconstruction efficiency) for the sensor included in the 
    NTUPLE.

    Parameters
    ----------
    fname: str
        The name of the ROOT ntuple produced by the EUTelTreeCreator
        processor (part of the EUTelescope package, which uses the 
        Marlin framework). The filename MUST follow the standard 
        notation defined through the class 
        `alibavaSkifftools.SPS2017TB_metadata.filename_parser`

    Raises
    ------
    IOError
        Whenever the ROOT file is not present
    """
    import ROOT
    import timeit
    import sys
    import os
    from .SPS2017TB_metadata import filename_parser 
    from .SPS2017TB_metadata import standard_sensor_name_map as name_converter

    # Get the name of the sensor and some other useful info
    fp = filename_parser(fname)
    # Get the root file and the tree
    f = ROOT.TFile.Open(fname)
    if not f or f.IsZombie():
        raise IOError("Invalid or not found ROOT file '{0}'".format(fname))
    try:
        t = f.Get("events")
    except ReferenceError:
        raise ReferenceError("Invalid NTUPLE format at '{0}'".format(fname))
    # And obtain some metadata (ID of the sensors and sensors
    # resolutions)
    metadata = metadata_container(t,name_converter[fp.sensor_name])

    # Set some globals
    #tree_inspector(t)

    # The hits and tracks 
    dut = hits_plane_accessor(t,metadata.dut_plane)
    ref = hits_plane_accessor(t,metadata.ref_plane)
    tracks = tracks_accessor(t,[0,1,2,3,4],metadata.ref_plane,metadata.dut_plane)
    
    # Process the data
    #start_time = timeit.default_timer()
    wtch = processor(metadata)
    
    nentries = t.GetEntries()
    pointpb = float(nentries)/100.0
    for i,_t in enumerate(t):
        sys.stdout.write("\r\033[1;34mINFO\033[1;m -- Processing file '"+\
                os.path.basename(fname)+" [ "+"\b"+str(int(float(i)/pointpb)+1).rjust(3)+"%]")
        sys.stdout.flush()
        wtch.process(_t,tracks,ref,dut)
    print
    #print timeit.default_timer()-start_time
    print wtch
    foutput = "sensor_maps_{0}_run000{1}.root".format(fp.sensor_name,fp.run_number)
    print "File created at '{0}'".format(foutput)
    wtch.store_histos(foutput)
