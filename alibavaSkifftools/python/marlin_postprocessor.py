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
# 3. Get all the metadata from SPS2017TB class
# 4. Be sure plots are in the LOCAL coordinates of the plane as well, i.e. channel (see the formula to convert in biblio)
# 5. Be sure the binning is taking into account the strip resolution (at least bin_width = 1/4*pitch or so) ??

DEBUG=True

# Be sure use the same alignment files: ct
ALIGN_FILE_FORMAT= "align_cts_ID_{0}.txt"

# -- Take care of the activated branches when reading 
#    the root file
class branch_list(object):
    """
    """
    def __init__(self):
        self._active_br = []
        pass
    def update(self,tree,brlist):
        """De-activate all branches but those
        described and introduced

        Parameters
        ----------
        tree: ROOT.TTree
        brlist: list(sts)
        """
        pass
        #for br in brlist:
        #    self._active_br.append(br)
        #    tree.SetBranchStatus(br,1)
        ## Deactivate all but the active
        #for b in filter(lambda x: not x.GetName() in self._active_br, \
        #        tree.GetListOfBranches()):
        #    tree.SetBranchStatus(b.GetName(),0)
        # Re-do the attributes list (??)

# The global instance
ACTIVE_BRS = branch_list()        

# XXX
# Not sure yet how to deal with that. So far see
# see declaration at sensor_map_production function
# XXX 
RESOLUTION= {} # in mm

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
        # Some useful extra info 
        # REF
        self.ref_pitchX = sensor_name_spec_map[self.ref_name].pitchX
        self.ref_pitchY = sensor_name_spec_map[self.ref_name].pitchY
        self.ref_sizeX  = sensor_name_spec_map[self.ref_name].sizeX
        self.ref_sizeY  = sensor_name_spec_map[self.ref_name].sizeY
        ##  DUT
        self.dut_pitchX = sensor_name_spec_map[self.dut_name].pitchX
        self.dut_pitchY = sensor_name_spec_map[self.dut_name].pitchY
        self.dut_sizeX  = sensor_name_spec_map[self.dut_name].sizeX
        self.dut_sizeY  = sensor_name_spec_map[self.dut_name].sizeY
        # Resolutions: REF and DUT (binary so far)
        _p = 1.0/math.sqrt(12.0)
        self.resolution = { self.ref_plane: self.ref_pitchX*_p/MM , 
                self.dut_plane: self.dut_pitchX*_p/MM }

class dummy_list(list):
    """A dummy list which returns always
    the value in the first element
    """
    def __getitem__(self,i):
        return list.__getitem__(self,0)
    
class alignment_cts(object):
    """Container for the alignment constants
    """
    def __init__(self):
        """Initialize the data-members
        """
        self.iteration = 0
        self.x_offset = 0.0
        self.y_offset = 0.0
        self.turn = 0.0
        self.tilt = 0.0
        self.rot  = 0.0
        self.dz   = 0.0

    def __str__(self):
        """Just the print out of the alignment constants
        """
        from math import degrees
        m  = "iteration: {0}\n".format(int(self.iteration))
        m += "x_offset : {0} mm\n".format(self.x_offset)
        m += "y_offset : {0} mm\n".format(self.y_offset)
        m += "turn     : {0} deg\n".format(degrees(self.turn))
        m += "tilt     : {0} deg\n".format(degrees(self.tilt))
        m += "rot      : {0} deg\n".format(degrees(self.rot))
        m += "dz       : {0} mm\n".format(self.dz)

        return m

    def __iadd__(self,other):
        """Updating the alignment constants

        Parameters
        ----------
        other: alignment_cts instance

        Returns
        -------
        A copy of self
        """
        for attr in self.__dict__.keys():
            prov = getattr(self,attr)+getattr(other,attr)
            setattr(self,attr,prov)
        return self
    
    def __eq__(self,other):
        """Compare two alignment constants, considered
        converged whenever the constants are equal within
        a TOLERANCE

        Parameters
        ----------
        other: alignment_cts instance

        Returns
        -------
        bool: whether or not converged
        """
        # The list of tolerances
        attr_tol = [ ('x_offset', 0.005), ('y_offset', 0.005),\
                ('turn', 0.018), ('tilt', 0.018), ('rot', 0.018),\
                ('dz', 4.0) ]
        # All attributes must be within tolerance
        for attr,tolerance in attr_tol:
            if abs(getattr(self,attr)-getattr(other,attr)) > tolerance:
                return False
        return True
        
    def parse_alignment(self,fobj):
        """Extract the information relative to the alignment
        from a file object

        Parameters
        ----------
        filename: file object
            The alignment file 

        Return
        ------
        A object class containing the relevant alignment
        attributes
        """
        lines=fobj.readlines()
        keywords = [ 'iteration', 'x_offset', 'y_offset', 'rot', 'tilt', 'turn', 'dz' ]
        for i in lines:
            try:
                align_word = filter(lambda x: i.find(x) == 0,keywords)[0]
            except IndexError:
                continue
            setattr(self,align_word,float(i.split(":")[-1]))

class memoize:
    """A simple memozation for the trigonometric 
    functions
    """
    def __init__(self,f):
        self.f = f
        self.memo = {}
    def __call__(self,*args):
        """ Note args = { id,  } 
        """
        if not args in self.memo:
            self.memo[args] = self.f(*args)
        return self.memo[args]

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
    charge: ROOT.std.vector(float)()
        The hit total charge in ADCs counts 
    n_cluster: ROOT.std.vector(int)()
        The total number of elements forming the cluster
    track_link: dict((int,int)
        The index of the hit is related with the track
        produced that hit (if any)
    """
    # The properties of the plane where the hit lives
    # The rotation:=beta (around z-axis), tilt:=alpha (around x-axis) 
    # and turn:=omega (around y-axis) [radians]
    # [A dictionary being the plane id the keys}
    align_constants = {}
    # Activate the update of the alignment constants
    resync = {}
    # vector normal to the plane
    normal   = {}
    
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
        import ROOT

        self.id = planeid
        # ----------------------------------
        global ACTIVE_BRS
        # -- De-activate those branches not needed
        pre_brnames = ["hit_X","hit_Y","hit_Z","hit_XLocal","hit_YLocal","hit_total_charge","hit_Ncluster"]
        brnames = map(lambda x: "{0}_{1}".format(x,self.id),pre_brnames)
        ACTIVE_BRS.update(tree,brnames)
        # ----------------------------------
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
        #self.x_strips = lambda i: (self.x_local[i]+sensor_xsize/2.0)/pitch+0.5
        self.charge = getattr(tree,"hit_total_charge_{0}".format(self.id))
        self.n_cluster = getattr(tree,"hit_Ncluster_{0}".format(self.id))

        # Get the maximum and minimum in x and y (to define the fiducial cuts)
        ROOT.gROOT.SetBatch()
        prv_c = ROOT.TCanvas()
        dummy = tree.Draw("hit_X_{0}>>hdum(100)".format(self.id))
        hdummy = ROOT.gDirectory.Get("hdum")
        self.xmin = hdummy.GetXaxis().GetBinLowEdge(1)
        self.xmax = hdummy.GetXaxis().GetBinUpEdge(100)
        # Dont have this coordinate!!!
        #dummy2 = tree.Draw("hit_Y_{0}>>hdum2".format(self.id))
        #hdummy2 = ROOT.gDirectory.Get("hdum2")
        #self.ymin = hdummy2.GetXaxis().GetBinLowEdge(1)
        #self.ymax = hdummy2.GetXaxis().GetBinUpEdge(hdummy2.GetNbinsX())

        # -- For the modulo definition, let's be sure we don't use any masked region
        #    or channels
        del prv_c

        # Initialize the link index 
        self.track_link = {}

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
        
        # Initialize resync static variable
        if not hits_plane_accessor.resync.has_key(self.id):
            hits_plane_accessor.resync[self.id] = True
        # Check the presence of an alignment file, and then update the alignment
        # constants. Assuming in the telescope frame system
        if hits_plane_accessor.resync[self.id]:
            if DEBUG:
                print
                print "DEBUG ALIGNMENT",self.id
            # Just do it if there weren't initialized previously
            if not hits_plane_accessor.align_constants.has_key(self.id):
                hits_plane_accessor.align_constants[self.id] = alignment_cts()
            filename = ALIGN_FILE_FORMAT.format(self.id)
            try:
                with open(filename) as f:
                    align = alignment_cts()
                    align.parse_alignment(f)
                    # Update the alignment constants
                    hits_plane_accessor.align_constants[self.id] = align
            except IOError:
                # First access, propagate the prealignment from the Marlin 
                # processor here:
                _k = 1
                # Just to be sure we have data
                while self.x_local.size() == 0:
                    dummy = tree.GetEntry(_k)
                    _k+=1
                hits_plane_accessor.align_constants[self.id].x_offset = self.x_local[0]-self.x[0]
            # Not re-synchronized until new order
            if DEBUG:
                print hits_plane_accessor.align_constants[self.id]
            hits_plane_accessor.resync[self.id] = False
        
    @property
    def n(self):
        """The number of hit elements

        Return
        ------
        int
        """
        return self.x.size()
    
    def new_event(self):
        """Call it every time a new event is going to
        be processed. Clean some variables.
        """
        # Clean up the links
        self.track_link = {}
    
    #@staticmethod
    def update_alignment(self,align_inst):
        """Update the alignment constants
        """
        hits_plane_accessor.align_constants[self.id] += align_inst
    
    def get_normal_vector(self):
        """Get the normal vector to the sensor plane, if perfectly 
        aligned it must be (0,0,-1) (note that the plane is defined
        following the opposite direction of the beam)

        Return
        ------
        vector normal: (float,float,float)
        """
        if not hits_plane_accessor.normal.has_key(self.id):
            # Calculate it and momorize
            hits_plane_accessor.normal[self.id] = (-self.cos_tilt(self.id)*self.sin_turn(self.id),\
                    self.sin_tilt(self.id),\
                    -self.cos_tilt(self.id)*self.cos_turn(self.id) )
        return hits_plane_accessor.normal[self.id]
    
    @staticmethod
    @memoize
    def cos_tilt(pid):
        """Memorize cos(tilt) angle

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        float: the cosinos of the tilt (rotation around x-axis)
        """
        from math import cos

        return cos(hits_plane_accessor.align_constants[pid].tilt)
    
    @staticmethod
    @memoize
    def sin_tilt(pid):
        """Memorize sin(tilt) angle

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        float: the sin of the tilt (rotation around x-axis)
        """
        from math import sin
        
        return sin(hits_plane_accessor.align_constants[pid].tilt)
    
    @staticmethod
    @memoize
    def cos_turn(pid):
        """Memorize cos(turn) angle

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        float: the cosinos of the turn (rotation around y-axis)
        """
        from math import cos
        
        return cos(hits_plane_accessor.align_constants[pid].turn)
    
    @staticmethod
    @memoize
    def sin_turn(pid):
        """Memorize sin(turn) angle

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        float: the sin of the turn (rotation around y-axis)
        """
        from math import sin
        
        return sin(hits_plane_accessor.align_constants[pid].turn)
    
    @staticmethod
    @memoize
    def cos_rot(pid):
        """Memorize cos(rot) angle

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        float: the cosinos of the turn (rotation around z-axis)
        """
        from math import cos
        
        return cos(hits_plane_accessor.align_constants[pid].rot)
    
    @staticmethod
    @memoize
    def sin_rot(pid):
        """Memorize sin(rot) angle

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        float: the sin of the turn (rotation around z-axis)
        """
        from math import sin
        
        return sin(hits_plane_accessor.align_constants[pid].rot)
    
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

    def match_isolated(self,track_acc,histos,is_ref=False,res=None):
        """Search for tracks matching (within a resolution) the list 
        of hits in the current event

        Parameters
        ----------
        track_acc: tracks_accessor instance
            The accessor to the telescope tracks branches
        histos: tuple(ROOT.TH)
            A tuple of histograms to be filled within this 
            function. The histograms must be in order:
            - A ROOT.TH2F histogram to store correlation between the 
              hit in the sensor and the hit in the track.
              See `residual_projection` in the processor class 
            - A ROOT.TH1F histogram to store the residuals between the 
              hit in the sensor and the predicted hit for the track
            - A ROOT.TH1F histogram to store the residuals between the 
              hit in the sensor and the predicted hit for the track of
              those passing the cuts (isolation, matching). See `dx_h`
              in the processor class.
            - A ROOT.TProfile to store the variation of the residual in y
              (see last histo) vs. the x-prediction. See `dx_y_h` in the 
              processor class. This histo is used for alignment (rotation)            
            - A ROOT.TProfile to store the variation of the residual in y
              (see last histo) vs. the y-prediction. See `dx_y_h` in the 
              processor class. This histo is used for alignment (tilt)
            - A ROOT.TProfile to store the variation of the residual in y
              (see last histo) vs. the y-slope. See `dx_tx_h` in the 
              processor class. This histo is used for alignment (z-shift)
        is_ref: bool, optional [XXX: Maybe define a data-member in the init)
            Whether or not the current plane is the REFERENCE sensor,
            in that case the isolation criteria should be applied
        res: float (optional)
            The maximum distance to consider a hit matched

        Return
        ------
        hitlist: dict(int,list(int))
            The list of indices of the hits matched with a track 
        """
        from math import sqrt

        if not res:
            res=RESOLUTION[self.id]
        ISOLATION = 0.6*MM

        # Get the histograms
        hcorr,hres,hdx,hdx_finer,hrot,hturn,hdz = histos

        matched_hits = []
        # non isolated tracks and already used: keep them listed 
        # to avoid re-processing
        used_tracks = []
        # Given the hits in that sensor,
        for ihit in xrange(self.n):
            closest = {}
            # Evaluate the distance from the hit to any track (not used already)
            for itrk in filter(lambda i: i not in used_tracks, xrange(track_acc.n)):
                # Isolation: be sure there is no other track surrounding this one
                xpred,ypred,zpred = track_acc.get_point(itrk,self.z[ihit])
                if len(filter(lambda (o_x,o_y,o_z): sqrt(o_x**2.0+o_y**2.0) < ISOLATION,
                    map(lambda other_i:  track_acc.get_point(other_i,zpred), xrange(itrk+1,track_acc.n)))) > 1:
                    # Another track has been found sourrounding this one, 
                    # not isolated, just continue
                    used_tracks.append(itrk)
                    continue
                # Note that the prediction is given in the sensor reference frame
                #xturn,xtilt,xrot,xsen=track_acc.get_point_in_sensor_frame(itrk,self)
                # Fill the alignment histograms
                #closest[abs(self.x_local[ihit]-xsen[0])] = itrk
                closest[abs(self.x[ihit]-xpred)] = itrk
            if len(closest) == 0:
                continue
            # Get the closest track (in x)
            distance_abs,trk_el =sorted(closest.iteritems())[0]
            # and tag it as used, not yet as we don't know if is below the cut!!
            #used_tracks.append(trk_el)
            # Link the hit with the track
            # self.track_link[ihit] = trk_el
            # Get again  the info in the reference frame of the sensor
            # XXX --- SEEMS not to need all the partial results
            #rturn,rtitl,rrot,rsensor = track_acc.get_point_in_sensor_frame(trk_el,self)
            rsensor = track_acc.get_point(trk_el,self.z[ihit])
            # And fill some histograms
            # Fill correlation histograms
            ### hcorr.Fill(self.x_local[ihit],rsensor[0])
            hcorr.Fill(self.x[ihit],rsensor[0])
            # -- resolution histogram
            ## -- dx = self.x_local[ihit]-rsensor[0]
            dx = self.x[ihit]-rsensor[0]
            hres.Fill(dx)
            # Alignment histos: x-offset
            # -- Pre-alignment, rot, tilt an turn but no shift
            ##hdx.Fill(self.x_local[ihit]-rrot[0])
            hdx.Fill(self.x[ihit]-rsensor[0])
            # -- finer
            hdx_finer.Fill(dx)
            # Fill only thos passing the cuts
            if distance_abs < res:
                # Fill the matched hits
                matched_hits.append( (ihit,trk_el) )
                # Tag the cut as used
                used_tracks.append(trk_el)
                # and link the hit with the track
                self.track_link[ihit] = trk_el
                ## Fill the alignment histograms after being accepted
                #  as providing from the hit
                ## -- rot
                hrot.Fill(rsensor[1],dx)
                # -- turn
                hturn.Fill(rsensor[0],dx)
                ## -- dz 
                hdz.Fill(track_acc.dxdz[trk_el],dx)
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
        # ----------------------------------
        global ACTIVE_BRS
        # -- De-activate those branches not needed
        pre_brnames = ["X","Y","Z","XLocal","YLocal","Z","index"]
        brnames = map(lambda x: "{0}_{1}_{2}".format(prefix,x,self.id),pre_brnames)
        ACTIVE_BRS.update(tree,brnames)
        # ----------------------------------
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
    
    def new_event(self):
        """Call it every time a new event is going to
        be processed. Clean some variables.
        Do nothing, just in case we need to extend this
        class.
        """
        pass
    
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
 
class memoize_hasher(object):
    """A simple memoization using a hash of 
    the inputs
    Use the track Id and a hash from the hit object
    as keys. 
    Inspired from
    http://code.activestate.com/recipes/577452-a-memoize-decorator-for-instance-methods/
    """
    def __init__(self,f):
        self.f = f
    def __get__(self,obj,objtype=None):
        """
        """
        from functools import partial
        if obj is None:
            return self.f
        return partial(self,obj)
    def __call__(self,*args):
        """Only need the z-plane, all the other stuff is the same
        """
        # Clean call 
        obj = args[0]
        if args[1] == "clean":
            self._clean(obj)
            return
        # -- REGULAR CALL
        try:
            cache = obj.__cache
        except AttributeError:
            cache = obj.__cache = {}
        if hasattr(args[2],'z'):
            args1 = int(args[2].z[0])
        else:
            args1 = int(args[2])
        key = (self.f,args[1],args1)
        try:
            res = cache[key]
        except KeyError:
            res = cache[key] = self.f(*args)
        return res
    def _clean(self,obj):
        """Clean the cache
        """
        obj.__cache ={}

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
        global ACTIVE_BRS
        # -- De-activate those branches not needed
        brnames = ["trk_hit_fit_X_{0}".format(dut_plane),"trk_hit_fit_X_{0}".format(ref_plane),\
                "trk_refPoint_X","trk_refPoint_Y","trk_refPoint_Z","trk_dxdz","trk_dydz"]
        ACTIVE_BRS.update(tree,brnames)
        # --------------------------------
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

    def new_event(self):
        """Call it every time a new event is going to
        be processed. Clean some variables.
        """
        # -- Clean the cache
        self.get_point("clean")
        self.get_point_in_sensor_frame("clean")
    
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
    
    @memoize_hasher
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
        (float,float,float): Predicted position at the z-plane
        """
        t = z-self.z0[i]
        return (self.x0[i]+self.dxdz[i]*t,self.y0[i]+self.dydz[i]*t,z)
    
    @memoize_hasher
    def get_point_in_sensor_frame(self,i,hit):
        """Give the predicted position (x,y) for the i-track using
        a transformation into the sensor frame (passive).
        Suitable when the sensors are rotated or titled

        Parameters
        ----------
        i: int
            The index of the track
        hit: hits_plane_accessor
            The hit instance
        
        Return
        ------
        (nt1,nt2,nt3,nt4): being nti=(float,float,float)
            The predicted position by the track in the sensor
            reference frame, after apply one transformation
            1: turn
            2: tilt 
            3: rot 
            4: sensor reference plane
        """
        # XXX Memoize using the itrk and the hit hash
        # Obtain the normal vector of the plane (defined in opposite direction to the beam)
        # With perfect alignment: (0,0,-1)
        n = hit.get_normal_vector()
        # Intersect the track with the sensor plane (to obtain the diferent z-values depending 
        # the x or y). Plane equation n (r-r_{p}),
        #  - the plane vector: n
        #  - the assume plane position r_{p}= (0,0,z_p-z0)
        #  - the predicted point       r = (x0+tx(z-z0), y0+ty(z-z0), z-z0)
        # Then, t=(z-z0) is obtained for an inclined plane (therefore 
        # will give different z-values depending on the x,y, instead of being at the fixed z_p
        # Note zc = z-z0 (the t-parameter equivalent but using a tilt-rot or turned plane)
        zc = (n[2]*(hit.z[0]-self.z0[i])-n[1]*self.y0[i]-n[0]*self.x0[i])/(n[0]*self.dxdz[i]+n[1]*self.dydz[i]+n[2])
        ### XXX --- print n,t,hit.z[0],hit.z[0]-self.z0[i]
        # Transfrom into the sensor plane frame (passive transformation)
        r_at_sens = (self.x0[i]+self.dxdz[i]*zc,self.y0[i]+self.dydz[i]*zc,zc)
        # Distance away from the sensor frame:
        #   the z in the telescope plane for the current track (zc), minus
        #   the z position of the sensor plane (hit.z[0]) corrected by the alignment cte
        dzc = zc  - (hit.z[0]+hit.align_constants[hit.id].dz)
        #print dzc, zc, self.z0[i], (hit.z[0]-hit.align_constants[hit.id].dz)
        ## print dzc,self.z0[i],zc

        ### XXX --- print r_at_sens
        # Euler angles
        # -- First: turn
        x1 = (r_at_sens[0]*hit.cos_turn(hit.id)-dzc*hit.sin_turn(hit.id),\
                r_at_sens[1],r_at_sens[0]*hit.sin_turn(hit.id)+dzc*hit.cos_turn(hit.id))
        # -- Second: tilt
        x2 = (x1[0],x1[1]*hit.cos_tilt(hit.id)+x1[2]*hit.sin_tilt(hit.id),-x1[1]*hit.sin_tilt(hit.id)+x1[2]*hit.cos_tilt(hit.id))
        # -- third: rot
        x3 = (x2[0]*hit.cos_rot(hit.id)+x2[1]*hit.sin_rot(hit.id),-x2[0]*hit.sin_rot(hit.id)+x2[1]*hit.cos_rot(hit.id),x2[2])
        # print r_at_sens,x3
        # -- And finally, the shift
        # Note that if x_offset is obtained from x_local - x 
        #if abs(hit.align_constants[hit.id].x_offset) < 1e-9:
        #    hit.align_constants[hit.id].x_offset = hit.x_local[i]
        return (x1,x2,x3,(x3[0]+hit.align_constants[hit.id].x_offset,x3[1]+hit.align_constants[hit.id].y_offset,x3[2]))

class processor(object):
    """Results extractor. Main class where statistics and histograms 
    are defined and filled. 

    Attributes
    ----------
    total_events: int
        The number of processed events
    events_with_dut: int
        The number of processed events containing at least a 
        hit in the DUT plane
    dut_hits: int
        The number of hits measured at the DUT sensor
    events_with_ref: int
        The number of processed events containing at least a 
        hit in the REF plane
    ref_hits: int
        The number of hits measured at the REF sensor
    events_with_dut_no_tracks: int
        The number of processed events containing at least a 
        hit in the DUT plane but with 0 reconstructed tracks
    events_with_ref_no_tracks: int
        The number of processed events containing at least a 
        hit in the REF plane but with 0 reconstructed tracks
    events_with_ref_dut_tracks: int
        The number of processed events containing at least a 
        hit in the REF, and at least a hit in the DUT and 
        at least a reconstructed track
    residual_projection: dict(int,ROOT.TH1F)
        Histograms of the difference between the hit measured at the REF/DUT 
        and the predicted point at the REF/DUT plane by a 
        telescope track. If an event has more than one track, the one used is
        the one with minimum difference between the prediction and the measurement
    residual_sensor: ROOT.TH2F
        Histogram of the difference between the hit predicted at the DUT and the
        hit predicted at the REF by the same, and matched to both sensors, track
    hcharge: dict(int, ROOT.TH2F)
        Histograms, per sensor, of the deposited charge in the track predicted 
        position at the plane of the sensor
    hcorr_trkX: dict(int,ROOT.TH2F)
        Histograms per sensor, of the difference between the measured hit and the
        predicted hits of the closest (defined as the x-distance) track of the event
    hcorrX(Y): dict(int,ROOT.TH2F)
        Histograms per sensor, of the difference between the track predicted position
        at the DUT and at the REF, using the same matched track (meaning that the same
        track is compatible with the REF and the DUT measured hits).
    """
    # Some static members to be used for the aligment
    # - whether or not the dz alignment should be performed
    #align_allow_dz = {}
    # - whether or not the rot alignment should be performed
    #align_allow_rot = {}

    def __init__(self,minst):
        """ XXX MISSING DOC
        Parameters
        ----------
        minst: metadata_container instance
        """
        import ROOT

        # Get some useful info to build the histograms
        self.dut_plane = minst.dut_plane
        self.ref_plane = minst.ref_plane

        # Asumming strips for the ref (so recall is going 
        # to be used as 2 x pitch)
        self.pitchX = { minst.dut_plane: minst.dut_pitchX, minst.ref_plane: minst.ref_pitchX }
        ## --- XXX Simulate pixels
        self.pitchY = { minst.dut_plane: minst.dut_pitchY, minst.ref_plane: minst.ref_pitchY }
        self.pitchY = self.pitchX
        # -- Check if is 3-D or strips
        #if minst.dut_name.lower().find("lgad") == -1:
        #    # simulate pixels
        #    self.pitchY[minst.dut_plane] = minst.dut_pitchX
        #else:
        #    # As is going to be use as 2 x pitch
        #    self.pitchY[minst.dut_plane] = 0.5*minst.dut_pitchY
        self.sizeX = { minst.dut_plane: minst.dut_sizeX, minst.ref_plane: minst.ref_sizeX }
        # Be careful than iLGAD and LGAD don't have the proper dimensions of the sensors
        if minst.dut_name.lower().find("ilgad") == 0:
            self.sizeX[minst.dut_plane]= 7.2*MM
        elif minst.dut_name.lower().find("lgad") == 0:
            self.sizeX[minst.dut_plane]= 5.2*MM
        self.sizeY = { minst.dut_plane: minst.dut_sizeY, minst.ref_plane: minst.ref_sizeY }
        # useful for the histos
        # half size plus some extra to keep alignment factors
        sxdut = self.sizeX[minst.dut_plane]/2.0+1.5
        sydut = self.sizeY[minst.dut_plane]/2.0+1.5
        sxref = self.sizeX[minst.ref_plane]/2.0+1.5
        syref = self.sizeY[minst.ref_plane]/2.0+1.5

        # some info numbers
        self.total_events    = 0
        self.events_with_dut = 0
        self.dut_hits        = 0
        self.events_with_ref = 0
        self.ref_hits        = 0
        self.events_with_dut_no_tracks = 0
        self.events_with_ref_no_tracks = 0
        self.events_with_ref_dut_tracks = 0
        self.number_matched_hits = {}
        # Alingment histos
        # -- the shift in x
        self.dx_h = { minst.dut_plane: ROOT.TH1F("dx_dut"," ; x_{DUT}-x_{trk}^{pred} [mm]; Entries",100,-sxdut,sxdut),
                minst.ref_plane: ROOT.TH1F("dx_ref"," ; x_{REF}-x_{trk}^{pred} [mm]; Entries",100,-sxref,sxref) }
        # -- finer alignment
        self.dx_finer_h = { minst.dut_plane: ROOT.TH1F("dx_finer_dut"," ; x_{DUT}-x_{trk}^{pred} [mm]; Entries",100,-0.5,0.5),
                minst.ref_plane: ROOT.TH1F("dx_finer_ref"," ; x_{REF}-x_{trk}^{pred} [mm]; Entries",100,-0.5,0.5) }
        # -- the rot (around z-axis)
        self.dx_y_h = { minst.dut_plane: ROOT.TProfile("dx_y_dut"," ;y_{trk}^{pred} [mm];#Deltax_{DUT} [mm]",50,-6.0,6.0,-0.2,0.2),
                minst.ref_plane: ROOT.TProfile("dx_y_ref"," ;y_{trk}^{pred} [mm];#Deltax_{REF} [mm]",50,-6.0,6.0,-0.2,0.2) }
        # -- the turn (around y-axis)
        self.dx_x_h = { minst.dut_plane: ROOT.TProfile("dx_x_dut"," ;x_{trk}^{pred} [mm];#Deltax_{DUT} [mm]",50,-6.0,6.0,-1.2,1.2),
                minst.ref_plane: ROOT.TProfile("dx_x_ref"," ;x_{trk}^{pred} [mm];#Deltax_{REF} [mm]",50,-6.0,6.0,-0.2,0.2) }
        self.dx_tx_h = { minst.dut_plane: ROOT.TProfile("dx_tx_dut"," ;#theta_{x}^{trk};#Deltax_{DUT} [mm]",50,-1e-4,1e-4,-0.2,0.2),
                minst.ref_plane: ROOT.TProfile("dx_tx_ref"," ;#theta_{x}^{trk};#Deltax_{REF} [mm]",50,-1e-4,1e-4,-0.2,0.2) }
        
        self._alignment_histos = self.dx_h.values()+self.dx_finer_h.values()+self.dx_y_h.values()+self.dx_x_h.values()+self.dx_tx_h.values()

        # Some statistic histograms
        # -------------------------
        self.residual_projection = { minst.dut_plane: ROOT.TH1F("res_projection_dut"," ; x_{DUT}-x_{trk}^{pred} [mm]; Entries",200,-3.5*MM,3.5*MM),
                minst.ref_plane: ROOT.TH1F("res_projection_ref"," ; x_{REF}-x_{trk}^{pred} [mm]; Entries",200,-3.5*MM,3.5*MM) }
        # Residuals between REF-DUT (matched tracks family)
        self.residual_sensor = ROOT.TH2F("res_sensor_projection"," ; x_{DUT}^{pred}-x_{REF}^{pred} [mm];y_{DUT}^{pred}-y_{REF}^{pred}",\
                100,-0.04*MM,0.04*MM,100,-3.5*MM,3.5*MM) 
        self.hcharge = { minst.dut_plane: ROOT.TProfile2D("charge_map_dut" ,";x_{DUT}^{pred} [mm];y_{DUT}^{pred} [mm];"\
                        "<charge cluster> [ADC]", 300,-sxdut,sxdut,300,-sydut,sydut),
                minst.ref_plane: ROOT.TProfile2D("charge_map_ref",";x_{REF}^{pred} [mm];y_{REF}^{pred} [mm];"\
                        "<charge cluster> [ADC]", 300,-sxref,sxref,300,-syref,syref) }
        self.hhitmap = { minst.dut_plane: ROOT.TH2F("hitmap_dut" ,";x_{DUT}^{pred} [mm];y_{DUT}^{pred} [mm];"\
                        "Entries", 300,-sxdut,sxdut,300,-sydut,sydut),
                minst.ref_plane: ROOT.TH2F("hitmap_ref",";x_{REF}^{pred} [mm];y_{REF}^{pred} [mm];"\
                        "Entries", 300,-sxref,sxref,300,-syref,syref) }
        # -- Module (2 x pitch X)
        self.hcluster_size_mod = { minst.dut_plane: ROOT.TProfile2D("cluster_size_mod_dut" ,";mod(x_{trk})_{2xpitch} [mm];mod(y_{trk})_{2xpitch}"\
                        "[mm];<cluster size>", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TProfile2D("cluster_size_mod_ref" ,";mod(x_{trk})_{2xpitch} [mm];mod(y_{trk})_{2xpitch} [mm];"\
                        "<cluster size>", 40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }
        self.hcharge_mod = { minst.dut_plane: ROOT.TProfile2D("charge_mod_dut",";mod(x_trk})_{2xpitch} [mm]; mod(y_{trk})_{2xpitch} [mm];"\
                        "<cluster charge> [ADC]", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TProfile2D("charge_mod_ref",";mod(x_trk})_{2xpitch} [mm]; mod(y_{trk})_{2xpitch} [mm];"\
                        "<cluster charge> [ADC]", 40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }
        self.hhitmap_mod = { minst.dut_plane: ROOT.TH2F("hitmap_mod_dut",";mod(x_trk})_{2xpitch} [mm]; mod(y_{trk})_{2xpitch} [mm];"\
                        "Entries", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TH2F("hitmap_mod_ref",";mod(x_trk})_{2xpitch} [mm]; mod(y_{trk})_{2xpitch} [mm];"\
                        "Entries", 40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }
        # -- efficiency mod
        self.heff_mod = ROOT.TProfile2D("eff_mod",";mod(x_trk})_{2xpitch} [mm]; mod(y_{trk})_{2xpitch} [mm];"\
                        "efficiency", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM,-1,2)
        # Correlations: 
        # -- Sensors-tracks
        self.hcorr_trkX = { minst.dut_plane: ROOT.TH2F("corr_trkX_dut",";x_{DUT} [mm]; x_{pred}^{trk} [mm]; Entries",200,-sxdut,sxdut,200,-sydut,sydut),
                minst.ref_plane: ROOT.TH2F("corr_trkX_ref",";x_{REF} [mm]; x_{pred}^{trk} [mm]; Entries",200,-sxdut,sxdut,200,-sydut,sydut)}
        # -- DUT-REF
        self.hcorrX = ROOT.TH2F("corrX_dut_ref",";x_{DUT} [mm]; x_{REF} [mm]; Entries",100,-sxdut,sxdut,100,-sxref,sxref)
        self.hcorrY = ROOT.TH2F("corrY_dut_ref",";y_{DUT} [mm]; y_{REF} [mm]; Entries",100,-sydut,sydut,100,-syref,syref)
        # efficiency
        self.heff = ROOT.TProfile2D("eff_map",";x_{trk}^{DUT} [mm]; y^{DUT}_{trk} [mm];#varepsilon",50,-sxdut,sxdut,50,-sydut,sydut)
        self.heff_entries = ROOT.TH2F("eff_entries",";x_{trk}^{DUT} [mm]; y^{DUT}_{trk} [mm];#varepsilon",50,-sxdut,sxdut,50,-sydut,sydut)

        self._allhistograms = self.residual_projection.values()+[self.residual_sensor]+\
                self.hcharge.values()+self.hhitmap.values()+self.hcorr_trkX.values()+\
                self.hcluster_size_mod.values()+self.hcharge_mod.values()+self.hhitmap_mod.values()+[self.heff_mod]+\
                self._alignment_histos+\
                [self.hcorrX,self.hcorrY,self.heff,self.heff_entries]
        dummy=map(lambda h: h.SetDirectory(0),self._allhistograms)

        # The alignment constants
        self.alignment = {}

    #def draw_all_in_canvas():
    #    """
    #    """
    #    import ROOT


    def store_histos(self,filename):
        """Actually write the defined histograms to 
        a root file

        Parameters
        ----------
        filename: str
            The name of file to store the results
        """
        import ROOT
        # Store all the histograms
        f=ROOT.TFile.Open(filename,"RECREATE")
        for h in self._allhistograms:
            h.Write()
        f.Close()
    
    def update_alignment_file(self):
        """Extract the information relative to the alignment
        and dump it into a file

        Explanation
        """
        from math import sin,cos,sqrt,asin
        MIN_ENTRIES = 1000

        # Just extract all the static info of the class,
        for pl_id in self.alignment.keys():
            # Update the counter
            self.alignment[pl_id].iteration += 1
            align_file_content='iteration: {0}\n'.format(int(self.alignment[pl_id].iteration))
            filename = ALIGN_FILE_FORMAT.format(pl_id)
            # check the new offset
            new_x_offset = self.alignment[pl_id].x_offset
            # Some checks first: do a gross alignment first (before the shift)
            if self.dx_h[pl_id].GetEntries() > MIN_ENTRIES:
                new_x_offset = get_x_offset(self.dx_h[pl_id])
                # Finer alignment -- including the shift, therefore to be added
                # to the old offset
                if self.alignment[pl_id].iteration > 1 \
                        and abs(new_x_offset-self.alignment[pl_id].x_offset) < 0.1:
                    new_x_offset = self.alignment[pl_id].x_offset+get_x_offset(self.dx_finer_h[pl_id],-0.1,0.1)
                # -- The x-offset 
                align_file_content += "x_offset: {0}\n".format(new_x_offset)
            else:
                align_file_content += "x_offset: {0}\n".format(self.alignment[pl_id].x_offset)
            # -- The y-offset (We don't have y-coordinate so far)
            align_file_content += "y_offset: {0}\n".format(0.0)
            # -- The rotation
            # XXX --- Maybe if the change is less than 0.1 degrees (or some % of the previous one)
            #         There is no more alignment in here
            if self.dx_y_h[pl_id].GetEntries() > MIN_ENTRIES:
                new_align_rot = self.alignment[pl_id].rot-get_linear_fit(self.dx_y_h[pl_id])
                # Check convergence, if the change are below 0.5 degree, just 
                # converged (not checked in the first iteration)
                if self.alignment[pl_id].iteration > 1 \
                        and abs(self.alignment[pl_id].rot-new_align_rot) < 0.0087:
                    #processor.align_allow_rot = False
                    align_file_content += "rot: {0}\n".format(self.alignment[pl_id].rot)
                else:
                    align_file_content += "rot: {0}\n".format(new_align_rot)
            else:
                align_file_content += "rot: {0}\n".format(self.alignment[pl_id].rot)
            # -- The tilt(??)
            #align_file_content[i].append("tilt: {0}\n".format(self.alignment[pl_id]tilt+get_linear_fit(self.dy_y[pl_id])))
            # -- The turn (only evaluate it if greater than 5 degrees: 0.087 rad.
            if abs(self.alignment[pl_id].turn) > 0.018 \
                    and self.dx_x_h[pl_id].GetEntries() > MIN_ENTRIES:
                new_turn = get_linear_fit(self.dx_x_h[pl_id])
                align_file_content += "turn: {0}\n".format(self.alignment[pl_id].turn+new_turn/sin(self.alignment[pl_id].turn))
            else:
                # --- XXX NOT TRUSTABLE XXX
                # Check the aproximation (asuming the angle is small): sin**2
                #new_turn = get_linear_fit(self.dx_x_h[pl_id])
                #sign_turn = int(new_turn/abs(new_turn))
                #align_file_content += "turn: {0}\n".format(self.alignment[pl_id].turn+sign_turn*asin(sqrt(abs(new_turn))))
                # --- XXX NOT TRUSTABLE XXX
                align_file_content += "turn: {0}\n".format(self.alignment[pl_id].turn)
            # -- the z-shift: Note the SPS beam hardly divergent (less than 0.1 mrad),
            #    therefore do not trust too much in this correction, not stable. 
            #    We trust in the initial measurement, if the diference is higher than a couple of mm, 
            #    deactivate it
            if self.dx_tx_h[pl_id].GetEntries() > MIN_ENTRIES:
                new_dz = get_linear_fit(self.dx_tx_h[pl_id],-1e-3,1e-3)
                # -- XXX Guard to avoid not estable results XXX
                if abs(new_dz) > 3.0:
                    align_file_content += "dz: {0}\n".format(self.alignment[pl_id].dz)
                else:
                    align_file_content += "dz: {0}\n".format(self.alignment[pl_id].dz+new_dz)
            else:
                align_file_content += "dz: {0}\n".format(self.alignment[pl_id].dz)
            # Now do the fits to the histograms
            with open(filename,"w") as f:
                f.write(align_file_content)
            # And prepare the alignment constants to be updated 
            # in the next iteration
            hits_plane_accessor.resync[pl_id] = True
        # Update the alignment loop here?

    def store_alignment(self,filename):
        """Actually write the alignment histograms to 
        a root file

        Parameters
        ----------
        filename: str
            The name of file to store the results
        """
        import ROOT
        # Store all the histograms
        f=ROOT.TFile.Open(filename,"RECREATE")
        for h in self._alignment_histos:
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
        
    def fill_statistics_matched(self,matched_hits):
        """Fill some statistic related with the number of 
        matched-to-track hits

        Parameters
        ----------
        matched_hits:  dict(int,[])
            The list of matched hits 
        """
        for plid,thelist in matched_hits.iteritems():
            if not self.number_matched_hits.has_key(plid):
                self.number_matched_hits[plid] = len(thelist)
            else:
                self.number_matched_hits[plid] += len(thelist)


    def process(self,t,trks,refhits,duthits,is_alignment=False):
        """Fill the relevant counters and the histograms

        Parameters
        ----------
        t: ROOT.TTree
            The tree
        trks: .tracks_ancessor
            The tracks
        refhits: list(.hits_plane_accessor )
            The list of measured hits at the REF sensor
        duthits: list(.hits_plane_accessor )
            The list of measured hits at the DUT sensor
        is_alignment: bool, optional
            Whether if the run is an alignment run, therefore
            no need for filling extra histograms
        """
        import math
        
        # Get the number of hits in the DUT and REF plane
        self.fill_statistics(duthits.n,refhits.n,t.ntracks)

        # Start new event
        trks.new_event()
        refhits.new_event()
        duthits.new_event()

        # Fill histograms for the DUT/REF
        # -- Charge map: using a matched track to the sensor, in order to
        # find the missing Y
        # -- Efficiency map: matched track to the reference, and extrapolate to
        # the DUT plane (total histogram), match to the DUT  (pass histogram)
        
        # - a dict to get the sensors by its plane ID
        sensor_hits = { refhits.id: refhits, duthits.id: duthits }
        
        # Update the sensors with the alignment results
        #for pl_id,hits in sensor_hits.iteritems():
        #    hits.update_alignment(self.alignment[pl_id])

        # Prepare the ntuple of histograms for DUT and REF
        histos_dut = (self.hcorr_trkX[duthits.id],self.residual_projection[duthits.id],self.dx_h[duthits.id],\
                self.dx_finer_h[duthits.id],self.dx_y_h[duthits.id],self.dx_x_h[duthits.id],self.dx_tx_h[duthits.id])
        histos_ref = (self.hcorr_trkX[refhits.id],self.residual_projection[refhits.id],self.dx_h[refhits.id],\
                self.dx_finer_h[refhits.id],self.dx_y_h[refhits.id],self.dx_x_h[refhits.id],self.dx_tx_h[refhits.id])
        histos = { refhits.id: histos_ref, duthits.id: histos_dut }
        # Get the matched hits { sensorID: [ (hit_index,track_index),. ..] , }
        #matched_hits = trks.matched_hits(duthits,refhits,self.hcorr_trkX)
        matched_hits = { duthits.id: duthits.match_isolated(trks,histos_dut),
                refhits.id: refhits.match_isolated(trks,histos_ref,True) }

        # If alignment, just return here, trying to avoid extra processing time
        if is_alignment:
            # Update the alignment constants 
            self.alignment = dict(map(lambda (i,d): (i,d.align_constants[i]) ,sensor_hits.iteritems()))
            # Get an estimation of the sensor efficiency ?
            self.fill_statistics_matched(matched_hits)
            return

        predictions = {}
        # Filling histograms
        for sensorID,pointlist in matched_hits.iteritems():
            hits = sensor_hits[sensorID]
            for hit_el,trk_index in pointlist:
                # Remember: it_el=index of the measured hit at the corrent sensor (sensorID)
                #           trk_index= index of the track
                xpred,ypred,zpred = trks.get_point(trk_index,hits.z[0])
                # Event and charge maps (using track predictions)
                # --- Some histograms should use the measured values (at least when possible)
                self.hcharge[sensorID].Fill(xpred,ypred,hits.charge[hit_el])
                self.hhitmap[sensorID].Fill(xpred,ypred)
                # Modulo (Assuming the center of the sensor: -0.5*pitch 
                # XXX provisional
                xmod = (32.+xpred)%(2.0*self.pitchX[sensorID])
                ymod = (32+ypred)%(2.0*self.pitchY[sensorID])
                if abs(xmod) > 2.0*self.pitchX[sensorID] \
                        or abs(ymod) > 2.0*self.pitchY[sensorID]:
                    raise RuntimeError("WHAT THE FUCK!??!?")
                self.hcluster_size_mod[sensorID].Fill(xmod*UM,ymod*UM,hits.n_cluster[hit_el])
                self.hcharge_mod[sensorID].Fill(xmod*UM,ymod*UM,hits.charge[hit_el])
                self.hhitmap_mod[sensorID].Fill(xmod*UM,ymod*UM)
                ## INCLUDE IT IN HERE
        # Matching beetwen DUT and REF: hits with the same track ID
        # Loop over all the measured hits in the REF, matching to a
        # track, using the track as common element in the DUT-REF match
        # XXX Just using those within the fiducial cut of the DUT XXX
        for ihit_ref,itrk in matched_hits[refhits.id]:
            # Get the track predicted points at the DUT plane
            xpred_dut,ypred_dut,zpred_dut = trks.get_point(itrk,duthits.z[0])
            # Only within fiducial  (note that we don't have coordinate Y)
            if xpred_dut > duthits.xmax or xpred_dut < duthits.xmin:
                continue
            # Modulo pitch
            xmod = (32.0+xpred_dut)%(2.0*self.pitchX[duthits.id])
            ymod = (32.0+ypred_dut)%(2.0*self.pitchY[duthits.id])
            # Get if any of the 
            dut_match = filter(lambda (ihit_dut,itrk_dut): itrk_dut == itrk,matched_hits[duthits.id])
            any_match = int(len(dut_match)>0)
            # Efficiency
            self.heff_mod.Fill(xmod*UM,ymod*UM,any_match)
            self.heff.Fill(xpred_dut,ypred_dut,any_match)
            if len(dut_match) > 0:
                self.heff_entries.Fill(xpred_dut,ypred_dut)
                # Get the track predicted points at the REF plane
                xpred_at_ref,ypred_at_ref,zpred_at_ref = trks.get_point(itrk,refhits.z[0])
                # And some correlation plots of the predicted values between
                # the DUT and REF sensors
                # Note that the expected gaussian could be not centered at zero,
                # in that case the difference will be related with the misaligned between them
                self.hcorrX.Fill(xpred_dut,xpred_at_ref)
                # Not measured elements at Y, so using the predictions
                self.hcorrY.Fill(ypred_dut,ypred_at_ref)
                ## Residuals only for the DUT-REF matched tracks
                self.residual_sensor.Fill(xpred_dut-xpred_at_ref,ypred_dut-xpred_at_ref)
            elif len(dut_match) > 1:
                raise RuntimeError("This must never happens!! Contact developer, error 3E43")

    def get_raw_sensors_efficiency(self):
        """Summarize the efficiency of the sensors:
        eff = Number of matched hits/number of events
        """
        m = ""
        for plid,n in self.number_matched_hits.iteritems():
            if plid == self.dut_plane:
                sensor_evts = self.events_with_dut
                sensor_name = "DUT"
            else:
                sensor_evts = self.events_with_ref
                sensor_name = "REF"
            m+= "{0} efficiency           : {1:.2f}%\n".format(sensor_name,float(sensor_evts)/float(self.total_events)*100.)
            m+= "{0} efficiency (MATCHED) : {1:.2f}%\n".format(sensor_name,float(n)/float(self.total_events)*100.)
        return m

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

def get_x_offset(h,xmin=-2.0,xmax=2.0):
    """XXX 
    """
    import ROOT
    ROOT.gROOT.SetBatch()
    
    cns=ROOT.TCanvas()

    gbg = ROOT.TF1("gbg","[0]*exp(-0.5*((x-[1])/[2])**2.0)+[3]",xmin,xmax)
    # Set the initial values
    # -- Amplitude
    gbg.SetParameter(0,h.GetMaximum())
    # -- Mean
    peak = h.GetBinCenter(h.GetMaximumBin())
    gbg.SetParameter(1,peak)
    # -- Sigma
    gbg.SetParameter(2,h.GetBinWidth(1))
    # -- Background: get the guess from far the peak point
    gbg.SetParameter(3,h.GetBinContent(h.FindBin(peak-0.4)))
    ## Do the fit
    status = h.Fit(gbg,"SQR","")
    if status == 0:
        # XXX
        print "KKITA!!!! FAILED THE X-OFFSET FIT"
    new_align_x = gbg.GetParameter(1)
    
    return new_align_x

def get_linear_fit(h,xmin=-4.0,xmax=4.0):
    """XXX 
    """
    import ROOT
    ROOT.gROOT.SetBatch()
    cns=ROOT.TCanvas()

    gbg = ROOT.TF1("linear_fit","pol1",xmin,xmax)
    ## Do the fit
    status = h.Fit(gbg,"SQR")
    if status == 0:
        # XXX
        print "KKITA!!!! FAILED THE LINEAR FIT"
    return gbg.GetParameter(1)

def sensor_alignment(fname,verbose):
    """Performs the alignment of the DUT and REF sensors
    with respect to the Telescope, using the tracks of 
    the telescope and assuming the telescope is already 
    aligned. The coordinate Y of both sensors are considered
    to be unmeasurable, but is assumed to be assigned by the
    track (in its first iteration?)

    Parameters
    ----------
    fname: str
        The name of the ROOT ntuple produced by the EUTelTreeCreator
        processor (part of the EUTelescope package, which uses the 
        Marlin framework). The filename MUST follow the standard 
        notation defined through the class 
        `alibavaSkifftools.SPS2017TB_metadata.filename_parser`
    verbose: bool
        Whether or not to print out the alignment parameters

    Returns
    -------
    status: int
        The code status describing the outcome of the alignment
        process: 0 means that the parameters describing the alignment 
        are below some predefined values, considering the sensor
        aligned. The value -1 is sent otherwise.

    Raises
    ------
    IOError
        Whenever the ROOT file is not present
    """
    import ROOT
    from math import sin

    # Get the data, build clusters, do the actual analysis, 
    # but only the first 50k events
    proc_inst = sensor_map_production(fname,entries_proc=60000,alignment=True)
    # Considered aligned if the aligned constants doesn't change (within some 
    # tolerance): then extract the alignment constants before
    aligned_sensors = 0
    for pl_id,current_align in proc_inst.alignment.iteritems():
        align_f = ALIGN_FILE_FORMAT.format(pl_id)
        old_align = alignment_cts()
        try:
            with open(align_f) as f:
                old_align.parse_alignment(f)
        except IOError:
            # not do anything at the begining
            continue
        if old_align == current_align:
            aligned_sensors += 1
    # Update the alignment files
    proc_inst.update_alignment_file()
    print proc_inst.get_raw_sensors_efficiency()
    # Check if all sensors are aligned: XXX Not correct  ?
    return (aligned_sensors == len(proc_inst.alignment.keys()))


def sensor_map_production(fname,entries_proc=-1,alignment=False):
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
    entries_proc: int, optional
        The number of entries to process, if -1 all entries available
    alignment: bool, optional
        Whether this run is used to align or not. Effectively means
        that less histograms are going to be fill and less entries processed

    Return
    ------
    A processor instance (if alignment mode)

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
    from .SPS2017TB_metadata import sensor_name_spec_map as specs
    # probably provisional XXX
    global RESOLUTION


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

    # FIXME provisional XXX : 2 x pitch for dut, same for ref?. Check the efficiency
    RESOLUTION= { metadata.dut_plane: 2.0*specs[name_converter[fp.sensor_name]].pitchX,\
            metadata.ref_plane : 2.0*specs["REF_0_b1"].pitchX } # in mm

    # Set some globals
    #tree_inspector(t)

    # The hits and tracks 
    dut = hits_plane_accessor(t,metadata.dut_plane)
    ref = hits_plane_accessor(t,metadata.ref_plane)
    tracks = tracks_accessor(t,[0,1,2,3,4],metadata.ref_plane,metadata.dut_plane)
    
    # Process the data
    #start_time = timeit.default_timer()
    wtch = processor(metadata)
    
    if entries_proc == -1:
        nentries = t.GetEntries()
    else:
        nentries = int(entries_proc)
    pointpb = float(nentries)/100.0
    for i,_t in enumerate(t):
        if i > nentries:
            break
        sys.stdout.write("\r\033[1;34mINFO\033[1;m -- Processing file '"+\
                os.path.basename(fname)+" [ "+"\b"+str(int(float(i)/pointpb)+1).rjust(3)+"%]")
        sys.stdout.flush()
        wtch.process(_t,tracks,ref,dut,alignment)
    print
    #print timeit.default_timer()-start_time
    if alignment:
        foutput = "alignment_plots_{0}_run000{1}.root".format(fp.sensor_name,fp.run_number)
        wtch.store_alignment(foutput)
        return wtch
    else:
        print wtch
        foutput = "sensor_maps_{0}_run000{1}.root".format(fp.sensor_name,fp.run_number)
        print "File created at '{0}'".format(foutput)
        wtch.store_histos(foutput)
