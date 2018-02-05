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
# 5. Be sure the binning is taking into account the strip resolution (at least bin_width = 1/4*pitch or so) ??

DEBUG=False

# Be sure use the same alignment files: ct
ALIGN_FILE_FORMAT= "aligncts_plane_{0}_{1}.txt"

# -- Take care of the activated branches when reading 
#    the root file
class branch_list(object):
    """Helper class containing all the active 
    branches list. To be used as a global and 
    unique instance
    """
    def __init__(self):
        """
        """
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
        for br in brlist:
            self._active_br.append(br)
            tree.SetBranchStatus(br,1)
        # Deactivate all but the active
        for b in filter(lambda x: not x.GetName() in self._active_br, \
                tree.GetListOfBranches()):
            tree.SetBranchStatus(b.GetName(),0)
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
    from math import pi
    # Dictionaries present into a alignment file
    keywords = [ 'iteration', 'x_offset', 'y_offset', 'rot', 'tilt', 'turn', 'dz' ]
    units    = { 'x_offset': (1.0,"mm"), 'y_offset':(1.0,'mm'), \
            'rot': (180.0/pi,'deg.'), 'tilt': (180.0/pi,'deg.'),\
            'turn':(180.0/pi,'deg.'), 'dz': (1.0,'mm') }
    # hyperplane definition : 
    #   x_mis-x_old = d_x - slope_x*d_z - y*slope_x*d_tilt + x*slope_x*d_turn - y*d_rot
    # xoffset, -dz, -titl, +turn, -rot
    hyp_par = [ [0,["x_offset",1.0]], \
                    [1,["dz",-1.0]], \
                    [2,["tilt",-1.0]],\
                    [3,["turn",1.0]],\
                    [4,["rot",-1.0]] ]
    def __init__(self,sensitive_direction):
        """Initialize the data-members
    
        Parameters
        ----------
        sensitive_direction: str
            The direction where the sensor is sensitive. Valid
            words: x|y
        """
        self.sensitive_direction=sensitive_direction
        self.iteration = int(0)
        self.x_offset = 0.0
        self.y_offset = 0.0
        self.turn = 0.0
        self.tilt = 0.0
        self.rot  = 0.0
        self.dz   = 0.0
        if self.sensitive_direction == "y":
            self.hyp_par[0][1][0] = "y_offset"

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
        iteration = self.iteration
        for attr in self.__dict__.keys():
            # Take into account the sign
            try:
                factor = filter(lambda (i,(name,fact)): name == attr,self.hyp_par)[0][0]
            except IndexError:
                # the other offset 
                factor = 1.0
            prov = getattr(self,attr)+getattr(other,attr)
            setattr(self,attr,prov)
        # next iterations, just overwritte the wrong value
        # when added from the other
        self.iteration = iteration+1
        # Also, not adding up the sensitive direction
        self.sensitive_direction = other.sensitive_direction
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
        fobj: file object
            The alignment file 

        Return
        ------
        A object class containing the relevant alignment
        attributes

        Raises
        ------
        IOError
            The file is not present

        KeyError
            The run number is not available in the dictionary
        """
        lines=fobj.readlines()
        for i in lines:
            try: 
                align_word = filter(lambda x: i.find(x) == 0,alignment_cts.keywords)[0]
            except IndexError:
                continue
            setattr(self,align_word,float(i.split("=")[-1]))

    def update(self,filename,run_number):
        """Write down the current snapshoot of the data-members
        into a file
        """
        import time
        import datetime

        # Dump back to the file
        file_content = '#!/usr/bin/env python\n'
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        file_content+= '"""Automaticaly updated by `alignment_cts` class {0}\n"""\n'.format(st)
        for key in alignment_cts.keywords:
            file_content+= "{0} = {1}\n".format(key,getattr(self,key))
        with open(filename,"w") as f:
            f.write(file_content)

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
    
def matmult(a,b):
    """PROVISIONAL
    Just copied from https://stackoverflow.com/questions/10508021/matrix-multiplication-in-python
    """
    zip_b = zip(*b)
    return [[sum(ele_a*ele_b for ele_a,ele_b in zip(row_a,col_b))\
            for col_b in zip_b] for row_a in a]

def matvectmult(a,b):
    """PROVISIONAL
    a matrix
    b vector assuming same number of elements than the matrix row length
    """
    return tuple([sum(el_a*b[j] for (j,el_a) in enumerate(row_a)) for row_a in a])

class hits_plane_accessor(object):
    """Access to the TBranches contained in a 
    TTree created by the EUTelTreeCreator class from
    the https://github.com/duartej/eutelescope/ package.
    
    Alignment notes
    ---------------
    When creating an instance, the code will look 
    for a (pre-formatted name) file containing
    the alignment variables. If found, the alignment
    constants are initialized by copying the content
    of that file. If not, a pre-alignment is performed
    just in the x-offset, by subtracting the difference
    between the local and global coordinates

    See the function tracks_accessor.get_point_in_sensor_frame
    to see how the alignment constants are used
    
    Attributes
    ----------
    id: int
        The telescope plane ID
    sensor_name: str
        The name of the sensor
    orientation: ntuple(float)
        The factors to be multiplied when 
        moving a telescope hit into the sensor coordinate system
        (see tracks_accessor.get_point_in_sensor_frame method)
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
    x_strip: ROOT.std.vector(float)()
        The vector of x local position of the hits 
        in channel number
    y_strip: ROOT.std.vector(float)()
        The vector of y local position of the hits
        in channel number
    x_channels: int
        The number of channels in the x-direction
    y_channels: int
        The number of channels in the y-direction
    sensitive_direction: str
        The name of the sensitive direction,
    sC: ROOT.std.vector(float)()
        The position in the sensitive direction
    sC_local: ROOT.std.vector(float)()
        The local position in the sensitive direction
    sC_channels: int
        The number of channels in the sensitive direction
    other_channels: int
        The number of channels in the non-sensitive direction
    charge: ROOT.std.vector(float)()
        The hit total charge in ADCs counts 
    n_cluster: ROOT.std.vector(int)()
        The total number of elements forming the cluster
    eta: ROOT.std.vector(float)()
        The eta value (charge distribution on the cluster)
    sizeX: float
        The sensor dimension in X
    sizeY: float
        The sensor dimension in Y
    pitchX: float
        The pitch in the x-direction
    pitchY: float
        The pitch in the y-direction
    orientation: (float,float,float)
        The orientation of the sensor plane w.r.t. to the
        local coordinate system
    run_number: int
        The run number of the accessed data
    track_link: dict((int,int)
        The index of the track which is related with the hit.
        A match condition between the track and the hit must 
        be fulfilled
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
    
    def __init__(self,tree,planeid,sensor_name):
        """Accessors for the branches related with Alibava
        sensors hits

        Parameters
        ----------
        tree: ROOT.TTree
            The ntuple 
        planeid: int
            The sensor plane ID
        sensor_name: str
            The name of the sensor
        """
        import ROOT
        import array
        from .SPS2017TB_metadata import get_orientation
        from .SPS2017TB_metadata import sensor_name_spec_map as specs

        self.id = planeid
        self.sensor_name = sensor_name
        # ----------------------------------
        global ACTIVE_BRS
        # -- De-activate those branches not needed
        pre_brnames = ["hit_X","hit_Y","hit_Z","hit_XLocal","hit_YLocal",\
                "hit_total_charge","hit_Ncluster","hit_cluster_eta"]
        brnames = map(lambda x: "{0}_{1}".format(x,self.id),pre_brnames)
        ACTIVE_BRS.update(tree,brnames+["RunNumber"])
        # ----------------------------------
        # Check if there is any problem
        if len(filter(lambda br: br.GetName() == "hit_X_{0}".format(self.id) ,tree.GetListOfBranches())) == 0:
            raise RuntimeError("The INDEX '{0}' does not identify any Alibava"\
                    " sensor plane".format(self.id))
        # Activate the tree if wasn't

        #self.run_number = tree.RunNumber
        self.x = ROOT.std.vector(float)()
        tree.SetBranchAddress("hit_X_{0}".format(self.id),self.x)
        self.y = ROOT.std.vector(float)()
        tree.SetBranchAddress("hit_Y_{0}".format(self.id),self.y)
        self.x_local = ROOT.std.vector(float)()
        tree.SetBranchAddress("hit_XLocal_{0}".format(self.id),self.x_local)
        self.y_local = ROOT.std.vector(float)()
        tree.SetBranchAddress("hit_YLocal_{0}".format(self.id),self.y_local)
        self.charge = ROOT.std.vector(float)()
        tree.SetBranchAddress("hit_total_charge_{0}".format(self.id),self.charge)
        self.n_cluster = ROOT.std.vector(int)()
        tree.SetBranchAddress("hit_Ncluster_{0}".format(self.id),self.n_cluster)
        self.eta = ROOT.std.vector(float)()
        tree.SetBranchAddress("hit_cluster_eta_{0}".format(self.id),self.eta)

        self.run_number = 0.0 
        tree.SetBranchAddress("RunNumber",self._run_number)

        # Get the maximum and minimum in x and y (to define the fiducial cuts)
        ROOT.gROOT.SetBatch()
        prv_c = ROOT.TCanvas()
        # Auto-obtain the sensitive direction (assuming x)
        self._sensitive_direction = "x"
        dummy = tree.Draw("hit_X_{0}>>hdum(100)".format(self.id))
        hdummy = ROOT.gDirectory.Get("hdum")
        # If everything is just in one bin -> No sensitive
        tot_entries = hdummy.GetEntries()
        for i in xrange(1,hdummy.GetNbinsX()+1):
            if tot_entries == hdummy.GetBinContent(i):
                # -- this is not the sensitive axis,
                #    found a bin with all the entries
                dummy = tree.Draw("hit_Y_{0}>>hdumY(100)".format(self.id))
                hdummy = ROOT.gDirectory.Get("hdumY")
                self._sensitive_direction = "y"
                break
        self.rmin = hdummy.GetXaxis().GetBinLowEdge(1)
        self.rmax = hdummy.GetXaxis().GetBinUpEdge(100)
        
        # -- For the modulo definition, let's be sure we don't use any masked region
        #    or channels
        del prv_c

        # Initialize the link index 
        self.track_link = {}

        # Let's assume that z doesn't change between events, as 
        # it is a fixed constrain, so find the first entry with
        # a valid value and then obtain
        #branch_z = getattr(tree,"hit_Z_{0}".format(self.id))
        branch_z = ROOT.std.vector(float)()
        tree.SetBranchAddress("hit_Z_{0}".format(self.id),branch_z)
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
        # De-activate the z-branch, not needed anymore
        tree.SetBranchStatus("hit_Z_{0}".format(self.id),0)

        # Initialize resync static variable
        if not hits_plane_accessor.resync.has_key(self.id):
            hits_plane_accessor.resync[self.id] = True
        # Check the presence of an alignment file, and then update the alignment
        # constants. Assuming in the telescope frame system
        if hits_plane_accessor.resync[self.id]:
            print "ALIGNMENT INPUT [PLANE:{0}]".format(self.id)
            # Just do it if there weren't initialized previously
            if not hits_plane_accessor.align_constants.has_key(self.id):
                hits_plane_accessor.align_constants[self.id] = alignment_cts(self.sensitive_direction)
            filename = ALIGN_FILE_FORMAT.format(self.id,self.run_number)
            try:
                with open(filename) as f:
                    align = alignment_cts(self.sensitive_direction)
                    align.parse_alignment(f)
                    # Update the alignment constants
                    hits_plane_accessor.align_constants[self.id] = align
            except IOError:
                # First access, propagate the prealignment from the Marlin 
                # processor here:
                _k = 1
                # Just to be sure we have data (if x_local is not empty, 
                # then y_local is not empty too) 
                while self.x_local.size() == 0:
                    dummy = tree.GetEntry(_k)
                    _k+=1
                if self.sensitive_direction == "x":
                    hits_plane_accessor.align_constants[self.id].x_offset = self.x_local[0]-self.x[0]
                elif self._sensitive_direction == "y":
                    hits_plane_accessor.align_constants[self.id].y_offset = self.y_local[0]-self.y[0]
            # Not re-synchronized until new order
            print "Inferred sensitive direction: {0}".format(self.sensitive_direction.upper())
            print hits_plane_accessor.align_constants[self.id]
            hits_plane_accessor.resync[self.id] = False
        
        # Dimensions of the sensor (note some issues with LGADS and iLGADs)
        self.sizeX = specs[self.sensor_name].sizeX 
        self.sizeY = specs[self.sensor_name].sizeY 
        self.pitchX = specs[self.sensor_name].pitchX 
        self.pitchY = specs[self.sensor_name].pitchY 
        if self.sensor_name.lower().find("ilgad") == 0:
            self.sizeY = 7.2*MM
            self.sizeX = 7.2*MM
            self.pitchX = self.sizeX
        elif self.sensor_name.lower().find("lgad") == 0:
            self.sizeY = 5.2*MM
            self.sizeX = 5.2*MM
            self.pitchX = self.sizeX
        self.x_channels = self.sizeX/self.pitchX
        self.y_channels = self.sizeY/self.pitchY

        # --The position in channel number
        self.x_strip = lambda i: (self.x_local[i]+0.5*self.sizeX)/self.pitchX+0.5
        self.y_strip = lambda i: (self.y_local[i]+0.5*self.sizeY)/self.pitchY+0.5

        # Defining some useful accessors which depends on the sensitive direction
        if self.sensitive_direction == "x":         
            self.sC = self.x
            self.sC_local = self.x_local
            self.resolution = self.pitchX
            self.strip_position = self.x_strip
            self.sC_channels = self.x_channels
            self.other_channels = self.y_channels
            self.get_sC_channel = self.get_x_channel
        elif self.sensitive_direction == "y":
            self.sC = self.y
            self.sC_local = self.y_local
            self.resolution = self.pitchY
            self.strip_position = self.y_strip
            self.sC_channels = self.y_channels
            self.other_channels = self.x_channels
            self.get_sC_channel = self.get_y_channel

        # The maximum distance to consider a hit matched a track
        self.matching_distance = 2.0*self.resolution
        if self.sensor_name.find("REF") == -1:
            # Loosing the cut
            self.matching_distance = 6.0*self.resolution
        
        # ----------------------------------
        #  Resolution

        # Last, define the orientation of this sensor, 
        # i.e. the factors to multiply each component when
        # moving a telescope hit into the sensor coordinate system
        # (see tracks_accessor.get_point_in_sensor_frame method)
        # |-------------------------------------------------------------------
        # As the sensor reference system is defined (wrt to the telescope 
        # plane) as  e0=(-1,0,0), e1=(0,1,0), n=(0,0,-1)pe gear file )
        # Therefore, for a sensor with the orientation -1  0 
        #                                               0 -1
        # the planes are like:     / z
        #                         /
        #                        ----> x
        #                        |
        #                        | y
        # 
        # The sensor planes are e0 = -x, e1 = y  and n = -z
        # In the case than the sensor have an orientation 1 0 
        #                                                 0 1
        # the planes are like:  y|  / z
        #                        | /
        #                  x <----/
        # 
        # and the sensor planes are e0 = x, e1 = -y and n = -z
        rotation1 = int(get_orientation(self.sensor_name,self.run_number))
        if rotation1 == 1:
            self.orientation = (1.0,-1.0,-1.0)
        else:
            self.orientation = (-1.0,1.0,-1.0)

    @property
    def n(self):
        """The number of hit elements

        Return
        ------
        int
        """
        return self.x.size()

    @property
    def sensitive_direction(self):
        """The sensitive direction of the sensor

        Return
        ------
        str: x|y
        """
        return self._sensitive_direction
    
    @property
    def run_number(self):
        """The run number

        Return
        ------
        int
        """
        return self._run_number[0]
    @run_number.setter
    def run_number(self,value):
        """The setter for the run number
        """
        import array
        self._run_number = array.array('i',[int(value)])

    def get_x_channel(self,x):
        """The equivalent channel number for the
        x-position

        Parameters
        ----------
        x: float
            The position in the local sensor frame

        Return
        ------
        The x position in channel number
        """
        return (x+0.5*self.sizeX)/self.pitchX-0.5
    
    def get_y_channel(self,y):
        """The equivalent channel number for the
        y-position 

        Parameters
        ----------
        y: float
            The position in the local sensor frame

        Return
        ------
        The y position in channel number
        """
        return (y+0.5*self.sizeY)/self.pitchY-0.5

    def is_within_fiducial(self,x,y):
        """Whether or not the particle is inside the defined
        fiducial region of the detector

        Parameters
        ----------
        x: float
            The x-component of the point, in sensor (corrected)
            components reference system
        y: float
            The y-component of the point, in sensor (corrected)
            components reference system
        """
        # Note that is very important that x and y are the described
        # in the sensor reference frame, and the origin of the 
        # reference system is in the center of the detector
        if abs(x) > self.sizeX/2.:
            return False
        if abs(y) > self.sizeY/2.:
            return False
        return True
    
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
    
    ### -- shoould --> @memoize
    # Maybe not needed because is only call once (as the get_point_in_sensor_frame
    # is memoize and is the only client..
    def get_normal_vector(self):
        """Get the normal vector to the sensor plane, if perfectly 
        aligned it must be (0,0,-1) (note that the plane is defined
        following the opposite direction of the beam)

        Return
        ------
        vector normal: (float,float,float)
        """
        if not hits_plane_accessor.normal.has_key(self.id):
            # Calculate it and memorize
            hits_plane_accessor.normal[self.id] = (self.cos_tilt(self.id)*self.sin_turn(self.id),\
                    -self.sin_tilt(self.id),\
                    -self.cos_tilt(self.id)*self.cos_turn(self.id))
            # Do calculations with matricially and compare
            #print self.id,hits_plane_accessor.normal[self.id]
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
    
    @staticmethod
    @memoize
    def matrix_turn(pid):
        """Memorize turn matrix angle in the telescope plane

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        [ [float,float,float],
          [float,float,float],
          [float,float,float]]
        """
        return  [ [hits_plane_accessor.cos_turn(pid), 0.0, hits_plane_accessor.sin_turn(pid)],\
                  [0.0,1.0,0.0],
                  [-hits_plane_accessor.sin_turn(pid), 0.0, hits_plane_accessor.cos_turn(pid)]]

    @staticmethod
    @memoize
    def matrix_tilt(pid):
        """Memorize tilt matrix angle in the telescope plane

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        [ [float,float,float],
          [float,float,float],
          [float,float,float]]
        """
        return  [ [1.0, 0.0, 0.0],\
                  [0.0, hits_plane_accessor.cos_tilt(pid),-hits_plane_accessor.sin_tilt(pid)],\
                  [0.0,hits_plane_accessor.sin_tilt(pid), hits_plane_accessor.cos_tilt(pid)] ]

    @staticmethod
    @memoize
    def matrix_rot(pid):
        """Memorize rot matrix angle in the telescope plane

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        [ [float,float,float],
          [float,float,float],
          [float,float,float]]
        """
        return  [ [hits_plane_accessor.cos_rot(pid), -hits_plane_accessor.sin_rot(pid),0.0],\
                  [hits_plane_accessor.sin_rot(pid), hits_plane_accessor.cos_rot(pid),0.0 ],\
                  [0.0,0.0,1.0] ]

    @staticmethod
    @memoize
    def matrix_turn_rec(pid):
        """Memorize turn matrix angle in the telescope plane

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        [ [float,float,float],
          [float,float,float],
          [float,float,float]]
        """
        return  [ [hits_plane_accessor.cos_turn(pid), 0.0, -hits_plane_accessor.sin_turn(pid)],\
                  [0.0,1.0,0.0],
                  [hits_plane_accessor.sin_turn(pid), 0.0, hits_plane_accessor.cos_turn(pid)]]

    @staticmethod
    @memoize
    def matrix_tilt_rec(pid):
        """Memorize tilt matrix angle in the telescope plane

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        [ [float,float,float],
          [float,float,float],
          [float,float,float]]
        """
        return  [ [1.0, 0.0, 0.0],\
                  [0.0,hits_plane_accessor.cos_tilt(pid),hits_plane_accessor.sin_tilt(pid)],\
                  [0.0,-hits_plane_accessor.sin_tilt(pid),hits_plane_accessor.cos_tilt(pid)] ]

    @staticmethod
    @memoize
    def matrix_rot_rec(pid):
        """Memorize rot matrix angle in the telescope plane

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        [ [float,float,float],
          [float,float,float],
          [float,float,float]]
        """
        return  [ [hits_plane_accessor.cos_rot(pid), hits_plane_accessor.sin_rot(pid),0.0],\
                  [-hits_plane_accessor.sin_rot(pid), hits_plane_accessor.cos_rot(pid),0.0 ],\
                  [0.0,0.0,1.0] ]

    @staticmethod
    @memoize
    def matrix_sensor(pid):
        """Memorize transformation matrix to the sensor reference frame

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        [ [float,float,float],
          [float,float,float],
          [float,float,float]]
        """
        #import numpy as np
        # np.matrix(hits_plane_accessor.matrix_turn(pid))*\
        #   np.matrix(hits_plane_accessor.matrix_tilt(pid))*\
        #   np.matrix(hits_plane_accessor.matrix_rot(pid))
        pr = matmult(hits_plane_accessor.matrix_tilt(pid),hits_plane_accessor.matrix_rot(pid))
        return matmult(hits_plane_accessor.matrix_turn(pid),pr)
    
    @staticmethod
    @memoize
    def matrix_sensor_rec(pid):
        """Memorize reciprocal transformation matrix to the sensor reference frame

        Parameters
        ----------
        pid: int
            The plane id

        Returns
        -------
        [ [float,float,float],
          [float,float,float],
          [float,float,float]]
        """
        #import numpy as np
        # np.matrix(hits_plane_accessor.matrix_turn(pid))*\
        #   np.matrix(hits_plane_accessor.matrix_tilt(pid))*\
        #   np.matrix(hits_plane_accessor.matrix_rot(pid))
        pr = matmult(hits_plane_accessor.matrix_tilt_rec(pid),hits_plane_accessor.matrix_turn_rec(pid))
        return matmult(hits_plane_accessor.matrix_rot_rec(pid),pr)
    
    def equal(self,i,other,j):
        """Allow comparation beetween hits, we assume equal
        if (x,z) are equal

        Parameters
        ----------
        i: int
            The hit index of this (self) instance to be 
            compare
        other: hits_plane_accessor instance (or any other instance
            with a method 'z.__getitem__')
            The other instance to compare with
        j: int
            The hit index of the 'other' instance

        Return
        ------
        bool
        """
        return abs(self.sC[i]-other.sC[j]) < 1e-9 \
                and abs(self.z[i]-other.z[j]) < 1e-9

    def close_hit(self,point):
        """Find the closest hit to the given point. 
        If there is no hits in the event, it returns -1.
        
        Parameters
        ----------
        point: (float,float,float)
            The point position where to search the closest hit 
        """
        # Loop over all available hits and found the closest
        ind_distance_list = map(lambda (i,r): (i,r-point[1]),enumerate(self.sC_local))
        if len(ind_distance_list) == 0 :
            return (-1,-9999)
        return min(ind_distance_list,key=lambda (i,d): abs(d))

    def is_within_matching_distance(self,dist_check):
        """Whether or not the distance introduced is within
        the considered matching distance for that sensor

        Parameters
        ----------
        dist_check: float
            The distance to check

        Return
        ------
        bool
        """
        return abs(dist_check) < self.matching_distance

    def fill_distance_histo(self,rtele,h):
        """Fill the histogram of the distance of the hit to the 
        rtele position (assumed to be given in the sensor plane 
        coordinate) for all hits in the event

        Parameter
        ---------
        rtele: float
            The (sensitiv direction) position in the sensor plane
        h: ROOT.TH1F
        """
        dummy = map(lambda r: h.Fill(r-rtele),self.sC_local)
    
    def fill_correlation_histo(self,rtele,h):
        """Fill the histogram of position of the hit versus the
        rtele position (assumed to be given in the sensor plane 
        coordinate) for all hits in the event

        Parameter
        ---------
        rtele: float
            The (sensitiv direction) position in the sensor plane
        h: ROOT.TH"F
        """
        dummy = map(lambda r: h.Fill(r,rtele),self.sC_local)

    def match_isolated(self,track_acc,histos,is_ref=False,res=None):
        """Search for tracks matching (within a resolution) the list 
        of hits in the current event.

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
              (see last histo) vs. the x-prediction. See `dy_x_h` in the 
              processor class. This histo is used for alignment (rotation)            
            - A ROOT.TProfile to store the variation of the residual in y
              (see last histo) vs. the y-prediction. See `dy_y_h` in the 
              processor class. This histo is used for alignment (tilt)
            - A ROOT.TProfile to store the variation of the residual in x
              (see last histo) vs. the x-prediction. See `dx_xtx_h` in the 
              processor class. This histo is used for alignment (turn)
            - A ROOT.TProfile to store the variation of the residual in y
              (see last histo) vs. the y-slope. See `dx_tx_h` in the 
              processor class. This histo is used for alignment (z-shift)
            - XXX MISSING HISTOS DESCRIPTION
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
        import array

        if not res:
            res=RESOLUTION[self.id]
        ISOLATION = 0.6*MM
    
        #if track_acc.n != 20 or self.n != 1: 
        #    return []

        # Get the histograms
        hcorr,hres,hdx,hdx_finer,hrot,htilt,hturn,hdz,hplane,\
                hmatch_eff,hy_eff,hntracks_eff,hnisotracks_eff = histos

        # The index of the needed coordinated and some
        # other data members dependent of the sensitive direction
        if self.sensitive_direction == "x":
            # -- the sensitive
            ic = 0
            # -- the complementary
            iccom = 1
            # -- the offset
            r_offset = self.align_constants[self.id].x_offset
            # -- the proper track slope
            track_slope = track_acc.dxdz
            
        elif self.sensitive_direction == "y":
            # -- the sensitive
            ic = 1
            # -- the complementary
            iccom = 0
            # -- the offset
            r_offset = self.align_constants[self.id].y_offset
            # -- the proper track slope
            track_slope = track_acc.dydz

        matched_hits = []
        # non isolated tracks and already used: keep them listed 
        # to avoid re-processing
        used_tracks = []
        # Given the hits in that sensor,
        for ihit in xrange(self.n):
            closest = {}
            # Evaluate the distance from the hit to any track (not used already)
            not_used_track_indices = filter(lambda i: i not in used_tracks,xrange(track_acc.n))
            for i_at_list,itrk in enumerate(not_used_track_indices):
                # ---  XXX --- WHAT??X
                #if int(track_acc.ndof[itrk]) != 6 \
                #        or track_acc.chi2[itrk]/track_acc.ndof[itrk] > 3.0:
                #    used_tracks.append(itrk)
                #    continue
                # Isolation: be sure there is no other track surrounding this one
                (rpred,rtel0) = track_acc.get_point_in_sensor_frame(itrk,self)
                # Before apply assign the alignment, otherwise biasing the distance
                # ---- using all tracks
                hcorr.Fill(self.sC_local[ihit],rpred[ic])
                hdx.Fill((self.sC_local[ihit]-(rpred[ic]-r_offset)))
                # note tha i_at_list catchs the indice at `not_used_track_indices`
                # list, therefore in order to evaluate only those tracks not previously
                # evaluated just use the remaining elements after the current track index
                #non_isolated = filter(lambda (oindex,((o_x,o_y,o_z),_tel)): sqrt((o_x-rpred[0])**2.0+(o_y-rpred[1])**2.0) < ISOLATION,\
                #        map(lambda other_i:  (other_i,track_acc.get_point_in_sensor_frame(other_i,self)),\
                #            not_used_track_indices[i_at_list+1:]))
                non_isolated=[]
                # -- Check if there is another track near-by
                for otrk in xrange(track_acc.n):
                    if otrk == itrk:
                        continue
                    if abs(track_acc.get_point_in_sensor_frame(otrk,self)[0][1]-rpred[1]) < ISOLATION:
                        non_isolated.append(itrk)
                        non_isolated.append(otrk)
                        break
                if len(non_isolated) > 0 :
                    # Another track has been found sourrounding this one, 
                    # not isolated, remove them from the available tracks
                    ##used_tracks.append(itrk)
                    ##dummy = map(lambda (i,other_stuff): used_tracks.append(i), non_isolated)
                    continue
                # -- Store the distance
                closest[abs(self.sC_local[ihit]-rpred[ic])] = itrk
            if len(closest) == 0:
                continue
            #V The hit was matched by at least one track
            #hmatch_eff.Fill(self.sC_local[ihit],1)
            # Get the closest track (in x)
            distance_abs,trk_el =sorted(closest.iteritems())[0]
            #rturn,rtitl,rrot,rsensor = track_acc.get_point_in_sensor_frame(trk_el,self)
            rsensor,rtel = track_acc.get_point_in_sensor_frame(trk_el,self)
            # -- A finer alignment, using only the closest tracks 
            hdx_finer.Fill(self.sC_local[ihit]-rsensor[ic])
            # Note that the prediction is given in the sensor reference frame (z should be zero)
            # Therefore, not to interesting this plot, better
            hplane.Fill(rtel[2]-self.z[0],rtel[0],rtel[1])
            # -- resolution histogram
            dc = self.sC_local[ihit]-rsensor[ic]
            hres.Fill(self.sC_local[ihit],dc)
            # Fill only those passing the cuts 
            if distance_abs < res:
                # Fill the matched hits
                matched_hits.append( (ihit,trk_el) )
                # Tag the cut as used
                used_tracks.append(trk_el)
                # and link the track with the hit
                self.track_link[trk_el] = ihit
                ## Fill the alignment histograms after being accepted
                #  as providing from the hit
                ## -- rot
                hrot.Fill(rsensor[iccom],dc)
                ## -- tilt
                #htilt.Fill(rsensor[iccom]*track_slope[trk_el]*UM,dc)
                htilt.Fill(rsensor[ic],dc)
                # -- turn
                hturn.Fill(rsensor[ic]*track_slope[trk_el]*UM,dc)
                ## -- dz 
                hdz.Fill(track_slope[trk_el]*UM,dc)
            # The hit was matched by at least one track
            hmatch_eff.Fill(self.sC_local[ihit],int(distance_abs<res))
            hy_eff.Fill(rsensor[iccom],int(distance_abs<res))
            hntracks_eff.Fill(track_acc.n,int(distance_abs<res))
            # Number of non-isolated tracs: 
            #    number of used tracks - number of matched hits (hit -> used track)
            hnisotracks_eff.Fill(len(used_tracks)-len(matched_hits),int(distance_abs)<res)
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
        import ROOT

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
        self.x = ROOT.std.vector(float)()
        tree.SetBranchAddress("{0}_X_{1}".format(prefix,self.id),self.x)
        self.y = ROOT.std.vector(float)()
        tree.SetBranchAddress("{0}_Y_{1}".format(prefix,self.id),self.y)
        self.x_local = ROOT.std.vector(float)()
        tree.SetBranchAddress("{0}_XLocal_{1}".format(prefix,self.id),self.x_local)
        self.y_local = ROOT.std.vector(float)()
        tree.SetBranchAddress("{0}_YLocal_{1}".format(prefix,self.id),self.y_local)
        self.z       = ROOT.std.vector(float)()
        tree.SetBranchAddress("{0}_Z_{1}".format(prefix,self.id),self.z)
        self.track_index = ROOT.std.vector(int)()
        tree.SetBranchAddress("{0}_index_{1}".format(prefix,self.id),self.track_index)
    
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
        other: track_hits_plane_accessor instance (or any other instance
            with a method 'z.__getitem__')
            The other instance to compare with
        j: int
            The hit index of the 'other' instance

        Return
        ------
        bool
        """
        return abs(self.x[i]-other.x[j]) < 1e-9 \
                and abs(self.y[i]-other.x[j]) < 1e-9 \
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
    telescope_measured_hits: list(track_hits_plane_accessor)
        The telescope (measured) hits organized by planes
    telescope_fitted_hits: list(track_hits_plane_accessor)
        The telescope (fitted) hits organized by planes
    dut_measured_hits: track_hits_plane_accessor
        The DUT (measured hits), note that this data member
        is only usable if the tracks were fitted in the 
        Marlin processor using the DUT sensor as telescope
        plane
    dut_fitted_hits: track_hits_plane_accessor
        The DUT (fitted hits), note that this data member
        is only usable if the tracks were fitted in the 
        Marlin processor using the DUT sensor as telescope
        plane
    ref_measured_hits: track_hits_plane_accessor
        The REF (measured hits), note that this data member
        is only usable if the tracks were fitted in the 
        Marlin processor using the DUT sensor as telescope
        plane
    ref_fitted_hits: track_hits_plane_accessor
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
        import ROOT
        global ACTIVE_BRS
        # -- De-activate those branches not needed
        brnames = ["trk_hit_fit_X_{0}".format(dut_plane),"trk_hit_fit_X_{0}".format(ref_plane),\
                "trk_refPoint_X","trk_refPoint_Y","trk_refPoint_Z","trk_dxdz","trk_dydz",\
                "trk_Ndf","trk_chi2"]
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
        self.x0 = ROOT.std.vector(float)()
        tree.SetBranchAddress("trk_refPoint_X",self.x0)
        self.y0 = ROOT.std.vector(float)()
        tree.SetBranchAddress("trk_refPoint_Y",self.y0)
        self.z0 = ROOT.std.vector(float)()
        tree.SetBranchAddress("trk_refPoint_Z",self.z0)
        # Director vector
        self.dxdz = ROOT.std.vector(float)()
        tree.SetBranchAddress("trk_dxdz",self.dxdz)
        self.dydz = ROOT.std.vector(float)()
        tree.SetBranchAddress("trk_dydz",self.dydz)
        # -- quality
        self.chi2 = ROOT.std.vector(float)()
        tree.SetBranchAddress("trk_chi2",self.chi2)
        self.ndof = ROOT.std.vector(float)()
        tree.SetBranchAddress("trk_Ndf",self.ndof)
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

        # -- The isolation condition
        self.isolation_condition = 0.6*MM

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
        # -- dict to be used to list of the non-isolated tracks
        self._cache_nonisolated = {}
    
    def is_isolated(self,i,hit_z_position):
        """Whether or not the track is isolated, i.e.
        if there is no other track in the plane defined by the 
        z position inside the isolation cone

        Parameters
        ----------
        i: int
            The index of the track
        hit_z_position: hit_accessor
            The hit accessor which defines the plane

        Return
        ------
        bool: False if there is any track inside the isolation
              cone
        """
        from math import sqrt
        
        # -- It was found previously..
        if self._cache_nonisolated.has_key(i):
            return False
        # -- Get the track at the plane
        ((xpred,ypred,zpred),rtel0) = self.get_point_in_sensor_frame(i,hit_z_position)
        # -- Get the list of non isolated (CONE OR SQUARES@??) XXX
        non_isolated = filter(lambda (oindex,((o_x,o_y,o_z),_tel)): sqrt((o_x-xpred)**2.0+(o_y-ypred)**2.0) < self.isolation_condition,\
                map(lambda other_i:  (other_i,self.get_point_in_sensor_frame(other_i,hit_z_position)),xrange(self.n)))
        if len(non_isolated) == 0:
            return True
        # -- if not isolated, update some stuff
        # -- Update the current event cache
        self._cache_nonisolated[i] = 1
        # -- And the others (reciprocal quality)
        for noniso_i,stuff in non_isolated:
            self._cache_nonisolated[noniso_i] = 1
        return False
    
    def fill_isolation_histograms(self,itrk,hitobj,h):
        """Fill the isolation histogram. The distance between
        all track pairs in the sensor plane are calculated and
        the corresponding histogram, filled

        Parameters
        ----------
        itrk: int
            Track index
        hitobj: hits_plane_accessor
            The sensor accesor
        h: ROOT.TH1F
            The isolation histogram
        """
        from math import sqrt
        
        ((xpred,ypred,zpred),rtel0) = self.get_point_in_sensor_frame(itrk,hitobj)
        dummy = map(lambda ((o_x,o_y,o_z),_tel): h.Fill(sqrt((o_x-xpred)**2.0+(o_y-ypred)**2.0)),\
                map(lambda other_i:  self.get_point_in_sensor_frame(other_i,hitobj),xrange(itrk+1,self.n)))

    def associate_hits(self,refhits,duthits,histos):
        """Perform the hits association to each track found,
        given that the tracks can only be associated if they
        are isolated (see isolation_condition datamember).
        The associating criteria is defined through the 
        hits accessors themselves.

        Parameters
        ----------
        refhits: hits_plane_accessor
            The accessor to the REF hits of the event
        duthits: hits_plane_accessor
            The accessor to the DUT hits of the event
        histos: dict(list(TObject))
            A list of histograms, profiles and maps per
            sensor, to fill during the association


        Return
        ------
        tracks_d: dict( (int,(int,int) )
            A dict with the indices of the tracks associated
            to the indices of the matched hits. If there is 
            not, a -1 is placed instead of the index.
            { track_index: (REF-hit-index,DUT-hit-index) ,... }
        """
        from math import sqrt

        # Get the histograms
        hcorr,hdx,hdx_finer,hdx_finer_wide,hplane= histos[refhits.id]
        hcorr_d,hdx_d,hdx_finer_d,hdx_finer_wide_d,hplane_d= histos[duthits.id]
        
        if refhits.sensitive_direction == "x":
            # -- the sensitive
            ic = 0
            # -- the complementary
            iccom = 1
            # -- the offset
            r_offset_r = refhits.align_constants[refhits.id].x_offset
            r_offset_d = refhits.align_constants[duthits.id].x_offset
            # -- the proper track slope
            #track_slope = self.dxdz
        elif refhits.sensitive_direction == "y":
            # -- the sensitive
            ic = 1
            # -- the complementary
            iccom = 0
            # -- the offset
            r_offset_r = refhits.align_constants[refhits.id].y_offset
            r_offset_d = refhits.align_constants[duthits.id].y_offset
            # -- the proper track slope
            #track_slope = self.dydz

        # Start track loop
        tracks_d = {}
        for itrk in xrange(self.n):
            (r_ref,t_ref) = self.get_point_in_sensor_frame(itrk,refhits)
            (r_dut,t_dut) = self.get_point_in_sensor_frame(itrk,duthits)
            # First check isolation (XXX SUBSTITUTE by is_isolated)
            # --- Fill the displacement offset histogram, in order to
            #     to align. First a coarse, raw 
            refhits.fill_distance_histo((r_ref[ic]-r_offset_r),hdx)
            refhits.fill_correlation_histo(r_ref[ic],hcorr)
            duthits.fill_distance_histo((r_dut[ic]-r_offset_d),hdx_d)
            duthits.fill_correlation_histo(r_dut[ic],hcorr_d)
            # --- SLOW ALGORITHM -- TRY TO IMPROVE IT
            is_isolated = True
            for otrk in filter(lambda i: i  != itrk,xrange(self.n)):
                (o_ref,o_ref) = self.get_point_in_sensor_frame(otrk,refhits)
                # Just asking for isolation in REF (XXX what about DUT?)
                if abs(r_ref[0]-o_ref[0]) < self.isolation_condition \
                        and abs(r_ref[1]-o_ref[1]) < self.isolation_condition:
                    is_isolated = False
                    break
            if not is_isolated:
                continue
            # -- Now find the closest hit to that track (-1 is returned if any)
            iref,distance_ref = refhits.close_hit(r_ref)
            idut,distance_dut = duthits.close_hit(r_dut)
            # -- Some plots: smoother residual evaluation (and alignment) 
            #    using closest hits to a track (event if the hit is re-used)
            if iref != -1:
                hdx_finer.Fill(refhits.sC_local[iref]-r_ref[ic])
                hdx_finer_wide.Fill(refhits.sC_local[iref]-r_ref[ic])
                hplane.Fill(t_ref[2]-refhits.z[0],t_ref[0],t_ref[1])
            if idut != -1:
                hdx_finer_d.Fill(duthits.sC_local[idut]-r_dut[ic])
                hdx_finer_wide_d.Fill(duthits.sC_local[idut]-r_dut[ic])
                hplane_d.Fill(t_dut[2]-duthits.z[0],t_dut[0],t_dut[1])
            # -- Association 
            if not refhits.is_within_matching_distance(distance_ref):
                iref = -1
            else:
                # equivalently, link the track with the hit
                refhits.track_link[itrk] = iref
            if not duthits.is_within_matching_distance(distance_dut):
                idut = -1
            else:
                # equivalently, link the track with the hit
                duthits.track_link[itrk] = idut
            # Fill the track dictionary with matching information
            tracks_d[itrk] = (iref,idut)
        return tracks_d

    
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
                    # Get the index of the hits_plane_accessor instance
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
        from math import sqrt

        # |--- Telescope coordinates -----------------------------------------|
        # Normal vector of the plane (defined in opposite direction to the 
        # beam). With perfect alignment would be: (0,0,-1)
        n = hit.get_normal_vector()
        #  - thea sensor position: (including alignment corrections on z)
        rpp = (0.0,0.0,hit.z[i]+hit.align_constants[hit.id].dz)
        # Helper vectors: the reference point of the tracks
        r0  = (self.x0[i],self.y0[i],self.z0[i])
        # Helper vectors: the track slopes, or the vector director of the track
        ts = (self.dxdz[i],self.dydz[i],1.0)
        # Track predicted points (given the distance from the center of the 
        # sensor: zc:=z-z0))
        track = lambda _zc: map(lambda (_r,_t): _r+_t*_zc, zip(r0,ts))
        # Let us to redefine zc := zcc-z0 
        # Being a: vector over the sensor pointing from the center to the hit,
        # can be expressed in function of the sensor position vector plus the 
        # track vector at the crossing point zc --> a = track-rpp, as a is on
        # the plane, is perpendicular to n, therefore solving:  
        #   (track(zc)-rpp) \cdot n = 0 
        # the crossing z point is obtained
        zc = -sum(map(lambda i: n[i]*(r0[i]-rpp[i]),xrange(len(n))))/sum(map(lambda i: n[i]*ts[i],xrange(len(n))))
        # The point 
        rpred = track(zc)

        # |--- Sensor plane coordinates --------------------------------------| 
        # As the sensor reference system is defined (wrt to the telescope 
        # plane) as  e0=(-1,0,0), e1=(0,1,0), n=(0,0,-1) if the telescope plane
        # was oriented as the coordinate system.
        # Note, however that the telescope planes are oriented as 
        #             
        #          / z  (beam direction)
        #   x  <--/           
        #         |
        #         | y
        #        
        # This orientation applies after the rotation of the plane
        # by the matrix  -1  0 
        #                 0 -1 
        # (See the rotation1[2,3,4] values of the telescope gear file )
        # Therefore, for a sensor with the same orientation than the telescope
        # planes: e0 = (1,0,0), e1 = (0,-1,0) and n=(0,0,1)
        sign_c = hit.orientation 

        # Getting the hit position vector of the prediction, in the sensor plane
        # (before correcting by the orientation)
        a0 = matvectmult(hit.matrix_sensor_rec(hit.id), \
                map(lambda (i,x): rpred[i]-x,enumerate(rpp)))
        # Correct the telescope prediction by the orientation and displace it 
        # to the position where should it be
        a = map(lambda (i,offset): sign_c[i]*a0[i]+offset, enumerate((hit.align_constants[hit.id].x_offset,\
                hit.align_constants[hit.id].y_offset, 0.0)) )

        return a,rpred 


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
    residual_associated: ROOT.TH2F
        Histogram of the difference between the hit predicted and the
        hit measured, for associated hits
    hcharge_associated: dict(int, ROOT.TH2F)
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

    def __init__(self,minst,run_number=-1):
        """ XXX MISSING DOC
        Parameters
        ----------
        minst: metadata_container instance
        """
        import ROOT

        # Get some useful info to build the histograms
        self.dut_plane = minst.dut_plane
        self.ref_plane = minst.ref_plane

        # The processed run number 
        self.run_number = run_number

        # Asumming strips for the ref (so recall is going 
        # to be used as 2 x pitch)
        ## --- XXX Simulate pixels
        self.pitchY = { minst.dut_plane: minst.dut_pitchY, minst.ref_plane: minst.ref_pitchY }
        ## ----
        self.pitchX = { minst.dut_plane: minst.dut_pitchX, minst.ref_plane: minst.ref_pitchX }
        #self.pitchX = self.pitchY
        # -- Check if is 3-D or strips
        #if minst.dut_name.lower().find("lgad") == -1:
        #    # simulate pixels
        #    self.pitchY[minst.dut_plane] = minst.dut_pitchX
        #else:
        #    # As is going to be use as 2 x pitch
        self.sizeX = { minst.dut_plane: minst.dut_sizeX, minst.ref_plane: minst.ref_sizeX }
        self.sizeY = { minst.dut_plane: minst.dut_sizeY, minst.ref_plane: minst.ref_sizeY }
        # Be careful than iLGAD and LGAD don't have the proper dimensions of the sensors
        if minst.dut_name.lower().find("ilgad") == 0:
            self.sizeY[minst.dut_plane]= 7.2*MM
            self.sizeX[minst.dut_plane]= 7.2*MM
        elif minst.dut_name.lower().find("lgad") == 0:
            self.sizeY[minst.dut_plane]= 5.2*MM
            self.sizeX[minst.dut_plane]= 5.2*MM
        # useful for the histos
        # half size plus some extra to keep alignment factors
        sxdut = self.sizeX[minst.dut_plane]*0.5+0.25
        sydut = self.sizeY[minst.dut_plane]*0.5+0.25
        sxref = self.sizeX[minst.ref_plane]*0.5+0.25
        syref = self.sizeY[minst.ref_plane]*0.5+0.25

        # some info numbers
        self.total_events    = 0
        self.total_events_tracks = 0
        self.events_with_dut = 0
        self.dut_hits        = 0
        self.events_with_ref = 0
        self.ref_hits        = 0
        self.events_with_dut_no_tracks = 0
        self.events_with_ref_no_tracks = 0
        self.events_with_ref_dut_tracks = 0
        self.number_matched_hits = {}
        self.events_with_matched_hits = {}
        # Alingment histos
        # -----------------
        # -- the shift in x
        self.dx_h = { minst.dut_plane: ROOT.TH1F("dx_dut"," ; y_{DUT}-y_{trk}^{pred} [mm]; Entries",200,-sydut,sydut),
                minst.ref_plane: ROOT.TH1F("dx_ref"," ; y_{REF}-y_{trk}^{pred} [mm]; Entries",200,-syref,syref) }
        # -- finer alignment
        self.dx_finer_h = { minst.dut_plane: ROOT.TH1F("dx_finer_dut"," ; x_{DUT}-x_{trk}^{pred} [mm]; Entries",100,-0.3,0.3),
                minst.ref_plane: ROOT.TH1F("dx_finer_ref"," ; x_{REF}-x_{trk}^{pred} [mm]; Entries",100,-0.3,0.3) }
        self.dx_finer_wide_h = { minst.dut_plane: ROOT.TH1F("dx_finer_wide_dut"," ; x_{DUT}-x_{trk}^{pred} [mm]; Entries",200,-1.0,1.0),
                minst.ref_plane: ROOT.TH1F("dx_finer_wide_ref"," ; x_{REF}-x_{trk}^{pred} [mm]; Entries",200,-1.0,1.0) }
        # -- the rot (around z-axis)
        self.dy_x_h = { minst.dut_plane: ROOT.TProfile("dy_x_dut"," ;x_{trk}^{pred} [mm];#Deltay_{DUT} [mm]",\
                        45,-sxdut,sxdut,-0.2,0.2),
                minst.ref_plane: ROOT.TProfile("dy_x_ref"," ;x_{trk}^{pred} [mm];#Deltay_{REF} [mm]",\
                        45,-sxref,sxref,-0.2,0.2) }
        # -- the tilt (around y-axis)
        self.dy_y_h = { minst.dut_plane: ROOT.TProfile("dy_y_dut"," ;y_{trk}^{pred} [mm];#Deltay_{DUT} [mm]",\
                        45,-sxdut,sxdut,-0.2,0.2),
                minst.ref_plane: ROOT.TProfile("dy_y_ref"," ;y_{trk}^{pred} [mm];#Deltay_{REF} [mm]",\
                        45,-sxref,sxref,-0.2,0.2) }
        # -- the turn (around y-axis)
        self.dx_xtx_h = { minst.dut_plane: ROOT.TProfile("dx_xtx_dut"," ;x_{trk}^{pred}*tan(#theta_{x}) [#mum];#Deltax_{DUT} [mm]",\
                        50,-0.1,0.1,-1.2,1.2),
                minst.ref_plane: ROOT.TProfile("dx_xtx_ref"," ;x_{trk}^{pred}*tan(#theta_{x}) [#mum];#Deltax_{REF} [mm]",\
                        50,-0.1,0.1,-1.2,1.2) }
        # -- the delta-z
        self.dx_tx_h = { minst.dut_plane: ROOT.TProfile("dx_tx_dut"," ;tan(#theta_{x}^{trk})#cdot1^{-3};#Deltax_{DUT} [mm]",50,-0.2,0.2,-0.2,0.2),
                minst.ref_plane: ROOT.TProfile("dx_tx_ref"," ;tan(#theta_{x}^{trk})#cdot1^{-3};#Deltax_{REF} [mm]",50,-0.2,0.2,-0.2,0.2) }
        # geometry: z
        #----------
        self.hplane = { minst.dut_plane: ROOT.TH3F("plane_dut",";dz^{pred} [mm];x^{pred} [mm];y^{pred} [mm]",\
                    51,-1.0,1.0,50,-sxdut*1.5,sxdut*1.5,50,-1.5*sydut,1.5*sydut),
                minst.ref_plane: ROOT.TH3F("plane_ref",";dz^{pred} [mm];x^{trk} [mm];y^{pred} Entries",\
                    51,-1.0,1.0,50,-syref*1.5,sxref*1.5,50,-1.5*syref,1.5*syref)}

        self._alignment_histos = self.dx_h.values()+self.dx_finer_h.values()+self.dx_finer_wide_h.values()+self.dy_x_h.values()+\
                self.dx_xtx_h.values()+self.dy_y_h.values()+self.dx_tx_h.values()+\
                self.hplane.values()
        
        # -- The number of channels of this dut
        nch_dut = self.sizeY[self.dut_plane]/minst.dut_pitchY
        nch_ref = self.sizeY[self.ref_plane]/minst.ref_pitchY
        # Associated histos (hits associated to a track)
        # -------------------------------------
        self.residual_associated = { minst.dut_plane: ROOT.TH2F("res_a_dut","y_{DUT} [mm];y_{DUT}-y_{DUT}^{pred} [mm];Entries",\
                    100,-1.1*sydut,1.1*sydut,100,-0.4*MM,0.4*MM), 
                minst.ref_plane: ROOT.TH2F("res_a_ref",";y_{REF} [mm] ;y_{REF}-y_{REF}^{pred} [mm];Entries",\
                    100,-1.1*syref,1.1*syref,100,-0.4*MM,0.4*MM) }
        self.hcharge_associated = { minst.dut_plane: ROOT.TProfile2D("charge_a_dut" ,\
                    ";x_{DUT}^{pred} [mm];y_{DUT}^{pred} [mm];<charge cluster> [ADC]", \
                    300,-1.1*sxdut,1.1*sxdut,300,-1.1*sydut,1.1*sydut),
                minst.ref_plane: ROOT.TProfile2D("charge_a_ref",\
                    ";x_{REF}^{pred} [mm];y_{REF}^{pred} [mm];<charge cluster> [ADC]",\
                    300,-1.1*sxref,1.1*sxref,300,-1.1*syref,1.1*syref) }
        self.hcharge1D_associated = { minst.dut_plane: ROOT.TH1F("charge1D_a_dut" ,\
                    ";charge cluster [ADC];Entries",300,0,600),
                minst.ref_plane: ROOT.TH1F("charge1D_a_ref",";charge cluster [ADC];Entries",300,0,600) }
        self.hclustersize_associated = { minst.dut_plane: ROOT.TProfile2D("clustersize_a_dut" ,\
                    ";x_{DUT}^{pred} [mm];y_{DUT}^{pred} [mm];<cluster size>", \
                    300,-1.1*sxdut,1.1*sxdut,300,-1.1*sydut,1.1*sydut),
                minst.ref_plane: ROOT.TProfile2D("clustersize_a_ref",\
                    ";x_{REF}^{pred} [mm];y_{REF}^{pred} [mm];<cluster size>",\
                    300,-1.1*sxref,1.1*sxref,300,-1.1*syref,1.1*syref) }
        self.hclustersize1D_associated = { minst.dut_plane: ROOT.TH1F("clustersize1D_a_dut" ,\
                    ";cluster size;Entries",10,-0.5,9.5),
                minst.ref_plane: ROOT.TH1F("clustersize1D_a_ref",";cluster size;Entries",10,-0.5,9.5) }
        self.hitmap_associated = { minst.dut_plane: ROOT.TH2F("hitmap_a_dut" ,";x_{DUT}^{pred} [mm];y_{DUT}^{pred} [mm];"\
                    "Entries", 300,-1.1*sxdut,1.1*sxdut,300,-1.1*sydut,1.1*sydut),
                minst.ref_plane: ROOT.TH2F("hitmap_a_ref",";x_{REF}^{pred} [mm];y_{REF}^{pred} [mm];"\
                    "Entries", 300,-1.1*sxref,1.1*sxref,300,-1.1*syref,1.1*syref) }
        # -- Module (2 x pitch)
        self.hclustersize_associated_mod = { minst.dut_plane: ROOT.TProfile2D("clustersize_a_mod_dut" ,\
                    ";mod(x_{trk})_{2pitch} [#mum];mod(y_{trk})_{2pitch}[#mum];<cluster size>", \
                    40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TProfile2D("clustersize_a_mod_ref" ,";mod(x_{trk})_{2pitch} [mm];"\
                    "mod(y_{trk})_{2xpitch} [mm];<cluster size>", \
                    40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }
        self.hcharge_associated_mod = { minst.dut_plane: ROOT.TProfile2D("charge_a_mod_dut",\
                    ";mod(x_{trk})_{2pitch} [#mum]; mod(y_{trk})_{2pitch} [#mum];<cluster charge> [ADC]",\
                    40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TProfile2D("charge_a_mod_ref",\
                    ";mod(x_{trk})_{2pitch} [#mum]; mod(y_{trk})_{2xpitch} [#mum];<cluster charge> [ADC]",\
                    40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }
        self.hitmap_associated_mod = { minst.dut_plane: ROOT.TH2F("hitmap_a_mod_dut",\
                    ";mod(x_{trk})_{2pitch} [#mum]; mod(y_{trk})_{2pitch} [#mum];Entries",\
                    40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TH2F("hitmap_a_mod_ref",\
                    ";mod(x_{trk})_{2pitch} [#mum]; mod(y_{trk})_{2pitch} [#mum];Entries", \
                    40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }
        
        
        associated_histos = self.residual_associated.values()+\
                self.hcharge_associated.values()+self.hcharge1D_associated.values()+\
                self.hitmap_associated.values()+\
                self.hclustersize_associated.values()+self.hclustersize1D_associated.values()+\
                self.hclustersize_associated_mod.values()+self.hcharge_associated_mod.values()+\
                self.hitmap_associated_mod.values()

        # Matched histos (particles (track+REF) associated to a DUT hit
        # -------------------------------------
        # Residuals between REF-DUT (matched tracks family)
        # Residuals, matched family
        self.residual_matched = { minst.dut_plane: ROOT.TH1F("res_m_dut"," ;r_{DUT}-r_{DUT}^{pred} [mm];Entries",\
                        100,-0.4*MM,0.4*MM), 
                minst.ref_plane: ROOT.TH1F("res_m_ref"," ;r_{REF}-r_{REF}^{pred} [mm];Entries",\
                        100,-0.4*MM,0.4*MM) }
        self.hcharge_matched = { minst.dut_plane: ROOT.TProfile2D("charge_m_dut" ,\
                    ";x_{DUT}^{pred} [mm];y_{DUT}^{pred} [mm];<charge cluster> [ADC]", \
                    300,-1.1*sxdut,1.1*sxdut,300,-1.1*sydut,1.1*sydut),
                minst.ref_plane: ROOT.TProfile2D("charge_m_ref",\
                    ";x_{REF}^{pred} [mm];y_{REF}^{pred} [mm];<charge cluster> [ADC]",\
                    300,-1.1*sxref,1.1*sxref,300,-1.1*syref,1.1*syref) }
        self.hcharge1D_matched = { minst.dut_plane: ROOT.TH1F("charge1D_m_dut" ,\
                    ";charge cluster [ADC];Entries",300,0,600),
                minst.ref_plane: ROOT.TH1F("charge1D_m_ref",";charge cluster [ADC];Entries",300,0,600) }
        self.hclustersize_matched = { minst.dut_plane: ROOT.TProfile2D("clustersize_m_dut" ,\
                    ";x_{DUT}^{pred} [mm];y_{DUT}^{pred} [mm];<cluster size>", \
                    300,-1.1*sxdut,1.1*sxdut,300,-1.1*sydut,1.1*sydut),
                minst.ref_plane: ROOT.TProfile2D("clustersize_m_ref",\
                    ";x_{REF}^{pred} [mm];y_{REF}^{pred} [mm];<cluster size>",\
                    300,-1.1*sxref,1.1*sxref,300,-1.1*syref,1.1*syref) }
        self.hclustersize1D_matched = { minst.dut_plane: ROOT.TH1F("clustersize1D_m_dut" ,\
                    ";cluster size;Entries",10,-0.5,9.5),
                minst.ref_plane: ROOT.TH1F("clustersize1D_m_ref",";cluster size;Entries",10,-0.5,9.5) }
        self.hitmap_matched = { minst.dut_plane: ROOT.TH2F("hitmap_m_dut" ,";x_{DUT}^{pred} [mm];y_{DUT}^{pred} [mm];"\
                    "Entries", 300,-1.1*sxdut,1.1*sxdut,300,-1.1*sydut,1.1*sydut),
                minst.ref_plane: ROOT.TH2F("hitmap_m_ref",";x_{REF}^{pred} [mm];y_{REF}^{pred} [mm];"\
                    "Entries", 300,-1.1*sxref,1.1*sxref,300,-1.1*syref,1.1*syref) }
        # -- Module (2 x pitch)
        self.hclustersize_matched_mod = { minst.dut_plane: ROOT.TProfile2D("clustersize_m_mod_dut" ,\
                    ";mod(x_{trk})_{2pitch} [#mum];mod(y_{trk})_{2pitch}[#mum];<cluster size>", \
                    40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TProfile2D("clustersize_m_mod_ref" ,";mod(x_{trk})_{2pitch} [mm];"\
                    "mod(y_{trk})_{2xpitch} [mm];<cluster size>", \
                    40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }
        self.hcharge_matched_mod = { minst.dut_plane: ROOT.TProfile2D("charge_m_mod_dut",\
                    ";mod(x_{trk})_{2pitch} [#mum]; mod(y_{trk})_{2pitch} [#mum];<cluster charge> [ADC]",\
                    40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TProfile2D("charge_m_mod_ref",\
                    ";mod(x_{trk})_{2pitch} [#mum]; mod(y_{trk})_{2xpitch} [#mum];<cluster charge> [ADC]",\
                    40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }
        self.hitmap_matched_mod = { minst.dut_plane: ROOT.TH2F("hitmap_m_mod_dut",\
                    ";mod(x_{trk})_{2pitch} [#mum]; mod(y_{trk})_{2pitch} [#mum];Entries",\
                    40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TH2F("hitmap_m_mod_ref",\
                    ";mod(x_{trk})_{2pitch} [#mum]; mod(y_{trk})_{2pitch} [#mum];Entries", \
                    40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }

        matched_histos = self.residual_matched.values()+\
                self.hcharge_matched.values()+self.hcharge1D_matched.values()+\
                self.hitmap_matched.values()+\
                self.hclustersize_matched.values()+self.hclustersize1D_matched.values()+\
                self.hclustersize_matched_mod.values()+self.hcharge_matched_mod.values()+\
                self.hitmap_matched_mod.values()

        # Correlations: 
        # -- Sensors-tracks
        self.hcorr_trkX = { minst.dut_plane: ROOT.TH2F("corr_trkX_dut",";x_{DUT} [mm]; x_{pred}^{trk} [mm]; Entries",200,-sxdut,sxdut,200,-sxdut,sxdut),
                minst.ref_plane: ROOT.TH2F("corr_trkX_ref",";x_{REF} [mm]; x_{pred}^{trk} [mm]; Entries",200,-sxref,sxref,200,-sxref,sxref)}
        ## Add it as alignment plot
        self._alignment_histos += self.hcorr_trkX.values()

        # -- DUT-REF
        self.hcorrX = ROOT.TH2F("corrX_dut_ref",";x_{DUT} [mm]; x_{REF} [mm]; Entries",100,-sxdut,sxdut,100,-sxref,sxref)
        self.hcorrY = ROOT.TH2F("corrY_dut_ref",";y_{DUT} [mm]; y_{REF} [mm]; Entries",100,-sydut,sydut,100,-syref,syref)
        
        # efficiency
        self.heff = ROOT.TProfile2D("eff_map",";x_{trk}^{DUT} [mm]; y^{DUT}_{trk} [mm];#varepsilon",\
                50,-1.1*sxdut,1.1*sxdut,int(2*nch_dut+3),-1.1*sydut,1.1*sydut)
        self.heff_ch = ROOT.TProfile2D("eff_map_ch",";x_{trk}^{DUT} [mm]; y^{DUT}_{trk} [channel];#varepsilon",\
                50,-1.1*sxdut,1.1*sxdut,int(2*nch_dut+3),-1.5,(nch_dut+1.0-0.5))
        self.heff_mod = ROOT.TProfile2D("eff_mod",";mod(x_{trk})_{2xpitch} [#mum]; mod(y_{trk})_{2xpitch} [#mum];"\
                        "efficiency", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM,0,1)
        dut_eff = [self.heff,self.heff_ch,self.heff_mod]

        # -- Eta distribution for matched, isolated
        self.heta = ROOT.TH1F("eta_2","#eta for isolated,matched cluster size 2;#eta;#frac{dN}{d#eta}",100,-0.1,1.1)
        self.heta_g = ROOT.TH1F("eta_g","#eta for isolated,matched cluster size #geq 2;#eta;#frac{dN}{d#eta}",100,-0.1,1.1)
        self.heta_csize = ROOT.TH2F("eta_cluster_size_dut","#eta vs. cluster size for isolated,matched;#eta;N_{cluster};Entries",\
                100,-0.1,1.1, 10,-0.5,9.5)
        self.hcl_size = { minst.dut_plane: ROOT.TH1F("cluster_size_dut","Isolated, REF-matched hits;N_{cluster};Entries",10,-0.5,9.5),
                minst.ref_plane: ROOT.TH1F("cluster_size_ref","Isolated-matched hits;N_{cluster};Entries",10,-0.5,9.5) }
        # -- Extra histos
        self.htrks_at_planes = { minst.dut_plane: ROOT.TH2F("trk_at_dut","Tracks at DUT;x_{DUT}^{trk} [mm]; y_{DUT}^{trk} [mm]; Entries",\
                        200,-10.0,10.0,200,-10.0,10.0),\
                minst.ref_plane: ROOT.TH2F("trk_at_ref","Tracks at REF;x_{REF}^{trk} [mm]; y_{REF}^{trk} [mm]; Entries",200,-10.0,10.0,200,-10.0,10.0)}
        self.evt_corr = { minst.dut_plane: ROOT.TProfile("dx_correlation_dut","Event correlation closest track;Event;(#Delta_x)_{DUT}",\
                        300000,0.0,299999),\
                    minst.ref_plane: ROOT.TProfile("dx_correlation_ref","Event correlation closest track;Event;(#Delta_x)_{REF}",\
                        300000,0.0,299999)}
        extra = self.htrks_at_planes.values()+self.evt_corr.values()
        
        # -- Diagnostics
        self.trk_iso = { minst.dut_plane: ROOT.TH1F("trkiso_dut","Distance between pair of tracks in the "\
                "same trigger-event;Track distance [mm];Triggers", 400,0,20.0*MM),\
                minst.ref_plane: ROOT.TH1F("trkiso_ref","Distance between pair of tracks in the same "\
                    "trigger-event;Track distance [mm];Triggers", 400,0,20.0*MM) }
        self.nhits_ntrks = { minst.dut_plane: ROOT.TH2F("nhits_ntracks_dut","at least 1 hit matched-isolated;N_{hits};N_{tracks};Triggers",10,-0.5,9.5,40,-0.5,39.5),
                minst.ref_plane: ROOT.TH2F("nhits_ntracks_ref","At least 1 hit matched-isolated ;N_{hits};N_{tracks};Triggers",10,-0.5,9.5,40,-0.5,39.5) }
        self.nhits_ntrks_all = { minst.dut_plane: ROOT.TH2F("nhits_ntracks_all_dut",";N_{hits};N_{tracks};Triggers",10,-0.5,9.5,40,-0.5,39.5),
                minst.ref_plane: ROOT.TH2F("nhits_ntracks_all_ref",";N_{hits};N_{tracks};Triggers",10,-0.5,9.5,40,-0.5,39.5) }
        self.hmatch_eff = { minst.dut_plane: ROOT.TProfile("match_eff_dut","Track-matching hit efficiency;x_{DUT};#varepsilon",\
                    100,-sxdut*1.01,sxdut*1.01,0,1),
                minst.ref_plane: ROOT.TProfile("match_eff_ref","Track-matching hit efficiency;x_{REF};#varepsilon",\
                    100,-sxref,sxref,0,1) }
        self.hy_eff = { minst.dut_plane: ROOT.TProfile("y_eff_dut","Track-matching hit efficiency vs. y-predicted;"\
                    "y_{DUT}^{pred}[mm];#varepsilon",100,-sydut,sydut,0,1),
                minst.ref_plane: ROOT.TProfile("y_eff_ref","Track-matching hit efficiency vs. y-predicted;"\
                    "y_{pred}[mm];#varepsilon",100,-syref,syref,0,1) }
        self.hntrk_eff = { minst.dut_plane: ROOT.TProfile("ntrk_eff_dut","Track-matching hit efficiency vs. number of tracks;"\
                    "N_{trk};#varepsilon",40,-0.5,39.5,0,1),
                minst.ref_plane: ROOT.TProfile("ntrk_eff_ref","Track-matching hit efficiency vs. number of tracks;"\
                    "N_{trk};#varepsilon",40,-0.5,39.5,0,1) }
        self.hnisotrk_eff = { minst.dut_plane: ROOT.TProfile("nisotrk_eff_dut","Track-matching hit efficiency vs. number of Non-isolated tracks;"\
                    "N_{trk} non-isolated;#varepsilon",8,-0.5,7.5,0,1),
                minst.ref_plane: ROOT.TProfile("nisotrk_eff_ref","Track-matching hit efficiency vs. number of tracks;"\
                    "N_{trk} non-isolated;#varepsilon",8,-0.5,7.5,0,1) }
        self.dutref_match_eff = ROOT.TProfile("match_dutref_eff","REF-matching DUT-hit efficiency;x_{DUT};#varepsilon",\
                                    100,-sxdut,sxdut,0,1)
        self.dutref_pure_match_eff = ROOT.TProfile("purematch_dutref_eff","REF-matching DUT-hit efficiency (REF hit present);x_{DUT};#varepsilon",\
                                    100,-sxdut,sxdut,0,1)
        self.hcorr_ref_dut = ROOT.TH2F("corr_ref_dut","Tracks matched with REF-hits and DUT-hits, at DUT position;"\
                "x_{DUT} [mm]; x_{REF} [mm]; Entries",100,-sydut,sydut,100,-syref,syref)
 
        diagnostics = []
        if DEBUG:
            diagnostics = self.trk_iso.values()+self.nhits_ntrks.values()+\
                    self.nhits_ntrks_all.values()+\
                    self.hmatch_eff.values()+self.hy_eff.values()+self.hntrk_eff.values()+self.hnisotrk_eff.values()+\
                    [self.dutref_match_eff,self.dutref_pure_match_eff,self.hcorr_ref_dut]


        # Keep track of all histograms which should be stored in the file
        self._allhistograms = self._alignment_histos+\
                associated_histos+\
                matched_histos+\
                [self.hcorrX,self.hcorrY]+\
                [self.heta,self.heta_g,self.heta_csize]+\
                dut_eff+diagnostics+\
                extra
        dummy=map(lambda h: h.SetDirectory(0),self._allhistograms)

        # The alignment constants
        self.alignment = {}

        # The eta values of the selected dut hits (isolated, matched with
        # a ref-hit), cluster size 2
        self.duthits_matched_eta = []
        # -- cluster size >=2
        self.duthits_matched_eta_h = []

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
        import ROOT
        from math import sin,cos,sqrt,asin,pi
        MIN_ENTRIES = 1000

        # Just extract all the static info of the class,
        for pl_id in self.alignment.keys():
            # an alignment instance to be populated
            new_align = alignment_cts(self.alignment[pl_id].sensitive_direction)

            if self.alignment[pl_id].sensitive_direction.lower() == "x":
                old_offset = self.alignment[pl_id].x_offset
            elif self.alignment[pl_id].sensitive_direction.lower() == "y":
                old_offset = self.alignment[pl_id].y_offset
            
            # -- The real offset
            coarse_offset =get_offset(self.dx_h[pl_id],-4.0,4.0)
            adding = False
            if self.alignment[pl_id].sensitive_direction.lower() == "x":
                old_offset = self.alignment[pl_id].x_offset
            elif self.alignment[pl_id].sensitive_direction.lower() == "y":
                old_offset = self.alignment[pl_id].y_offset

            if self.alignment[pl_id].iteration > 0 \
                    and (abs(coarse_offset-old_offset) < 0.1 \
                            and self.dx_finer_h[pl_id].Integral() > MIN_ENTRIES) :
                # PRe-alignment done, finer alignment now
                finer_offset = get_offset(self.dx_finer_h[pl_id],xmin=-0.2,xmax=0.2,coarse=False)
                if self.alignment[pl_id].sensitive_direction.lower() == "x":
                    new_align.x_offset = finer_offset
                elif self.alignment[pl_id].sensitive_direction.lower() == "y":
                    new_align.y_offset = finer_offset
                adding = True
                # -- Others things here --> 
                ## -- > rotation ()
                rot = -get_linear_fit(self.dy_x_h[pl_id],-3,3)
                # Asume resolution about 5 deg 
                if abs(rot) > 0.087:
                    new_align.rot = rot
                ## --> tilt, just if there is an initial inclination (at least 1 degree)
                #      otherwise, not evaluate nothing
                if hits_plane_accessor.sin_tilt(pl_id) > 0.84:
                    tilt = get_linear_fit(self.dy_y_h[pl_id],-3.,3)/hits_plane_accessor.sin_tilt(pl_id)
                    if abs(titl) > 0.087:
                        new_align.tilt = tilt
                ## --> turn, just if there is an initial inclination (at least 1 degree)
                #      otherwise, not evaluate nothing
                #if hits_plane_accessor.sin_turn(pl_id) > 0.84:
                #    tilt = get_linear_fit(self.dx_tx_h[pl_id],-3.,3)/hits_plane_accessor.sin_tilt(pl_id)
                #    if abs(titl) > 0.087:
                #        new_align.tilt = tilt
                ## --> dz
                dz = get_linear_fit(self.dx_tx_h[pl_id],-3,-3)*1e3
                ## --> Note that some detectors are not between telescope planes,
                #      meaning the dz is not well determined
                if abs(dz) < 10.0:
                    new_align.dz = dz

            # Update the alignment, correcting the old constants
            self.alignment[pl_id] += new_align
            
            if adding:
                print
                print "Extra alignment (PL:{0})".format(pl_id)
                print new_align
            else:
                # and overwrite the displacement
                if self.alignment[pl_id].sensitive_direction.lower() == "x":
                    self.alignment[pl_id].x_offset = coarse_offset
                elif self.alignment[pl_id].sensitive_direction.lower() == "y":
                    self.alignment[pl_id].y_offset = coarse_offset
                print
                print "Total alignment (PL:{0})".format(pl_id)
                print self.alignment[pl_id]
            
            # Update the module
            filename = ALIGN_FILE_FORMAT.format(pl_id,self.run_number)
            #new_align.update(filename,self.run_number)
            self.alignment[pl_id].update(filename,self.run_number)
            # And prepare the alignment constants to be updated 
            # in the next iteration
            hits_plane_accessor.resync[pl_id] = True

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

    def fill_alignment_histos(self,associated,trks,refs,duts):
        """Fill all the alignment histograms for the associated
        tracks and hits. 

        Parameters
        ----------
        associated: dict((int,(int,int))
            A track index associated with a reference and dut hit
            { track index: (ref-hit index,dut-hit index) }
        trks: tracks_accessor
            The tracks accessor instance
        refs: hits_plane_accessor
            The REF hits accessor instance
        duts: hits_plane_accessor
            The DUT hits accessor instance
        """
        for (itrk,(iref,idut)) in associated.iteritems():
            # Just do it for each sensor
            for i,hits in [(iref,refs),(idut,duts)]:
                # -- just matched 
                if i == -1:
                    continue
                # The sensitive axis
                if hits.sensitive_direction == "x":
                    ic = 0
                    iccom = 1
                    trk_drdz = trks.dxdz
                elif hits.sensitive_direction == "y":
                    ic = 1
                    iccom = 0
                    trk_drdz = trks.dydz
                (rpred,tel) = trks.get_point_in_sensor_frame(itrk,hits)
                dc = hits.sC_local[i]-rpred[ic]
                # -- rotation 
                self.dy_x_h[hits.id].Fill(rpred[iccom],dc)
                ## -- tilt
                self.dy_y_h[hits.id].Fill(rpred[ic],dc)
                # -- turn
                self.dx_xtx_h[hits.id].Fill(rpred[ic]*trk_drdz[itrk]*UM,dc)
                ## -- dz 
                self.dx_tx_h[hits.id].Fill(trk_drdz[itrk]*UM,dc)

    def fill_associated_hit_histos(self,(itrk,trks),(ihit,hits),ip):
        """Fill those histograms related with associated hit quantities
        (associated hit=a track is associated to the hit, see 
        self.associate_hits method)

        Parameters
        ----------
        itrk: int
            The index of the track to use
        trks: tracks_accessor
        ihit: int
            The index of the hit to use
        hits: hits_plane_accessor
        ip: int
            The index of the sensitive direction
        """
        # First of all check tha the index of the hit
        # is valid
        if ihit < 0:
            return
        r,tr = trks.get_point_in_sensor_frame(itrk,hits)
        # -- hit map
        self.hitmap_associated[hits.id].Fill(r[0],r[1])
        # -- residuals
        dr = hits.sC_local[ihit]-r[ip]
        self.residual_associated[hits.id].Fill(hits.sC_local[ihit],dr)
        # -- charge
        # -- Get the inverse of the sensitive direction:
        ip_compl = 1^ip
        self.hcharge_associated[hits.id].Fill(r[0],r[1],hits.charge[ihit])
        self.hcharge1D_associated[hits.id].Fill(hits.charge[ihit])
        # -- cluster size
        self.hclustersize_associated[hits.id].Fill(r[0],r[1],hits.n_cluster[ihit])
        self.hclustersize1D_associated[hits.id].Fill(hits.n_cluster[ihit])

        # -- Modulo 2 times the pitch
        xmod = ((r[0])%(2.0*hits.pitchX))
        ymod = ((r[1])%(2.0*hits.pitchY))
        # Charge, cluster and hit map for a region of the sensor
        self.hclustersize_associated_mod[hits.id].Fill(xmod*UM,ymod*UM,hits.n_cluster[ihit])
        self.hcharge_associated_mod[hits.id].Fill(xmod*UM,ymod*UM,hits.charge[ihit])
        self.hitmap_associated_mod[hits.id].Fill(xmod*UM,ymod*UM)

    def fill_matched_hit_histos(self,(itrk,trks),(ihit,hits),ip):
        """Fill those histograms related with matched hit quantities
        (matched hit=a track associated to a REF is associated to 
        the DUT hit , self.associate_hits method)

        Parameters
        ----------
        itrk: int
            The index of the track to use
        trks: tracks_accessor
        ihit: int
            The index of the hit to use
        hits: hits_plane_accessor
        ip: int
            The index of the sensitive direction
        """
        r,tr = trks.get_point_in_sensor_frame(itrk,hits)
        # -- hit map
        self.hitmap_matched[hits.id].Fill(r[0],r[1])
        # -- residuals
        dr = hits.sC_local[ihit]-r[ip]
        self.residual_matched[hits.id].Fill(hits.sC_local[ihit],dr)
        # -- charge
        # -- Get the inverse of the sensitive direction:
        ip_compl = 1^ip
        self.hcharge_matched[hits.id].Fill(r[0],r[1],hits.charge[ihit])
        self.hcharge1D_matched[hits.id].Fill(hits.charge[ihit])
        # -- cluster size
        self.hclustersize_matched[hits.id].Fill(r[0],r[1],hits.n_cluster[ihit])
        self.hclustersize1D_matched[hits.id].Fill(hits.n_cluster[ihit])
        
        # -- Modulo 2 times the pitch
        xmod = ((r[0])%(2.0*hits.pitchX))
        ymod = ((r[1])%(2.0*hits.pitchY))
        # Charge, cluster and hit map for a region of the sensor
        self.hclustersize_matched_mod[hits.id].Fill(xmod*UM,ymod*UM,hits.n_cluster[ihit])
        self.hcharge_matched_mod[hits.id].Fill(xmod*UM,ymod*UM,hits.charge[ihit])
        self.hitmap_matched_mod[hits.id].Fill(xmod*UM,ymod*UM)

    def fill_diagnosis_histos(self,trks,refhits,duthits,track_dict):
        """
        """
        # -- tracks at sensors
        for itrk in xrange(trks.n):
            # Track isolation
            trks.fill_isolation_histograms(itrk,refhits,self.trk_iso[refhits.id])
            trks.fill_isolation_histograms(itrk,duthits,self.trk_iso[duthits.id])
            ((xpred_ref,ypred_ref,zpred_ref),tel_ref) = trks.get_point_in_sensor_frame(itrk,refhits)
            self.htrks_at_planes[refhits.id].Fill(xpred_ref,ypred_ref)
            ((xpred_dut,ypred_dut,zpred_dut),tel_ref) = trks.get_point_in_sensor_frame(itrk,duthits)
            self.htrks_at_planes[duthits.id].Fill(xpred_dut,ypred_dut)
        # Evaluate REF-DUT matching
        for (itrk,(iref,idut)) in track_dict.iteritems():
            # -- Check what's the probability of REF association, i.e. 
            if idut == -1:
                continue
            self.dutref_match_eff.Fill(duthits.sC_local[idut],int(iref != -1))
        #    # Track efficiencies incorporated
        #    if len(matched_hits[refhits.id]) > 0:
        #        self.dutref_pure_match_eff.Fill(duthits.sC_local[ihit_dut],(len(ref_match) > 0))
        #    # Whats the position the track matched with the REF gives in the DUT?
        #    #for ihit_ref,itrk_r in matched_hits[refhits.id]:
        #    #    self.dutref_distance.Fill(trks.get_point_in_sensor_frame(itrk_d,duthits)[0][ic]-\
        #    #            trks.get_point_in_sensor_frame(itrk_r,refhits)[0][ic])
    
    def fractionary_position_plot(self):
        """Using the eta-distribution obtained 
        in the processing, calculate the fractionary
        position for size 2 clusters and beyond
        """
        import ROOT
        from .analysis_functions import frac_position
        # Create a new histogram
        self.hfrac_pos_2 = ROOT.TH2F("hfrac_pos_2","Fractionary cluster position for"\
                " cluster-2 size;mod(x)_{pitch};charge [ADC];clusters",100,-0.01,1.01,100,0,1000)
        self.hfrac_pos_2.SetDirectory(0)
        self.hfrac_pos = ROOT.TH2F("hfrac_pos","Fractionary cluster position for"\
                " cluster #geq 2 size;mod(x)_{pitch};charge [ADC];clusters",100,-0.01,1.01,100,0,1000)
        self.hfrac_pos.SetDirectory(0)
        self._allhistograms.append(self.hfrac_pos_2)
        self._allhistograms.append(self.hfrac_pos)
        # Get the loop
        dummy = map(lambda (eta,charge): self.hfrac_pos_2.Fill(frac_position(eta,self.heta),charge),\
                self.duthits_matched_eta)
        dummy = map(lambda (eta,charge): self.hfrac_pos.Fill(frac_position(eta,self.heta),charge),\
                self.duthits_matched_eta_h)
        # done :)


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
        if tracks != 0:
            self.total_events_tracks += 1
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
                self.events_with_matched_hits[plid] = int(len(thelist)>0)
            else:
                self.number_matched_hits[plid] += len(thelist)
                self.events_with_matched_hits[plid] += int(len(thelist)>0)

    def process(self,t,trks,refhits,duthits,is_alignment=False):
        """Process an event. Per each track in the event which is 
        isolated, and tries to associate a hit in the REF sensor 
        (in that case defines a particle) and in the DUT sensor.
        Fill the relevant counters and histograms.

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
        self.fill_statistics(duthits.n,refhits.n,trks.n)

        # Start new event
        trks.new_event()
        refhits.new_event()
        duthits.new_event()
        
        # ------------------------------------------------------------
        # RECALL definitions:           
        # - Associated track: a track with a matched (track-hit 
        #   distance lower than a threshold) hit at the sensor plane
        # - Matched track == particle: a track associated to a REF hit,
        #   and associated to a DUT hit
        # ------------------------------------------------------------
        
        # - a dict to get the sensors by its plane ID
        sensor_hits = { refhits.id: refhits, duthits.id: duthits }

        # --- Obtain the number of hits and tracks in the event
        self.nhits_ntrks_all[refhits.id].Fill(refhits.n,trks.n)
        self.nhits_ntrks_all[duthits.id].Fill(duthits.n,trks.n)
        
        # Prepare the ntuple of histograms for DUT and REF
        histos_ref = (self.hcorr_trkX[refhits.id],self.dx_h[refhits.id],\
                self.dx_finer_h[refhits.id],self.dx_finer_wide_h[refhits.id],self.hplane[refhits.id])
        histos_dut = (self.hcorr_trkX[duthits.id],self.dx_h[duthits.id],\
                self.dx_finer_h[duthits.id],self.dx_finer_wide_h[duthits.id],self.hplane[duthits.id])
        histos = { refhits.id: histos_ref, duthits.id: histos_dut }
        # -- Get the sensitive axis, first check an evidence
        assert(duthits.sensitive_direction == refhits.sensitive_direction)
        if duthits.sensitive_direction == "x":
            ic = 0
        elif duthits.sensitive_direction == "y":
            ic = 1

        # -- The workhorse method: obtain one hit per track
        track_dict = trks.associate_hits(refhits,duthits,histos)

        # -- If alignment just perform the histo filling and finish
        if is_alignment:
            # -- fill the relevant histograms 
            self.fill_alignment_histos(track_dict,trks,refhits,duthits)
            # -- update the alignment constants
            self.alignment = dict(map(lambda (i,d): (i,d.align_constants[i]) ,sensor_hits.iteritems()))
            #Get an estimation of the sensor efficiency ?
            #self.fill_statistics_matched(matched_hits)
            return
        elif len(self._alignment_histos) > 0: 
            # Remove the aligment histos (only do it the first time)
            # (except some of them) and from the generic counter list
            keepthem =[]
            for hname in ["corr_trkX_{0}","dx_{0}","dx_finer_{0}","dx_finer_wide_{0}","plane_{0}"]:
                for sname in [ "dut","ref" ]:
                    keepthem.append(hname.format(sname)) 
            dummy = map(lambda h: (h.Delete(),self._allhistograms.remove(h)),\
                        filter(lambda _h: _h.GetName() not in keepthem,self._alignment_histos))
            # Empty the list to avoid enter again
            self._alignment_histos = []
        
        # -- Analysis efficiency and fill histograms:
        #    Associated and matched 
        for (itrk,(iref,idut)) in track_dict.iteritems():
            # -- Fill extra histograms: associated hits 
            self.fill_associated_hit_histos((itrk,trks),(iref,refhits),ic)
            self.fill_associated_hit_histos((itrk,trks),(idut,duthits),ic)
            # -- No tracks (i.e., REF-matched particle) presence 
            if iref == -1:
                continue
            # -- Fill matched histograms
            r_at_dut,tel_at_dut = trks.get_point_in_sensor_frame(itrk,duthits)
            # Just within acceptance
            if not duthits.is_within_fiducial(r_at_dut[0],r_at_dut[1]):
                continue
            # -- Efficiency plots
            self.heff.Fill(r_at_dut[0],r_at_dut[1],(idut != -1))
            self.heff_ch.Fill(r_at_dut[0],duthits.get_y_channel(r_at_dut[1]),(idut != -1))
            # Modulo pitch
            xmod = ((r_at_dut[0])%(2.0*duthits.pitchX))*UM
            ymod = ((r_at_dut[1])%(2.0*duthits.pitchY))*UM
            self.heff_mod.Fill(xmod,ymod,(idut != -1))
            # --- Fill /
            if idut == -1:
                continue
            # Fill the matched histograms
            self.fill_matched_hit_histos((itrk,trks),(iref,refhits),ic)
            self.fill_matched_hit_histos((itrk,trks),(idut,duthits),ic)
            # --- Some histograms exclusivelly for particles at DUT
            # -- Eta distribution for the matched cluster size=2
            clustersize = duthits.n_cluster[duthits.track_link[itrk]]
            if clustersize >= 2:
                self.heta_g.Fill(duthits.eta[idut])
                self.heta_csize.Fill(duthits.eta[idut],clustersize)
                self.duthits_matched_eta_h.append((duthits.eta[idut],duthits.charge[idut]))
                if clustersize == 2:
                    self.heta.Fill(duthits.eta[idut])
                    self.duthits_matched_eta.append(self.duthits_matched_eta_h[-1])
        # Extra histograms: 
        if DEBUG:
            self.fill_diagnosis_histos(trks,refhits,duthits,track_dict)

        #self.fill_statistics_matched(matched_hits)

    def get_raw_sensors_efficiency(self):
        """Summarize the efficiency of the sensors:
        eff = Number of matched hits/number of events with tracks
        """
        m = ""
        for plid,n in self.number_matched_hits.iteritems():
            evts_w_match = self.events_with_matched_hits[plid]
            if plid == self.dut_plane:
                sensor_evts = self.events_with_dut
                sensor_name = "DUT"
            else:
                sensor_evts = self.events_with_ref
                sensor_name = "REF"
            #m+= "{0} efficiency (Isolated track-matched): {1:.2f}%\n".format(sensor_name,float(n)/float(self.total_events_tracks)*100.)
            m+= "{0} efficiency (Isolated track-matched): {1:.2f}%\n".format(sensor_name,float(evts_w_match)/float(self.total_events_tracks)*100.)
        return m

    def __str__(self):
        """Summarize the information
        """
        m = "|-------------------------------------------|\n"
        m+= " Events with DUT: {0} (eff: {1:.2f}%)\n".\
                format(self.events_with_dut,float(self.events_with_dut)/float(self.total_events)*100.)
        m+= " Events with REF: {0} (eff: {1:.2f}%)\n".\
                format(self.events_with_ref,float(self.events_with_ref)/float(self.total_events)*100.)
        m+= "\n Events with DUT but no tracks: {0}  (eff: {1:.2f}%)\n".\
                format(self.events_with_dut_no_tracks,float(self.events_with_dut_no_tracks)/float(self.total_events)*100.0)
        m+= " Events with REF but no tracks: {0} (eff: {1:.2f}%)\n".\
                format(self.events_with_ref_no_tracks,float(self.events_with_ref_no_tracks)/float(self.total_events)*100.0)
        m+= " Events with REF and DUT and tracks: {0} (eff:{1:.2f}%)\n".\
                format(self.events_with_ref_dut_tracks,float(self.events_with_ref_dut_tracks)/float(self.total_events)*100.0)
        m+= "\n Total processed events: {0}\n\n".format(self.total_events)
        m+= self.get_raw_sensors_efficiency()
        m+= "|-------------------------------------------|"

        return m

def get_offset(h,xmin=-2.0,xmax=2.0,coarse=True):
    """Fit the histogram to a gaussian plus a linear background
    using a Chebyshev polynomial of 2nd order.

    Parameters
    ----------
    h: ROOT.TH1F
        The histogram to fit
    xmin: float, optional
        The fit range minimum
    xmax: float, optional
        The fit range maximum
    coarse: bool    
        Whether or not perform a previous subtraction of a 2nd
        order polynomial background. This assumes a huge background
        contribution which could give higher frequencies away from 
        the offset

    Return
    ------
    new_align_x: flot
        The peak of the histogram where the offset is
    """
    import ROOT
    ROOT.gROOT.SetBatch()
    
    cns=ROOT.TCanvas()

    if h.GetEntries() == 0:
        raise RuntimeError("Empty histogram '{0}'".format(h.GetName()))

    # Just do it when needed
    if coarse:
        # -- There is cases where the peak is difficult to obtain 
        #    with the background presence: first evaluate the background 
        #    and remove it
        #    (o gaussian?)
        bkg = ROOT.TF1("bkg","cheb2",xmin,xmax)
        #bkg = ROOT.TF1("bkg","gausn",xmin,xmax)
        status_bkg = h.Fit(bkg,"SQ","",xmin,xmax)
        try:
            if status_bkg.Status() != 0:
                raise RuntimeError("\033[1;33mWARNING\033[1;m FAILED "\
                        "THE X-OFFSET FIT for '{0}': "\
                        "Status {1}".format(h.GetName(),status.Status()))
        except ReferenceError:
            raise RuntimeError("\033[1;32mERROR\033[1;m FAILED THE X-OFFSET FIT for '{0}': "\
                    "No data was fitted".format(h.GetName()))
        # -- Extract the background to obtain the peak
        hsub = ROOT.TH1F("h_subtracted_"+str(hash(cns)),"",h.GetNbinsX(),xmin,xmax)
        #hsub.SetDirectory(0)
        for i in xrange(1,hsub.GetNbinsX()+1):
            x = hsub.GetBinCenter(i)
            ysub = h.GetBinContent(i)-bkg.Eval(x)
            if ysub > 0:
                hsub.SetBinContent(i,ysub)
    else:
        # - fine adjustment, 
        hsub = h
    # -- Now find the peak (if needed)
    gbg = ROOT.TF1("gbg","gausn(0)+cheb2(3)",xmin,xmax)
    # Set the initial values
    # -- Amplitude
    gbg.SetParameter(0,h.GetMaximum())
    # -- Mean (find the peak from the background subtracted
    if coarse:
        peak = h.GetBinCenter(hsub.GetMaximumBin())
    else:
        # It supose to be centered
        peak = 0.0
        # -- Some extra requirements: at least 50 entries per 
        #    bin in average
        while float(h.Integral())/float(h.GetNbinsX()) < 50:
            h.Rebin(2)

    gbg.SetParameter(1,peak)
    # -- Sigma
    gbg.SetParameter(2,hsub.GetBinWidth(1))
    # -- Background: guess it from a region far away from the peak
    #    only when coarse
    if coarse:
        gbg.SetParameter(3,hsub.GetBinContent(h.FindBin(peak-0.4)))
    ## Do the fit
    status = h.Fit(gbg,"SQR","")
    try:
        if status.Status() != 0:
            print "\033[1;33mWARNING\033[1;m FAILED THE X-OFFSET FIT for '{0}': "\
                    "Status {1}".format(h.GetName(),status.Status())
    except ReferenceError:
        print "\033[1;32mERROR\033[1;m FAILED THE X-OFFSET FIT for '{0}': "\
                "No data was fitted".format(h.GetName())
    new_align_x = gbg.GetParameter(1)

    return new_align_x

def get_linear_fit(h,xmin=-2.0,xmax=2.0,robval=0.75):
    """Perform a linear (robust) fit of the histogram (profile)

    Parameters
    ----------
    h: ROOT.TH1
        The histogram or the profile 1D
    xmin: float, optional (None)
        The minimum x-value to consider in the fit
    xmax: float, optional (None)
        The maximum x-value to consider in the fit
    robval: float, optional (0.75)
        The percentage of data points that are good points
    """
    import ROOT
    import array

    ROOT.gROOT.SetBatch()
    cns=ROOT.TCanvas()
    
    if h.GetEntries() == 0:
        raise RuntimeError("Empty histogram '{0}'".format(h.GetName()))
    gbg = ROOT.TF1("kkita","pol1")
    status = h.Fit(gbg,"RQS","",xmax,xmin)
    if status.Status() != 0:
        print "\033[1;33mWARNING\033[1;m FAILED THE LINEAR FIT for '{0}': "\
                "Status {1}".format(h.GetName(),status.Status())
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
    proc_inst = sensor_map_production(fname,entries_proc=60000,alignment=True,verbose=verbose)
    # Considered aligned if the aligned constants doesn't change (within some 
    # tolerance): then extract the alignment constants before
    aligned_sensors = 0
    for pl_id,current_align in proc_inst.alignment.iteritems():
        align_f = ALIGN_FILE_FORMAT.format(pl_id,proc_inst.run_number)
        old_align = alignment_cts(current_align.sensitive_direction)
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
    print
    print proc_inst.get_raw_sensors_efficiency()
    # Check if all sensors are aligned: XXX---Not correct---XXX
    return (aligned_sensors == len(proc_inst.alignment.keys()))


def sensor_map_production(fname,entries_proc=-1,alignment=False,verbose=False):
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
    verbose: bool
        Whether to activate the debug mode, with more print-outs

    Return
    ------
    A processor instance (if alignment mode)

    Raises
    ------
    IOError
        Whenever the ROOT file is not present
    """
    import ROOT
    import sys
    import os
    from .SPS2017TB_metadata import filename_parser 
    from .SPS2017TB_metadata import standard_sensor_name_map as name_converter
    from .SPS2017TB_metadata import sensor_name_spec_map as specs
    # probably provisional XXX

    global DEBUG
    DEBUG=verbose
    if DEBUG:
        import timeit

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
    dut = hits_plane_accessor(t,metadata.dut_plane,sensor_name=name_converter[fp.sensor_name])
    ref = hits_plane_accessor(t,metadata.ref_plane,"REF_0_b1")
    tracks = tracks_accessor(t,[0,1,2,3,4],metadata.ref_plane,metadata.dut_plane)
    
    # Process the data
    if DEBUG:
        start_time = timeit.default_timer()
    wtch = processor(metadata,dut.run_number)
    
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
    if DEBUG:
        print timeit.default_timer()-start_time
    if alignment:
        if DEBUG:
            foutput = "alignment_plots_{0}_run000{1}_{2}.root".format(fp.sensor_name,fp.run_number,int(wtch.alignment.values()[0].iteration))
        else:
            foutput = "alignment_plots_{0}_run000{1}.root".format(fp.sensor_name,fp.run_number)
        wtch.store_alignment(foutput)
        return wtch
    else:
        print wtch
        # Calculate the fractionary position 
        wtch.fractionary_position_plot()
        foutput = "sensor_maps_{0}_run000{1}.root".format(fp.sensor_name,fp.run_number)
        print "File created at '{0}'".format(foutput)
        wtch.store_histos(foutput)
