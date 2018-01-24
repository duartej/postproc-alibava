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

DEBUG=False

# Be sure use the same alignment files: ct
ALIGN_FILE_FORMAT= "aligncts_plane_{0}_{1}.txt"

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
    hyp_par = ( (0,("x_offset",1.0)), \
                    (1,("dz",-1.0)), \
                    (2,("tilt",-1.0)),\
                    (3,("turn",1.0)),\
                    (4,("rot",-1.0)) )
    def __init__(self):
        """Initialize the data-members
        """
        self.iteration = int(0)
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
        iteration = self.iteration
        for attr in self.__dict__.keys():
            # Take into account the sign
            try:
                factor = filter(lambda (i,(name,fact)): name == attr,self.hyp_par)[0][0]
            except IndexError:
                # y_offset
                factor = 1.0
            prov = getattr(self,attr)+getattr(other,attr)
            setattr(self,attr,prov)
        # next iterations, just overwritte the wrong value
        # when added from the other
        self.iteration = iteration+1
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
    eta: ROOT.std.vector(float)()
        The eta value (charge distribution on the cluster)
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
        import array

        self.id = planeid
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
        #dummy = tree.GetEntry(0)
        ## Fill the elements
        #self.x = getattr(tree,"hit_X_{0}".format(self.id))
        #self.y = getattr(tree,"hit_Y_{0}".format(self.id))
        #self.x_local = getattr(tree,"hit_XLocal_{0}".format(self.id))
        #self.y_local = getattr(tree,"hit_YLocal_{0}".format(self.id))
        ##self.x_strips = lambda i: (self.x_local[i]+sensor_xsize/2.0)/pitch+0.5
        #self.charge = getattr(tree,"hit_total_charge_{0}".format(self.id))
        #self.n_cluster = getattr(tree,"hit_Ncluster_{0}".format(self.id))
        #self.eta = getattr(tree,"hit_cluster_eta_{0}".format(self.id))

        #self.run_number = tree.RunNumber
        self.x = ROOT.std.vector(float)()
        tree.SetBranchAddress("hit_X_{0}".format(self.id),self.x)
        self.y = ROOT.std.vector(float)()
        tree.SetBranchAddress("hit_Y_{0}".format(self.id),self.y)
        self.x_local = ROOT.std.vector(float)()
        tree.SetBranchAddress("hit_XLocal_{0}".format(self.id),self.x_local)
        self.y_local = ROOT.std.vector(float)()
        tree.SetBranchAddress("hit_YLocal_{0}".format(self.id),self.y_local)
        #self.x_strips = lambda i: (self.x_local[i]+sensor_xsize/2.0)/pitch+0.5
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
                hits_plane_accessor.align_constants[self.id] = alignment_cts()
            filename = ALIGN_FILE_FORMAT.format(self.id,self.run_number)
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
                # And give an extra turn allowing to perform an initial 
                # adjustment (1.15 degrees)
                #hits_plane_accessor.align_constants[self.id].turn = 0.021
            # Not re-synchronized until new order
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

    def x_channel(self,i,sizeX,pitch):
        """The equivalent channel number for the
        i-hit

        Parameters
        ----------
        i: int
            The index of the hit
        sizeX: float
            The size of the sensor in mm
        pitch: float
            The pitch of the segmented channels in mm

        Return
        ------
        The x position in channel number
        """
        return (self.x_local[i]+0.5*sizeX)/pitch-0.5
    
    def y_channel(self,i,sizeY,pitch):
        """The equivalent channel number for the
        i-hit

        Parameters
        ----------
        i: int
            The index of the hit
        sizeX: float
            The size of the sensor in mm
        pitch: float
            The pitch of the segmented channels in mm

        Return
        ------
        The y position in channel number
        """
        return (self.y_local[i]+0.5*sizeY)/pitch-0.5
    
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
        other: hit_plane_accessor instance (or any other instance
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
              (see last histo) vs. the x-prediction. See `dx_y_h` in the 
              processor class. This histo is used for alignment (rotation)            
            - A ROOT.TProfile to store the variation of the residual in x
              (see last histo) vs. the y-prediction times the slope of the
              track (x). See `dx_ytx_h` in the processor class. This histo
              is used for alignment (tilt)
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
                hmatch_eff,hy_eff,hntracks_eff,hnisotracks_eff,fitter = histos

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
                # Isolation: be sure there is no other track surrounding this one
                ((xpred,ypred,zpred),rtel0) = track_acc.get_point_in_sensor_frame(itrk,self)
                dummy = map(lambda ((o_x,o_y,o_z),_tel): hiso.Fill(sqrt((o_x-xpred)**2.0+(o_y-ypred)**2.0)),\
                        map(lambda other_i:  track_acc.get_point_in_sensor_frame(other_i,self),\
                        not_used_track_indices[i_at_list+1:]))
                # note tha i_at_list catchs the indice at `not_used_track_indices`
                # list, therefore in order to evaluate only those tracks not previously
                # evaluated just use the remaining elements after the current track index
                non_isolated = filter(lambda (oindex,((o_x,o_y,o_z),_tel)): sqrt((o_x-xpred)**2.0+(o_y-ypred)**2.0) < ISOLATION,\
                        map(lambda other_i:  (other_i,track_acc.get_point_in_sensor_frame(other_i,self)),\
                            not_used_track_indices[i_at_list+1:]))
                if len(non_isolated) > 0 :
                    # Another track has been found sourrounding this one, 
                    # not isolated, remove them from the available tracks
                    used_tracks.append(itrk)
                    dummy = map(lambda (i,other_stuff): used_tracks.append(i), non_isolated)
                    continue
                # Fill the alignment histograms
                closest[abs(self.x_local[ihit]-xpred)] = itrk
            if len(closest) == 0:
                continue
            ### The hit was matched by at least one track
            ##hmatch_eff.Fill(self.x_local[ihit],1)
            # Get the closest track (in x)
            distance_abs,trk_el =sorted(closest.iteritems())[0]
            #rturn,rtitl,rrot,rsensor = track_acc.get_point_in_sensor_frame(trk_el,self)
            rsensor,rtel = track_acc.get_point_in_sensor_frame(trk_el,self)
            # Note that the prediction is given in the sensor reference frame (z should be zero)
            # Therefore, not to interesting this plot, better
            hplane.Fill(rtel[2]-self.z[0],rtel[0],rtel[1])
            # Now we can check in the telescope plane, the discrepancies is due
            # to the misalignment
            #rsensor = track_acc.get_point(trk_el,self.z[ihit])
            # And fill some histograms
            # Fill correlation histograms using only the closest isolated tracks
            hcorr.Fill(self.x_local[ihit],rsensor[0])
            # -- resolution histogram
            dx = self.x_local[ihit]-rsensor[0]
            ###dx = self.x[ihit]-rsensor[0]
            hres.Fill(self.x_local[ihit],dx)
            # Alignment histos: x-offset
            # -- Pre-alignment, rot, tilt an turn but no shift
            hdx.Fill((self.x_local[ihit]+self.align_constants[self.id].x_offset)-rsensor[0])
            ###hdx.Fill(self.x[ihit]-rsensor[0])
            # -- finer
            hdx_finer.Fill(dx)
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
                hrot.Fill(rsensor[1],dx)
                ## -- tilt
                htilt.Fill(rsensor[1]*track_acc.dxdz[trk_el]*UM,dx)
                # -- turn
                hturn.Fill(rsensor[0]*track_acc.dxdz[trk_el]*UM,dx)
                ## -- dz 
                hdz.Fill(track_acc.dxdz[trk_el]*UM,dx)
                ## -- to perform the alignment fit
                coord = array.array('d',[-track_acc.dxdz[trk_el],-track_acc.dxdz[trk_el]*rsensor[1],\
                        track_acc.dxdz[trk_el]*rsensor[0],-rsensor[1]])
                fitter.AddPoint(coord,dx)
            # The hit was matched by at least one track
            hmatch_eff.Fill(self.x_local[ihit],int(distance_abs<res))
            hy_eff.Fill(rsensor[1],int(distance_abs<res))
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
        ##self.x = getattr(tree,"{0}_X_{1}".format(prefix,self.id))
        ##self.y = getattr(tree,"{0}_Y_{1}".format(prefix,self.id))
        ##self.x_local = getattr(tree,"{0}_XLocal_{1}".format(prefix,self.id))
        ##self.y_local = getattr(tree,"{0}_YLocal_{1}".format(prefix,self.id))
        ##self.z     = getattr(tree,"{0}_Z_{1}".format(prefix,self.id))
        ##self.track_index = getattr(tree,"{0}_index_{1}".format(prefix,self.id))
    
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
        ##self.x0 = tree.trk_refPoint_X
        ##self.y0 = tree.trk_refPoint_Y
        ##self.z0 = tree.trk_refPoint_Z
        ### Director vector
        ##self.dxdz = tree.trk_dxdz
        ##self.dydz = tree.trk_dydz
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
        ##self.dydz = tree.trk_dydz
        ##self.dydz = tree.trk_dydz
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
        # -- dict to be used to list of the non-isolated tracks
        self._cache_nonisolated = {}
    
    def is_isolated(self,i,hit_z_position,isolation_cone=0.6*MM):
        """Whether or not the track is isolated, i.e.
        if there is no other track in the plane defined by the 
        z position inside the isolation cone

        Parameters
        ----------
        i: int
            The index of the track
        hit_z_position: hit_accessor
            The hit accessor which defines the plane
        isolation_code: float, optional
            The size of the isolation cone

        Return
        ------
        bool: False if there is any track inside the isolation
              cone
        """
        from math import sqrt

        if self._cache_nonisolated.has_key(i):
            return False
        # -- Extract all the indices of isolated tracks (so far)
        not_used_track_indices = filter(lambda _itr: i != _itr and _itr not in self._cache_nonisolated.keys(), xrange(self.n))
        # -- Get the track at the plane
        ((xpred,ypred,zpred),rtel0) = self.get_point_in_sensor_frame(i,hit_z_position)
        # -- Get the list of non isolated
        non_isolated = filter(lambda (oindex,((o_x,o_y,o_z),_tel)): sqrt((o_x-xpred)**2.0+(o_y-ypred)**2.0) < isolation_cone,\
                map(lambda other_i:  (other_i,self.get_point_in_sensor_frame(other_i,hit_z_position)),\
                not_used_track_indices))
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
        """Fill the isolation histograms, and for this need the x and y 
        prediction of the track at the plane of the sensor, therefore
        as the calculation must be done, it returns the predicted point
        XXX -- SURE ? Mixing two different things ... Remember 
        get_point_
        """
        from math import sqrt

        ((xpred,ypred,zpred),rtel0) = self.get_point_in_sensor_frame(itrk,hitobj)
        dummy = map(lambda ((o_x,o_y,o_z),_tel): h.Fill(sqrt((o_x-xpred)**2.0+(o_y-ypred)**2.0)),\
                map(lambda other_i:  self.get_point_in_sensor_frame(other_i,hitobj),xrange(itrk+1,self.n)))
    
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
        from math import sqrt

        # |--- Telescope coordinates -----------------------------------------|
        # Normal vector of the plane (defined in opposite direction to the 
        # beam). With perfect alignment would be: (0,0,-1)
        n = hit.get_normal_vector()
        #  - the sensor position: (including alignment corrections)
        rpp = (hit.align_constants[hit.id].x_offset,\
                hit.align_constants[hit.id].y_offset,\
                hit.z[i]+hit.align_constants[hit.id].dz)
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
        # Re-write rpp in zc coordinates (zcc-z0)
        #rppc = (rpp[0],rpp[1],rpp[2]-zc)
        # The point 
        rpred = track(zc)

        # |--- Sensor plane coordinates --------------------------------------| 
        # As the sensor reference system is defined (wrt to the telescope 
        # plane) as  e0=(-1,0,0), e1=(0,1,0), n=(0,0,-1), we will need to 
        # correct the sign in the transformations 
        sign_c = (-1.0,1.0,-1.0)
        
        # Getting the hit position vector in the sensor coordinates:
        # a = track-rpp
        a = matvectmult(hit.matrix_sensor_rec(hit.id), \
                map(lambda (i,x): sign_c[i]*(rpred[i]-x),enumerate(rpp)))
        
        ### -- Now, let's obtain the point (which is in the telescope reference plane)
        ###    into the sensor reference plane in order to compare this point with
        ###    the hit in the plane

        ### -- 1st. find the displacement w.r.t the sensor in z (dz)
        #dzc = zc-rpp[2]

        ## -- 2nd. turn around y (-omega)
        #r_turn = (hit.cos_turn(hit.id)*rpred[0]-hit.sin_turn(hit.id)*dzc,\
        #        rpred[1],\
        #        hit.sin_turn(hit.id)*rpred[0]+hit.cos_turn(hit.id)*dzc)
        ## -- 3rd. tilt around x (-alpha) (Note that the r_tilt[2] should be zero,
        ##    in the dut plane)
        #r_tilt = (r_turn[0], \
        #        hit.cos_tilt(hit.id)*r_turn[1]+hit.sin_tilt(hit.id)*r_turn[2],\
        #        -hit.sin_tilt(hit.id)*r_turn[1]+hit.cos_tilt(hit.id)*r_turn[2])
        #if abs(r_tilt[2]) > 1e-9:
        #    raise RuntimeError("Z-component is not ZERO!!!! INCONSISTENT ERROR.")


        ## -- 4th. rotation around z
        #r_rot = (hit.cos_rot(hit.id)*r_tilt[0]+hit.sin_rot(hit.id)*r_tilt[1],\
        #        -hit.sin_rot(hit.id)*r_tilt[0]+hit.cos_rot(hit.id)*r_tilt[1],\
        #        r_tilt[2])

        ## -- 5th. Shifts
        #r = (r_rot[0]-hit.align_constants[hit.id].x_offset,\
        #        r_rot[1]-hit.align_constants[hit.id].y_offset,\
        #        r_rot[2])
        ###print 
        ###print r,rpred

        ##print
        ###print "Telescope , Sensor frame:"
        ### The rpp
        #print "KK",a,r,map(lambda (i,(_a,_b)): _a-sign_c[i]*_b,enumerate(zip(a,r)))
    
        ##print "N:",n,n_u
        ##print "Sensor position:",rpp,rpp_u
        ##print "Track reference point:",r0,r0_u
        ##print "Track slope",ts,ts_u
        ##print "Crossing point z",zc,zc_u
        ##print "Position:",rpred,tuple(map(lambda i: r0_u[i]+ts_u[i]*zc_u,xrange(len(ts_u))))
        ##print "Position the same should be:",r,tuple(map(lambda i: r0_u[i]+ts_u[i]*zc_u,xrange(len(ts_u))))
        
        return a,rpred #r,rpred

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
        self.total_events_tracks = 0
        self.events_with_dut = 0
        self.dut_hits        = 0
        self.events_with_ref = 0
        self.ref_hits        = 0
        self.events_with_dut_no_tracks = 0
        self.events_with_ref_no_tracks = 0
        self.events_with_ref_dut_tracks = 0
        self.number_matched_hits = {}
        # Alingment histos
        # -----------------
        # -- the shift in x
        self.dx_h = { minst.dut_plane: ROOT.TH1F("dx_dut"," ; x_{DUT}-x_{trk}^{pred} [mm]; Entries",100,-sxdut,sxdut),
                minst.ref_plane: ROOT.TH1F("dx_ref"," ; x_{REF}-x_{trk}^{pred} [mm]; Entries",100,-sxref,sxref) }
        # -- finer alignment
        self.dx_finer_h = { minst.dut_plane: ROOT.TH1F("dx_finer_dut"," ; x_{DUT}-x_{trk}^{pred} [mm]; Entries",100,-0.5,0.5),
                minst.ref_plane: ROOT.TH1F("dx_finer_ref"," ; x_{REF}-x_{trk}^{pred} [mm]; Entries",100,-0.5,0.5) }
        # -- the rot (around z-axis)
        self.dx_y_h = { minst.dut_plane: ROOT.TProfile("dx_y_dut"," ;y_{trk}^{pred} [mm];#Deltax_{DUT} [mm]",\
                        50,-6.0,6.0,-0.2,0.2),
                minst.ref_plane: ROOT.TProfile("dx_y_ref"," ;y_{trk}^{pred} [mm];#Deltax_{REF} [mm]",\
                        50,-6.0,6.0,-0.2,0.2) }
        # -- the tilt (around x-axis)
        self.dx_ytx_h = { minst.dut_plane: ROOT.TProfile("dx_ytx_dut"," ;y_{trk}^{pred}*tan(#theta_{x}) [#mum];#Deltax_{DUT} [mm]",\
                        50,-0.1,0.1,-0.2,0.2),
                minst.ref_plane: ROOT.TProfile("dx_ytx_ref"," ;y_{trk}^{pred}*tan(#theta_{x}) [#mum];#Deltax_{REF} [mm]",\
                        50,-0.1,0.1,-0.2,0.2) }
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

        # The linear fit to the hyperplane (Delta x)
        self.dx_hyp = { minst.dut_plane: ROOT.TLinearFitter(4), minst.ref_plane: ROOT.TLinearFitter(4) }
        dummy = map(lambda f: f.SetFormula("hyp4"), self.dx_hyp.values())
        
        self._alignment_histos = self.dx_h.values()+self.dx_finer_h.values()+self.dx_y_h.values()+\
                self.dx_xtx_h.values()+self.dx_ytx_h.values()+self.dx_tx_h.values()+\
                self.hplane.values()

        # Analysis histos using isolated tracks
        # -------------------------------------
        # Residuals between REF-DUT (matched tracks family)
        self.residual_sensor = ROOT.TH2F("res_sensor_projection"," ; x_{DUT}^{pred}-x_{REF}^{pred} [mm];y_{DUT}^{pred}-y_{REF}^{pred}",\
                100,-0.04*MM,0.04*MM,100,-3.5*MM,3.5*MM) 
        self.hcharge = { minst.dut_plane: ROOT.TProfile2D("charge_map_dut" ,";x_{DUT}^{pred} [mm];y_{DUT}^{pred} [mm];"\
                        "<charge cluster> [ADC]", 300,-1.1*sxdut,1.1*sxdut,300,-1.1*sydut,1.1*sydut),
                minst.ref_plane: ROOT.TProfile2D("charge_map_ref",";x_{REF}^{pred} [mm];y_{REF}^{pred} [mm];"\
                        "<charge cluster> [ADC]", 300,-1.1*sxref,1.1*sxref,300,-1.1*syref,1.1*syref) }
        self.hhitmap = { minst.dut_plane: ROOT.TH2F("hitmap_dut" ,";x_{DUT}^{pred} [mm];y_{DUT}^{pred} [mm];"\
                        "Entries", 300,-1.1*sxdut,1.1*sxdut,300,-1.1*sydut,1.1*sydut),
                minst.ref_plane: ROOT.TH2F("hitmap_ref",";x_{REF}^{pred} [mm];y_{REF}^{pred} [mm];"\
                        "Entries", 300,-1.1*sxref,1.1*sxref,300,-1.1*syref,1.1*syref) }
        # -- Module (2 x pitch X)
        self.hcluster_size_mod = { minst.dut_plane: ROOT.TProfile2D("cluster_size_mod_dut" ,";mod(x_{trk})_{2xpitch} [#mum];mod(y_{trk})_{2xpitch}"\
                        "[#mum];<cluster size>", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TProfile2D("cluster_size_mod_ref" ,";mod(x_{trk})_{2xpitch} [mm];mod(y_{trk})_{2xpitch} [mm];"\
                        "<cluster size>", 40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }
        self.hcharge_mod = { minst.dut_plane: ROOT.TProfile2D("charge_mod_dut",";mod(x_{trk})_{2xpitch} [#mum]; mod(y_{trk})_{2xpitch} [#mum];"\
                        "<cluster charge> [ADC]", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TProfile2D("charge_mod_ref",";mod(x_{trk})_{2xpitch} [#mum]; mod(y_{trk})_{2xpitch} [#mum];"\
                        "<cluster charge> [ADC]", 40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }
        self.hhitmap_mod = { minst.dut_plane: ROOT.TH2F("hitmap_mod_dut",";mod(x_{trk})_{2xpitch} [#mum]; mod(y_{trk})_{2xpitch} [#mum];"\
                        "Entries", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM),
                minst.ref_plane: ROOT.TH2F("hitmap_mod_ref",";mod(x_{trk})_{2xpitch} [#mum]; mod(y_{trk})_{2xpitch} [#mum];"\
                        "Entries", 40,0,2.0*self.pitchX[minst.ref_plane]*UM,40,0.0,2.0*self.pitchY[minst.ref_plane]*UM) }
        
        # -- Diagnostics
        self.hhit_raw = { minst.dut_plane: ROOT.TH1F("hitraw_dut",";x_{DUT} [mm];Clusters;",\
                    300,-1.1*sxdut,1.1*sxdut),
                minst.ref_plane: ROOT.TH1F("hitraw_ref",";x_{REF}[mm];Clusters",\
                    300,-1.1*sxref,1.1*sxref) }
        self.residual_projection = { minst.dut_plane: ROOT.TH2F("res_projection_dut",";x_{DUT} [mm];x_{DUT}-x_{trk}^{pred} [mm]",\
                    200,-sxdut,sxdut,200,-3.5*MM,3.5*MM),
                minst.ref_plane: ROOT.TH2F("res_projection_ref",";x_{REF} [mm];x_{REF}-x_{trk}^{pred} [mm]",\
                    200,-sxref,sxref,200,-3.5*MM,3.5*MM) }
        self.trk_iso = { minst.dut_plane: ROOT.TH1F("trkiso_dut","Distance between pair of tracks in the "\
                "same trigger-event;Track distance [mm];Triggers", 400,0,20.0*MM),\
                minst.ref_plane: ROOT.TH1F("trkiso_ref","Distance between pair of tracks in the same "\
                    "trigger-event;Track distance [mm];Triggers", 400,0,20.0*MM) }
        self.nhits_ntrks = { minst.dut_plane: ROOT.TH2F("nhits_ntracks_dut","at least 1 hit matched-isolated;N_{hits};N_{tracks};Triggers",10,-0.5,9.5,40,-0.5,39.5),
                minst.ref_plane: ROOT.TH2F("nhits_ntracks_ref","At least 1 hit matched-isolated ;N_{hits};N_{tracks};Triggers",10,-0.5,9.5,40,-0.5,39.5) }
        self.nhits_ntrks_all = { minst.dut_plane: ROOT.TH2F("nhits_ntracks_all_dut",";N_{hits};N_{tracks};Triggers",10,-0.5,9.5,40,-0.5,39.5),
                minst.ref_plane: ROOT.TH2F("nhits_ntracks_all_ref",";N_{hits};N_{tracks};Triggers",10,-0.5,9.5,40,-0.5,39.5) }
        self.ntracks = ROOT.TH1F("ntracks",";N_{tracks};Triggers",40,-0.5,39.5)
        self.nmatched = ROOT.TH2F("nmatched","Number of isolated track-matched hits per trigger/event;N_{DUT};N_{REF};Triggers",9,-0.5,8.5,\
                9,-0.5,8.5)
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
        self.dutref_distance = ROOT.TH1F("dutref_distance","Distance between tracks matched with DUT and REF;"\
                "x^{pred,DUT}_{DUT}-x^{pred,REF}_{DUT} [mm];a.u.",200,-6.5,6.5)

        diagnostics = self.hhit_raw.values()+self.residual_projection.values()+self.trk_iso.values()+self.nhits_ntrks.values()+\
                self.nhits_ntrks_all.values()+[self.ntracks,self.nmatched]+\
                self.hmatch_eff.values()+self.hy_eff.values()+self.hntrk_eff.values()+self.hnisotrk_eff.values()+\
                [self.dutref_match_eff,self.dutref_pure_match_eff,self.dutref_distance]

        # Correlations: 
        # -- Sensors-tracks
        self.hcorr_trkX = { minst.dut_plane: ROOT.TH2F("corr_trkX_dut",";x_{DUT} [mm]; x_{pred}^{trk} [mm]; Entries",200,-sxdut,sxdut,200,-sxdut,sxdut),
                minst.ref_plane: ROOT.TH2F("corr_trkX_ref",";x_{REF} [mm]; x_{pred}^{trk} [mm]; Entries",200,-sxref,sxref,200,-sxref,sxref)}
        ## Add it as alignment plot
        self._alignment_histos += self.hcorr_trkX.values()

        # -- DUT-REF
        self.hcorrX = ROOT.TH2F("corrX_dut_ref",";x_{DUT} [mm]; x_{REF} [mm]; Entries",100,-sxdut,sxdut,100,-sxref,sxref)
        self.hcorrY = ROOT.TH2F("corrY_dut_ref",";y_{DUT} [mm]; y_{REF} [mm]; Entries",100,-sydut,sydut,100,-syref,syref)
        
        self.hcharge_mod_m = ROOT.TProfile2D("charge_mod_dut_matched",";mod(x_{trk})_{2xpitch} [mm]; mod(y_{trk})_{2xpitch} [#mum];"\
                        "<cluster charge> [ADC]", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM)
        self.hhitmap_mod_m = ROOT.TH2F("hitmap_mod_dut_matched",";mod(x_trk})_{2xpitch} [mm]; mod(y_{trk})_{2xpitch} [#mum];"\
                        "Entries", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM)
        self.hcluster_size_mod_m = ROOT.TProfile2D("cluster_size_mod_dut_m" ,";mod(x_{trk})_{2xpitch} [mm];mod(y_{trk})_{2xpitch} [#mum];"\
                        "<cluster charge> [ADC]", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM)
        # efficiency
        self.heff = ROOT.TProfile2D("eff_map",";x_{trk}^{DUT} [mm]; y^{DUT}_{trk} [mm];#varepsilon",50,-1.1*sxdut,1.1*sxdut,50,-1.1*sydut,1.1*sydut)
        self.heff_entries = ROOT.TH2F("eff_entries",";x_{trk}^{DUT} [mm]; y^{DUT}_{trk} [mm];#varepsilon",50,-1.1*sxdut,1.1*sxdut,50,-1.1*sydut,1.1*sydut)
        self.heff_mod = ROOT.TProfile2D("eff_mod",";mod(x_{trk})_{2xpitch} [#mum]; mod(y_{trk})_{2xpitch} [#mum];"\
                        "efficiency", 40,0,2.0*self.pitchX[minst.dut_plane]*UM,40,0.0,2.0*self.pitchY[minst.dut_plane]*UM,-1,2)
        dut_eff = [self.heff,self.heff_mod,self.heff_entries]

        # -- Eta distribution for matched, isolated
        self.heta = ROOT.TH1F("eta_2","#eta for isolated,matched cluster size 2;#eta;#frac{dN}{d#eta}",100,-0.1,1.1)
        self.heta_g = ROOT.TH1F("eta_g","#eta for isolated,matched cluster size #geq 2;#eta;#frac{dN}{d#eta}",100,-0.1,1.1)
        self.heta_csize = ROOT.TH2F("eta_cluster_size_dut","#eta vs. cluster size for isolated,matched;#eta;N_{cluster};Entries",\
                100,-0.1,1.1, 10,-0.5,9.5)
        self.hcl_size = { minst.dut_plane: ROOT.TH1F("cluster_size_dut","Isolated, REF-matched hits;N_{cluster};Entries",10,-0.5,9.5),
                minst.ref_plane: ROOT.TH1F("cluster_size_ref","Isolated-matched hits;N_{cluster};Entries",10,-0.5,9.5) }
        dut_matched_iso = [self.hcharge_mod_m,self.hhitmap_mod_m,self.hcluster_size_mod_m,\
                self.heta,self.heta_g,self.heta_csize]+self.hcl_size.values()
        # -- Extra histos
        self.htrks_at_planes = { minst.dut_plane: ROOT.TH2F("trk_at_dut","Tracks at DUT;x_{DUT}^{trk} [mm]; y_{DUT}^{trk} [mm]; Entries",\
                        200,-10.0,10.0,200,-10.0,10.0),\
                minst.ref_plane: ROOT.TH2F("trk_at_ref","Tracks at REF;x_{REF}^{trk} [mm]; y_{REF}^{trk} [mm]; Entries",200,-10.0,10.0,200,-10.0,10.0)}
        self.evt_corr = { minst.dut_plane: ROOT.TProfile("dx_correlation_dut","Event correlation closest track;Event;(#Delta_x)_{DUT}",\
                        300000,0.0,299999),\
                    minst.ref_plane: ROOT.TProfile("dx_correlation_ref","Event correlation closest track;Event;(#Delta_x)_{REF}",\
                        300000,0.0,299999)}
        extra = self.htrks_at_planes.values()+self.evt_corr.values()

        # Keep track of all histograms which should be stored in the file
        self._allhistograms = [self.residual_sensor]+\
                self.hcharge.values()+self.hhitmap.values()+\
                self.hcluster_size_mod.values()+self.hcharge_mod.values()+self.hhitmap_mod.values()+\
                self._alignment_histos+\
                [self.hcorrX,self.hcorrY]+\
                dut_matched_iso+dut_eff+diagnostics+\
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
            ### an alignment instance
            ##new_align = alignment_cts()
            ### Update the counter
            ##self.alignment[pl_id].iteration += 1
            ##new_align.iteration=self.alignment[pl_id].iteration
            ### check the new offset
            ##new_x_offset = self.alignment[pl_id].x_offset
            ### Some checks first: do a gross alignment first (before the shift)
            ##if self.dx_h[pl_id].GetEntries() > MIN_ENTRIES:
            ##    new_x_offset = get_x_offset(self.dx_h[pl_id])
            ##    # Finer alignment -- including the shift, therefore to be added
            ##    # to the old offset
            ##    if self.alignment[pl_id].iteration > 1 \
            ##            and abs(new_x_offset-self.alignment[pl_id].x_offset) < 0.1:
            ##        new_x_offset = self.alignment[pl_id].x_offset+get_x_offset(self.dx_finer_h[pl_id],-0.1,0.1)
            ##    # -- The x-offset 
            ##    new_align.x_offset = new_x_offset
            ##else:
            ##    new_align.x_offset = self.alignment[pl_id].x_offset
            ### -- The y-offset (We don't have y-coordinate so far)
            ##new_align.y_offset = 0.0
            ### -- The rotation
            ### XXX --- Maybe if the change is less than 0.1 degrees (or some % of the previous one)
            ###         There is no more alignment in here
            ##if self.dx_y_h[pl_id].GetEntries() > MIN_ENTRIES:
            ##    new_align_rot = self.alignment[pl_id].rot-get_linear_fit(self.dx_y_h[pl_id])
            ##    # Check convergence, if the change are below 0.5 degree, just 
            ##    # converged (not checked in the first iteration)
            ##    if self.alignment[pl_id].iteration > 1 \
            ##            and abs(self.alignment[pl_id].rot-new_align_rot) < 0.0087:
            ##        new_align.rot = self.alignment[pl_id].rot
            ##    else:
            ##        new_align.rot = new_align_rot
            ##else:
            ##    new_align.rot = self.alignment[pl_id].rot
            ### -- The tilt: we don't have coordinate to deal with, just using whatever 
            ###    the user introduced
            ### self.dx_ytx_h --->
            ##if self.dx_ytx_h[pl_id].GetEntries() > MIN_ENTRIES:
            ##    # Slope: microns/radian --> 
            ##    new_tilt = get_linear_fit(self.dx_ytx_h[pl_id],xmin=-0.6,xmax=0.6)*1e3
            ##    # GUARD for stability problems
            ##    if abs(new_tilt) > 1e-3:
            ##        new_align.tilt = 0.0
            ##    else:
            ##        #new_align.turn = self.alignment[pl_id].turn+new_turn/sin(self.alignment[pl_id].turn)
            ##        new_align.tilt = self.alignment[pl_id].tilt-new_tilt
            ###new_align.tilt = self.alignment[pl_id].tilt
            ### -- The turn (only evaluate it if greater than 1 degrees: 0.018 rad.
            ###if abs(self.alignment[pl_id].turn) > 0.02 \
            ###        and self.dx_xtx_h[pl_id].GetEntries() > MIN_ENTRIES:
            ##if self.dx_xtx_h[pl_id].GetEntries() > MIN_ENTRIES:
            ##    # In mmrad
            ##    new_turn = get_linear_fit(self.dx_xtx_h[pl_id],xmin=-0.6,xmax=0.6)*1e3
            ##    # GUARD for stability problems
            ##    if abs(new_turn) > 1e-3:
            ##        new_align.turn = 0.0
            ##    else:
            ##        #new_align.turn = self.alignment[pl_id].turn+new_turn/sin(self.alignment[pl_id].turn)
            ##        new_align.turn = self.alignment[pl_id].turn+new_turn
            ####else:
            ##    # Compatible with zero (so, the resolution we have in this angle is 1 degree)
            ##    # Going down, it could drive to unsense results (because of the division)
            ####    new_align.turn = 0.0
            ### -- the z-shift: Note the SPS beam hardly divergent (less than 0.1 mrad),
            ###    therefore do not trust too much in this correction, not stable. 
            ###    We trust in the initial measurement, if the diference is higher than a couple of mm, 
            ###    deactivate it
            ##if self.alignment[pl_id].iteration > 1 \
            ##        and self.dx_tx_h[pl_id].GetEntries() > MIN_ENTRIES:
            ##    # Slope in um/mm
            ##    new_dz = get_linear_fit(self.dx_tx_h[pl_id],xmin=-0.1,xmax=0.1,robval=0.6)*1e3
            ##    # -- XXX Guard to avoid not estable results: changes in 10 mm makes no sense
            ##    #        Note that the fit turns out not to be too stable...
            ##    if abs(new_dz) > 100.0:
            ##        new_align.dz = self.alignment[pl_id].dz
            ##    else:
            ##        #new_align.dz = self.alignment[pl_id].dz+new_dz
            ##        new_align.dz = self.alignment[pl_id].dz-new_dz
            ##else:
            ##    new_align.dz = self.alignment[pl_id].dz
            #### Update the module
            ###filename = ALIGN_FILE_FORMAT.format(pl_id,self.run_number)
            ###new_align.update(filename,self.run_number)
            #### And prepare the alignment constants to be updated 
            #### in the next iteration
            ###hits_plane_accessor.resync[pl_id] = True
            
            params = ROOT.TVectorD()
            errors = ROOT.TVectorD()
            self.dx_hyp[pl_id].Eval()
            #self.dx_hyp[pl_id].EvalRobust(0.7)
            self.dx_hyp[pl_id].GetParameters(params)
            #self.dx_hyp[pl_id].GetErrors(errors)

            # an alignment instance
            new_align = alignment_cts()

            for (argi,(name,pre)) in new_align.hyp_par:
                v = pre*params[argi]
                # Finer alignment, check after the x_offset is evaluated
                # if the correction is small enough
                if argi == 0 and \
                        (abs(new_align.x_offset) > 0.1 or \
                        self.dx_h[pl_id].GetEntries() < MIN_ENTRIES):
                    new_align.x_offset = v
                    break
                if name not in ["x_offset", "y_offset"]:
                    print "="*80
                    print "\033[1;33mWARNING\033[1;m Alignment is de-activated"\
                            " (except for the x_offset). Still  some  issues"
                    print "\033[1;33mWARNING\033[1;m are affecting the alignment "\
                            "alogorithm. Work in progress..."
                    print "="*80
                    # For instance: Changing the dz estimation makes no difference
                    #               almost in the dx_y_h histogram, this is the 
                    #               symptom of something weird going on, however
                    #               not sure how to debug it
                    continue
                # Either turn or tilt only evaluated if they are higher than
                # 5 deg.
                # XXX- The user should be able to include the angle 
                if name in ["tilt","turn"] and \
                        abs(getattr(self.alignment[pl_id],name)) < 0.088:
                    continue
                # If the difference is too big, ignore it
                #if abs(v/getattr(self.alignment[pl_id],name)) > 0.5:
                #    continue
                setattr(new_align,name,v)
                # Check the size of the errors
                #if errors
                #print " - extra {0:10}={1:.5f} [{2}] (error:{3})".format(name,\
                #        v*new_align.units[name][0],\
                #        new_align.units[name][1],\
                #        new_align.units[name][0]*errors[argi])
                #print "   --- t-value:{0}, significance:{1}".format(self.dx_hyp[pl_id].GetParTValue(argi),self.dx_hyp[pl_id].GetParSignificance(argi))
            print
            print "Extra alignment (PL:{0})".format(pl_id)
            print new_align
            # Update the alignment, correcting the old constants
            self.alignment[pl_id] += new_align
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
        self.fill_statistics(duthits.n,refhits.n,trks.n)

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

        # --- XX
        self.nhits_ntrks_all[refhits.id].Fill(refhits.n,trks.n)
        self.nhits_ntrks_all[duthits.id].Fill(duthits.n,trks.n)
        
        # Prepare the ntuple of histograms for DUT and REF
        histos_dut = (self.hcorr_trkX[duthits.id],self.residual_projection[duthits.id],self.dx_h[duthits.id],\
                self.dx_finer_h[duthits.id],self.dx_y_h[duthits.id],self.dx_ytx_h[duthits.id],\
                self.dx_xtx_h[duthits.id],self.dx_tx_h[duthits.id],self.hplane[duthits.id],\
                self.hmatch_eff[duthits.id],self.hy_eff[duthits.id],self.hntrk_eff[duthits.id],self.hnisotrk_eff[duthits.id],\
                self.dx_hyp[duthits.id])
        histos_ref = (self.hcorr_trkX[refhits.id],self.residual_projection[refhits.id],self.dx_h[refhits.id],\
                self.dx_finer_h[refhits.id],self.dx_y_h[refhits.id],self.dx_ytx_h[refhits.id],\
                self.dx_xtx_h[refhits.id],self.dx_tx_h[refhits.id],self.hplane[refhits.id],\
                self.hmatch_eff[refhits.id],self.hy_eff[refhits.id],self.hntrk_eff[refhits.id],self.hnisotrk_eff[refhits.id],\
                self.dx_hyp[refhits.id])
        histos = { refhits.id: histos_ref, duthits.id: histos_dut }
        # Get the matched hits { sensorID: [ (hit_index,track_index),. ..] , }
        #matched_hits = trks.matched_hits(duthits,refhits,self.hcorr_trkX)
        matched_hits = { duthits.id: duthits.match_isolated(trks,histos_dut),
                refhits.id: refhits.match_isolated(trks,histos_ref,True) }

        # Fill some track histograms
        for itr in xrange(trks.n):
            for ho in sensor_hits.values():
                trks.fill_isolation_histograms(itr,ho,self.trk_iso[ho.id])
                ((xpred,ypred,zpred),tel) = trks.get_point_in_sensor_frame(itr,ho)
                self.htrks_at_planes[ho.id].Fill(xpred,ypred)
        # Some diagnostic for sensor hits
        for plane,oh in sensor_hits.iteritems():
            for ih in xrange(oh.n):
                self.hhit_raw[plane].Fill(oh.x_local[ih])
        # --
        # [H]--for pl,plist in matched_hits.iteritems():
        # [H]--    print "PLANE",pl,"CHANNEL"
        # [H]--    for hj,ht in plist:
        # [H]--        print self.sizeX[pl],self.pitchX[pl]
        # [H]--        print sensor_hits[pl].x_local[hj],sensor_hits[pl].x_channel(hj,self.sizeX[pl],self.pitchX[pl])

        
        # Some diagnostic plots
        # -- number of tracks in the event
        self.ntracks.Fill(trks.n)
        # -- number of matched,isolated hits in the event
        self.nmatched.Fill(len(matched_hits[duthits.id]),len(matched_hits[refhits.id]))
        
        # If alignment, just return here, trying to avoid extra processing time
        if is_alignment:
            # Update the alignment constants 
            self.alignment = dict(map(lambda (i,d): (i,d.align_constants[i]) ,sensor_hits.iteritems()))
            # Get an estimation of the sensor efficiency ?
            self.fill_statistics_matched(matched_hits)
            return
        # Filling histograms
        for sensorID,pointlist in matched_hits.iteritems():
            hits = sensor_hits[sensorID]
            for hit_el,trk_index in pointlist:
                # Remember: it_el=index of the measured hit at the corrent sensor (sensorID)
                #           trk_index= index of the track
                ((xpred,ypred,zpred),rtel) = trks.get_point_in_sensor_frame(trk_index,hits)
                # XXX TO CODE: Fiducial cut: some percentange from the edges
                # -- track-hits correlation
                self.evt_corr[sensorID].Fill(self.total_events,hits.x_local[hit_el]-xpred)
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
                self.nhits_ntrks[sensorID].Fill(hits.n,trks.n)
                ## INCLUDE IT IN HERE
        # Matching beetwen DUT and REF: hits with the same track ID
        # Loop over all the measured hits in the REF, matching to a
        # track, using the track as common element in the DUT-REF match
        # XXX Just using those within the fiducial cut of the DUT XXX
        # First check the ref matching efficiency
        for ihit_dut,itrk_d in matched_hits[duthits.id]:
            ref_match = filter(lambda (ir,itrk_r): itrk_d == itrk_r,matched_hits[refhits.id])
            self.dutref_match_eff.Fill(duthits.x_local[ihit_dut],(len(ref_match) > 0))
            # Track efficiencies incorporated
            if len(matched_hits[refhits.id]) > 0:
                self.dutref_pure_match_eff.Fill(duthits.x_local[ihit_dut],(len(ref_match) > 0))
            # Whats the position the track matched with the REF gives in the DUT?
            for ihit_ref,itrk_r in matched_hits[refhits.id]:
                self.dutref_distance.Fill(trks.get_point_in_sensor_frame(itrk_d,duthits)[0][0]-\
                        trks.get_point_in_sensor_frame(itrk_r,refhits)[0][0])
        for ihit_ref,itrk in matched_hits[refhits.id]:
            # Cluster size for track matched ref hits
            self.hcl_size[refhits.id].Fill(refhits.n_cluster[ihit_ref])
            # Get the track predicted points at the DUT plane
            ((xpred_dut,ypred_dut,zpred_dut),rtel) = trks.get_point_in_sensor_frame(itrk,duthits)
            # Only within fiducial  (note that we don't have coordinate Y)
            if xpred_dut > duthits.xmax or xpred_dut < duthits.xmin:
                continue
            # Modulo pitch
            xmod = (32.0+xpred_dut)%(2.0*self.pitchX[duthits.id])
            ymod = (32.0+ypred_dut)%(2.0*self.pitchY[duthits.id])
            # For the given track (matched with the ref), see if the track
            # matched as well with the DUT
            dut_match = filter(lambda (ihit_dut,itrk_dut): itrk_dut == itrk,matched_hits[duthits.id])
            any_match = int(len(dut_match)>0)
            # Efficiency
            self.heff_mod.Fill(xmod*UM,ymod*UM,any_match)
            self.heff.Fill(xpred_dut,ypred_dut,any_match)
            self.heff_entries.Fill(xpred_dut,ypred_dut)
            if len(dut_match) > 0:
                # Get the track predicted points at the REF plane
                ((xpred_at_ref,ypred_at_ref,zpred_at_ref),rptel) = trks.get_point_in_sensor_frame(itrk,refhits)
                # And some correlation plots of the predicted values between
                # the DUT and REF sensors
                # Note that the expected gaussian could be not centered at zero,
                # in that case the difference will be related with the misaligned 
                # between them
                self.hcorrX.Fill(xpred_dut,xpred_at_ref)
                # Not measured elements at Y, so using the predictions
                self.hcorrY.Fill(ypred_dut,ypred_at_ref)
                ## Residuals only for the DUT-REF matched tracks, i.e gives the 
                ## relative position between ref and dut
                ## XXX
                self.residual_sensor.Fill(xpred_dut-xpred_at_ref,ypred_dut-xpred_at_ref)
                # Cluster size 
                self.hcl_size[duthits.id].Fill(duthits.n_cluster[duthits.track_link[itrk]])
                self.hcluster_size_mod_m.Fill(xmod*UM,ymod*UM,duthits.n_cluster[duthits.track_link[itrk]])
                # -- charge
                self.hcharge_mod_m.Fill(xmod*UM,ymod*UM,hits.charge[hit_el])
                self.hhitmap_mod_m.Fill(xmod*UM,ymod*UM)
                # -- Eta distribution for the matched cluster size=2
                clustersize = duthits.n_cluster[duthits.track_link[itrk]]
                if clustersize >= 2:
                    dhit = duthits.track_link[itrk]
                    self.heta_g.Fill(duthits.eta[dhit])
                    self.heta_csize.Fill(duthits.eta[dhit],clustersize)
                    self.duthits_matched_eta_h.append((duthits.eta[dhit],duthits.charge[dhit]))
                    if clustersize == 2:
                        self.heta.Fill(duthits.eta[dhit])
                        self.duthits_matched_eta.append(self.duthits_matched_eta_h[-1])
            elif len(dut_match) > 1:
                raise RuntimeError("This must never happens!! Contact developer, error 3E43")
        self.fill_statistics_matched(matched_hits)

    def get_raw_sensors_efficiency(self):
        """Summarize the efficiency of the sensors:
        eff = Number of matched hits/number of events with tracks
        """
        m = ""
        for plid,n in self.number_matched_hits.iteritems():
            if plid == self.dut_plane:
                sensor_evts = self.events_with_dut
                sensor_name = "DUT"
            else:
                sensor_evts = self.events_with_ref
                sensor_name = "REF"
            m+= "{0} efficiency (Isolated track-matched): {1:.2f}%\n".format(sensor_name,float(n)/float(self.total_events_tracks)*100.)
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

def get_x_offset(h,xmin=-2.0,xmax=2.0):
    """XXX 
    """
    import ROOT
    ROOT.gROOT.SetBatch()
    
    cns=ROOT.TCanvas()

    if h.GetEntries() == 0:
        raise RuntimeError("Empty histogram '{0}'".format(h.GetName()))

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

    ##gbg = ROOT.TF1("linear_fit_{0}".format(h.GetName()),"pol1",xmin,xmax)
    ### Convert to a TGraphErrors to perform the robust fit
    ##x    = array.array('f',map(lambda i: h.GetBinCenter(i),xrange(1,h.GetNbinsX()+1)))
    ##xerr = array.array('f',map(lambda i: 0.0,xrange(1,h.GetNbinsX()+1)))
    ##y    = array.array('f',map(lambda i: h.GetBinContent(i),xrange(1,h.GetNbinsX()+1)))
    ##yerr = array.array('f',map(lambda i: h.GetBinError(i),xrange(1,h.GetNbinsX()+1)))

    ##gr = ROOT.TGraphErrors(len(x),x,y,xerr,yerr)
    ###  Do the fit: Note that is using a robust fit (assuming 
    ###  the 75% of the points at least are not outliers -- I need a TGraph to do that
    ##status = gr.Fit(gbg,"QS+ROB={0}".format(robval))
    ##if status.Status() != 0:
    ##    print "\033[1;33mWARNING\033[1;m FAILED THE LINEAR FIT for '{0}': "\
    ##            "Status {1}".format(h.GetName(),status.Status())
    ### Add the function to the histogram  XXX NOT WORKING -- at some point 
    ### the functions are deleted
    ###h.GetListOfFunctions().Add(gbg)
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

    global DEBUG
    DEBUG=verbose


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
    # Note, loose cut for the dut, while tight cut for the ref link
    RESOLUTION= { metadata.dut_plane: 6.0*specs[name_converter[fp.sensor_name]].pitchX,\
            metadata.ref_plane : 2.0*specs["REF_0_b1"].pitchX } # in mm
    # 

    # Set some globals
    #tree_inspector(t)

    # The hits and tracks 
    dut = hits_plane_accessor(t,metadata.dut_plane)
    ref = hits_plane_accessor(t,metadata.ref_plane)
    tracks = tracks_accessor(t,[0,1,2,3,4],metadata.ref_plane,metadata.dut_plane)
    
    # Process the data
    #start_time = timeit.default_timer()
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
    #print timeit.default_timer()-start_time
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
