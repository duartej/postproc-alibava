#!/usr/bin/env python
"""Post-processor module to be applied to 
the output of the EUTelescope Marlin jobs
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
import math
# TODO OBTAIN FROM SPS
REFPLANE=7
DUTPLANE=6
RESOLUTION= { DUTPLANE: 0.2, REFPLANE : 0.1 } # in mm
SENSOR_NAME = { REFPLANE: "REF", DUTPLANE: "DUT" }

# UNITS
mm = 1.0
um = 1e3

class metadata_container(object):
    """Container of some globals which are going to be
    used along the module; these globals are run/file
    dependent.
    """
    def __init__(self,tree,dut_name):
        """Use the tree structure to define the DUT 
        planes ids, their resolutions, etc

        Raises
        ------
        RuntimeError:
            When the tree doesn't contain the Branch hit_X_#
        """
        from .SPS2017TB_metadata import sensor_name_spec_map
        import math

        # The REF plane is explicitely set to 7 ALWAYS
        self.ref_plane = 7
        # Get the DUT by using the hit branches
        try:
            self.dut_plane = map(lambda br: int(br.GetName().replace("hit_X","")),
                filter(lambda x: x.GetName().find("hit_X_") == 0, tree.GetListOfBranches()))
        except IndexError:
            raise RuntimeError("Unexpected tree structure")
        # Resolutions: REF and DUT (binary plust
        _p = 1.0/sqrt(12.0)        
        self.resolution = { self.ref_plane: sensor_name_spec_map["REF_0_b1"].pitchX*_p, 
                self.dut_plane: sensor_name_spec_map[dut_name]*_p }

class dummy_list(list):
    """A dummy list which returns always
    the value in the firsts element
    """
    def __getitem__(self,i):
        return list.__getitem__(self,0)

class hits_plane_accessor(object):
    """Access to the TBranches contained in a 
    TTree created by the EUTelTreeCreator class from
    the https://github.com/duartej/eutelescope/ package.
    """
    def __init__(self,tree,planeid):
        """Accessor
        """
        # Improve!! Associate just a TBranch
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
        return self.x.size()
    
    def equal(self,i,other,j):
        """Allow comparation beetween hits, we assume equal
        if (x,z) are equal
        """
        return abs(self.x[i]-other.x[j]) < 1e-9 \
                and abs(self.z[i]-other.z[j]) < 1e-9


    def matched_hits(self,track_acc,histo):
        """Search for tracks matching (within a resolution) the list 
        of hits in the current event

        Parameters
        ----------
        track_acc: 

        Return
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
    """
    def __init__(self,tree,planeid,ismeas):
        """
        """
        if ismeas:
            prefix = "trk_hit_meas"
        else:
            prefix = "trk_hit_fit"
        # Improve!! Associate just a TBranch
        self.id = planeid
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
        return self.x.size()
    
    def equal(self,i,other,j):
        """Allow comparation beetween hits, we assume equal
        if (x,z) are equal
        """
        return abs(self.x[i]-other.x[j]) < 1e-9 \
                and abs(self.z[i]-other.z[j]) < 1e-9

    def get_hits_with_same_track_indices(self,itrk):
        """Return the list of indices for those hits belonging 
        to the itrk-track
        """
        return filter(lambda i: self.track_index[i] == itrk,xrange(self.n))

    def get_hit_with_same_track_index(self,i,other):
        """Return the index corresponding to 'other'
        that has the same 'track_index' data member that the i-hit

        Parameters
        ----------
        i: int
            The index of the hit to check
        other: track_hits_plane_accessor
            The instance to check for a hit with the same track 
            index than the i-hit
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
        """
        """
        return self.id == DUTPLANE
    
    def is_ref(self):
        """
        """
        return self.id == REFPLANE
    

class tracks_accessor(object):
    """Access to the TBranches contained in a 
    TTree created by the EUTelTreeCreator class from
    the https://github.com/duartej/eutelescope/ package.
    """
    def __init__(self,tree,tel_planes,ref_plane,dut_plane):
        """
        Parameters
        ----------
         hits
        """
        # Improve!! Associate just a TBranch
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
        self.telescope_measured_hits = map(lambda i: track_hits_plane_accessor(tree,i,True),self.tel_planes)
        self.telescope_fitted_hits = map(lambda i: track_hits_plane_accessor(tree,i,False),self.tel_planes)
        # --- dut
        self.dut_measured_hits = track_hits_plane_accessor(tree,self.dutid,True)
        self.dut_fitted_hits = track_hits_plane_accessor(tree,self.dutid,False)
        # --- ref
        self.ref_measured_hits = track_hits_plane_accessor(tree,self.refid,True)
        self.ref_fitted_hits = track_hits_plane_accessor(tree,self.refid,False)
        
    @property
    def n(self):
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
        """
        t = z-self.z0[i]
        return (self.x0[i]+self.dxdz[i]*t,self.y0[i]+self.dydz[i]*t)

class processor(object):
    def __init__(self):
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
        self.residual_projection = { DUTPLANE: ROOT.TH1F("res_projection_dut"," ; x_{DUT}-x_{trk}^{pred} [mm]; Entries",100,-RESOLUTION[DUTPLANE],RESOLUTION[DUTPLANE]),
                REFPLANE: ROOT.TH1F("res_projection_ref"," ; _x_{REF}-x_{trk}^{pred} [mm]; Entries",100,-RESOLUTION[REFPLANE],RESOLUTION[REFPLANE]) }
        # -- histograms DUTS/REF
        self.hmap = { DUTPLANE: ROOT.TH2F("local_map_dut" ,";x [mm]; y [mm]; Entries", 130,-10.0,10.0,100,-10.0,10.0),
                REFPLANE: ROOT.TH2F("local_map_ref",";x [mm]; y [mm]; Entries", 130,-10.0,10.0,100,-10.0,10.0) }
        self.hcharge = { DUTPLANE: ROOT.TH2F("charge_map_dut",";x [mm]; y [mm]; charge cluster [ADC]",130,-10.0,10.0,100,-10.0,10.0),
                REFPLANE: ROOT.TH2F("charge_map_ref",";x [mm]; y [mm]; charge cluster [ADC]",130,-10.0,10.0,100,-10.0,10.0) }
        # Correlations
        self.hcorr_trkX = { DUTPLANE: ROOT.TH2F("corr_trkX_dut",";x_{DUT} [mm]; x_{pred}^{trk} [mm]; Entries",200,-10.0,10.0,200,-10.0,10.0),
                REFPLANE: ROOT.TH2F("corr_trkX_ref",";x_{REF} [mm]; x_{pred}^{trk} [mm]; Entries",200,-10.0,10.0,200,-10.0,10.0)}
        self.hcorrX = ROOT.TH2F("corrX_dut_ref",";x_{DUT} [mm]; x_{REF} [mm]; Entries",100,-10.0,10.0,100,-10.0,10.0)
        self.hcorrY = ROOT.TH2F("corrY_dut_ref",";y_{DUT} [mm]; y_{REF} [mm]; Entries",100,-10.0,10.0,100,-10.0,10.0)
        # efficiendy
        self.hpass = ROOT.TH2F("pass_map",";x [mm]; y [mm]; #varepsilon",130,-10.0,10.0,100,-10.0,10.0)
        self.htotal = ROOT.TH2F("total_map",";x [mm]; y [mm]; #varepsilon",130,-10.0,10.0,100,-10.0,10.0)

        self._allhistograms = self.residual_projection.values()+self.hmap.values()+\
                self.hcharge.values()+self.hcorr_trkX.values()+\
                [self.hcorrX,self.hcorrY,self.hpass,self.htotal]
        dummy=map(lambda h: h.SetDirectory(0),self._allhistograms)

    def draw_all_in_canvas():
        """
        """
        import ROOT

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
        m = "Events with DUT: {0} (eff: {1:.2f}%)\n".format(self.events_with_dut,float(self.events_with_dut)/float(self.total_events)*100.)
        m+= "Events with REF: {0} (eff: {1:.2f}%)\n".format(self.events_with_ref,float(self.events_with_ref)/float(self.total_events)*100.)
        m+= "\nEvents with DUT but no tracks: {0}\n".format(self.events_with_dut_no_tracks)
        m+= "Events with REF but no tracks: {0}\n".format(self.events_with_ref_no_tracks)
        m+= "Events with REF and DUT and tracks: {0}".format(self.events_with_ref_dut_tracks)

        return m

def run(fname):
    import ROOT
    
    import timeit

    f =ROOT.TFile.Open(fname)
    t = f.Get("events")
    # Set some globals
    tree_inspector(t)

    # The hits and tracks
    dut = hits_plane_accessor(t,DUTPLANE)
    ref = hits_plane_accessor(t,REFPLANE)
    tracks = tracks_accessor(t,[0,1,2,3,4],REFPLANE,DUTPLANE)

    start_time = timeit.default_timer()
    wtch = processor()
    for _t in t:
        wtch.process(_t,tracks,ref,dut)
    print timeit.default_timer()-start_time
    return wtch


if __name__ == "__main__":
    import sys
    run(sys.argv[1])
