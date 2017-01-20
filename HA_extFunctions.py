#!/usr/bin/python
# ==================================================================== #
### some little helpers ###
#
# FK2014: summary of little helpers

#############################################
## See LICENSE for licensing and copyright ##
#############################################

########################################################################

# ==================================================================== #

## Convert phase values to cosine to get better hypal ##
# Conversion functions orignal by Michael Hanke # 
# need for back to deg for further presentation or statistics.
# deg2cos is ONLY for better hypal projection, cause it has problems 
# with circular data
def _deg2cosine(deg, channels=10):
    from mvpa2.datasets import Dataset
    import numpy as np
    s = np.deg2rad(deg.samples[0])
    cosine = np.array([np.cos(s + (i * 2 * np.pi / channels))
                        for i in range(channels)])
    DS = Dataset(cosine, fa=deg.fa, a=deg.a)
    DS.a['degORcosORts'] = "cos"
    return DS

def _shifted_cosine(x, phase_shift):
    import numpy as np
    return np.cos(x + phase_shift)

def _decode_cosine_phase(ds):
    from mvpa2.misc.fx import least_sq_fit
    from mvpa2.datasets import Dataset
    import numpy as np
    from HA_extFunctions import _shifted_cosine 
    
    x = [i * 2 * np.pi / len(ds) for i in range(len(ds))]
    def fit_helper(d):
        fit = least_sq_fit(_shifted_cosine, [0.0], d, x)
        if not fit[1] in [1,2,3,4]:
            #raise RuntimeError("leastsq fit failed")
            print "DBG: leastsq fit failed for a voxel"
            print 'DBG: this is not a cosine:', d
            return 0.0
            # TODO: put me back: np.nan
        return fit[0]
    angles = np.apply_along_axis(fit_helper, 0, ds.samples)
    deg = np.rad2deg(angles).squeeze()
    # make all positive angle -- rotate around if necessary
    deg[deg < 0] = deg[deg < 0] + 360
    return Dataset([deg], fa=ds.fa, a=ds.a)

def _get_avg_euc_dist(maps):
    from scipy.stats import circmean, circvar, circstd
    import numpy as np
    rmss = []
    
    # go through all input maps
    for i, m in enumerate(maps):
        
        # compute the average map while holding out the current one
        avg_map= np.array([circmean(ds,low=0,high=360) for ds in np.array([x.samples for j, x in enumerate(maps) if not j == i]).squeeze().T],ndmin=2)
        
        # compute avg absolute per voxel difference from the avg map
        rms = np.sum(np.abs(m - avg_map)) / avg_map.shape[1]
        rmss.append((rms,
                     avg_map,
                     np.corrcoef(np.cos(np.deg2rad(avg_map)), np.cos(np.deg2rad(m)))[0,1]))
    return rmss

# ==================================================================== #

def _ds2deg(backproj):
    """ convert the ds back to deg """
    
    from HA_extFunctions import _decode_cosine_phase
    from mvpa2.suite import vstack
    back = []
    nscansProj = backproj.shape[0]/backproj.a.deg2cos_channels
    for scan in range(nscansProj):
        # scans are defined by thier nVolumes with start (vv) and 
        # end (bb) Volume #
        vv = (scan * backproj.a.deg2cos_channels)
        bb = vv + backproj.a.deg2cos_channels
        bp_tmp = _decode_cosine_phase(backproj[vv:bb,:])
        print "DEBUG: .. convert volumes of scan",scan,": vv,bb :", vv,bb
        print "DEBUG: .. bp_tmp.shape, bp_tmp.max(), bp_tmp.min(): ", bp_tmp.shape, bp_tmp.S.max(), bp_tmp.S.min()
        back.append(bp_tmp)
    bp = vstack(back)
    # save attributes
    bp.a = backproj.a
    backproj = bp
    backproj.a['degORcosORts'] = "deg"
    return backproj

def _getCosCdist(DS1, DS2, channels=10, doBack2deg=True):
    """ 
    correlate cos values to get Cdist between two datasets phase values 
    doBack2deg : False  - expects _ds2cos transformed datasets
                 True   - convert to cos values to correlate them
    """
    
    import scipy.spatial.distance as sd
    import numpy as np
    from HA_extFunctions import _ds2deg, _circdescribe
    
    # each phase (ECC, POL) is converted into nChannels values
    # i.e. 1 phaseval (eg 45 deg) is converted into 10 cos values
    # TODO: now takes nChannel for both inputDS -> make it that both can be different
    if not 'channels' in locals(): channels = DS1.a.deg2cos_channels
    
    # just ignore if no ds as input (use it for plottools..)
    DS1.a['deg2cos_channels'] = channels
    DS2.a['deg2cos_channels'] = channels

    cosCdist = []
    for mul in [0,1]:
        ds1 = DS1[(mul*channels):(mul+1)*channels]
        ds2 = DS2[(mul*channels):(mul+1)*channels]
        if doBack2deg:
            ## transform the X-channel-cos-values back to degree and correlate the rad-value ##
            # correlate phase values, thus transform back to deg [0 .. 360]
            ds1 =  _ds2deg(ds1)
            ds2 =  _ds2deg(ds2)
            
            # if POL values, use cos to get "circular" stats
            if mul == 1: 
                ds1 = np.cos(np.deg2rad(ds1))
                ds2 = np.cos(np.deg2rad(ds2))
            
            tmp = np.mean(np.diag(sd.cdist(ds1,ds2,'co')))
        else:
            # correlate transformed cos values (nChannels for each ph value)
            # this would implicate circular stats on ECC and POL
            ds1 = ds1.S.reshape(1,ds1.shape[0]*ds1.shape[1])
            ds2 = ds2.S.reshape(1,ds2.shape[0]*ds2.shape[1])
            tmp = sd.cdist(ds1,ds2,'co').squeeze()
        cosCdist.append(tmp)

    return cosCdist

def _circdescribe(arr):
    """
    expect values between 0 and 360 deg (see "low" and "high" boundaries)
     len, (min,max), avg, var
    """
    
    from scipy.stats import circmean, circvar, circstd
    return len(arr),(arr.min(),arr.max()),circmean(arr,low=0,high=360), \
            circvar(arr,low=0,high=360),circstd(arr,low=0,high=360)

def _roitmpl2roiname(roitmpl):
    """
    generate roiNames out of roitempl 
    
    written in an amazing function to easyly change pattern for roiNames
    in all further function at once. This way keep the 
    featureAttributes (fa) consistent
    """
    
    roiRname    = 'Right'+roitmpl
    roiLname    = 'Left'+roitmpl 
    return roiLname,roiRname

def _findfile(path, pattern):
    """ find file depending on pattern """
    
    import os, fnmatch
    result = []
    for root, dirs, files in os.walk(path):
        #DEBUG: print root, dirs, files
        for name in files:
            if fnmatch.fnmatch(name, pattern): 
                result.append(os.path.join(root, name))
    return result

def _dimreducer(inNii, keepdim=None):
    """ remove AFNI-dimenstion (brick ...) """
    
    import nibabel as nb
    import os
    # exclude AFNI sub-brick dimension to make it 
    # readable for mvspa2.suite.fmri_dataset
    fileorg = nb.load(inNii)
    if keepdim is None: 
        datared = fileorg.get_data().squeeze()
    else: 
        datared = fileorg.get_data().squeeze()[:,:,:,keepdim]
    print 'Warning: using .nii as splitter of filename -> could cause ERRORS'
    outname = ''.join([inNii.split('.nii')[0],'_reddim'+str(keepdim)+'.nii.gz'])
    # use header of loaded inNii
    outnifti = nb.Nifti1Image(datared, fileorg.get_affine())
    
    # check if file exist and delete in this case because nb.save can't overwrite
    if os.path.isfile(outname): os.remove(outname) 
    
    nb.save(outnifti, outname)
    return outname

def _getMeanCS(CSexlToSubj,cs_fromSubj):
    """
    get mean map in CS: mean(other), thus excl. toSubj 
     
    all DS need to have the same shape -> works for CS or 
    data in grpSpace with same mask
    """
    
    import numpy as np
    from HA_extFunctions import _ds2deg
    from mvpa2.suite import vstack
    
    # mean in CS
    Ecc_tmp,Pol_tmp = [],[]
    half = CSexlToSubj[0].shape[0]/2
    for csTmp in CSexlToSubj:
        Ecc_tmp.append(csTmp.S[:half,:])
        Pol_tmp.append(csTmp.S[half:,:])
    Pol_mean = np.mean(np.array(Pol_tmp),axis=0)
    Ecc_mean = np.mean(np.array(Ecc_tmp),axis=0)

    # merge Ecc and Pol to be able to use the rest of the scripts
    meanVals = vstack([Ecc_mean,Pol_mean])

    # write mean values into the samples of a CS to get the attributes back
    cs_mean = cs_fromSubj.copy()
    cs_mean.samples[:,:] = meanVals        
    return cs_mean


def _get_seed_mean(seedRegionDS,sparse_radius, nproc, method="spherepacking"):
    """
    generate spheres outside mask for Connectome
    
    do all for different radii of spheres of connectome and for
    all subjects.
    
    (Use the biggest brain as reference to get sphereCenter)
    (scoords, idx)
    >> it is indepenent from ROIvoxel
    Now can choose between method of spherepacking (by FloB) and 
    scatter_neighborhoods by HYPAL 
    """
    
    from mvpa2.misc.neighborhood import Sphere, scatter_neighborhoods
    from mvpa2.measures.searchlight import sphere_searchlight
    from mvpa2.mappers.fx import mean_group_feature
    from spherepack import packing
    import numpy as np
    from mvpa2.suite import fmri_dataset
  
    # get coords to place spheres for searchlight
    # - coordinates may vary between runs
    # - sphere center dont overlap but sphere elements may overlap 
    if method == "spherepacking":
        sidx, scoords = packing(seedRegionDS,radius=sparse_radius)
    elif method == "spheres_overlapping":
        sidx, scoords = scatter_neighborhoods(Sphere(sparse_radius),
                                          seedRegionDS.fa.voxel_indices,
                                          deterministic=True)
    # get mean of tseries within these spheres
    seed_mean = sphere_searchlight(mean_group_feature(attrs=None),
                                                radius=sparse_radius,
                                                nproc=nproc,
                                                center_ids=sidx)

    return seed_mean, sidx, scoords

def _back2nifti(headerDs, nVol, ds, nscansProj, nscansMapper, fromSubj, toSubj,nameadd=""):
    """
    maskpath: take that niftis header to project in this space.
    nVol:      use that to split dataset into each parts. nVol is nVols of 
              scans loaded for dsAll
    fromSubject:    subject to 
    headerDs: use the header of the dataset for the map2nifti reverse-mapping
    """
    
    from mvpa2.suite import map2nifti
    import os
    import nibabel as nb

    niftiout = []
    # number of scans left to project
    for scan in range(nscansProj):
        # scans are defined by thier nVolumes with start (vv) and 
        # end (bb) Volume #
        vv = (scan*nVol)
        bb = vv + nVol
        
        # make a backup of the ds to keep only the Volumes for the actual 
        # scan
        dsback = ds.copy()
        dsback.samples = dsback.samples[vv:bb]
        
        dsback2nifti = map2nifti(headerDs, dsback)
        savekey = \
            'used'+str(nscansMapper)+'Scans_part'+str(scan+1)+ \
            '_SubjFrom'+str(fromSubj)+'to'+str(toSubj)+'.nii.gz'
        dsback2niftisave = 'Backproj'+nameadd+'_'+savekey
        
        # check if file exist and delete in this case because nb.save can't overwrite
        if os.path.isfile(dsback2niftisave): os.remove(dsback2niftisave)
     
        nb.save(dsback2nifti, dsback2niftisave)
        print 'saved: ',dsback2niftisave
        niftiout.append(os.path.abspath(dsback2niftisave))
    return niftiout

def _getIntersectionBrainMask(DSlist):
    """
    uses brainmasks from ds.fa to generate dataset with BRAINmask to save as nifti
    """
    
    from mvpa2.suite import map2nifti, fmri_dataset, h5load
    import numpy as np
    import nibabel as nb
    import os
    import copy

    dsback2niftiPath = 'IntersectionBrainmask.nii.gz'
    
    # get brainmask to start with
    if type(DSlist[0]) == str: dstmp = h5load(DSlist[0])
    else: dstmp = DSlist[0]
    
    def get_BMwithGM(ds):
        bm = fmri_dataset(samples=ds.a.brainmask,
                            add_fa = {  'rightGM': ds.a.rightGM,
                                        'leftGM': ds.a.leftGM
                                     })
        ## get intersection of brainmask and GM (gray matter)
        # mark where GM
        idx = np.logical_or(bm.fa.leftGM,bm.fa.rightGM)
        bm.samples[:,idx] =  bm.samples[:,idx] + 1 
        # keep only where GM+braimask
        bm.samples[:,bm.samples[0,:] != 2] = 0 
        bm.samples[:,bm.samples[0,:] == 2] = 1
        return bm
    
    brainmask = get_BMwithGM(dstmp)
    for filename in DSlist:
        if type(DSlist[0]) == str: bm_in = h5load(filename)
        else: bm_in = filename
        
        bm = get_BMwithGM(bm_in)
        #print "DEBUG: bm.shape, brainmask.shape:", bm.shape, brainmask.shape
        brainmask.S[:,:] = brainmask.S * bm.S
    
    #reduce dimension to use it as mask (shape (1,XX) can't be used)
    brainmask.samples = brainmask.samples.squeeze() 
    dsback2nifti = map2nifti(brainmask)
    
    # check if file exist and delete in this case because nb.save can't overwrite
    if os.path.isfile(dsback2niftiPath): os.remove(dsback2niftiPath) 
    
    nb.save(dsback2nifti, dsback2niftiPath)
    print "new maskfile: ",  dsback2niftiPath
    return os.path.abspath(dsback2niftiPath), brainmask 

def _passAttr(refDS, outDS, attrlist):
    """
    just pass all the given attributes to the other ds, IF they exist
    .. there should be a function doining this already .. hmmmmm
    """
    
    for attr in attrlist:
        if hasattr(refDS.a,attr):
            outDS.a[attr] = refDS.a[attr]
    return outDS

def _getSumROIMask(DSlist, intersection=False):
    """
    1.) sums up ROImask from subjects 
    2.) also give ds0 and ds1 back masked
    """
    
    from mvpa2.suite import map2nifti, fmri_dataset
    import nibabel as nb
    import os
    import numpy as np
    from HA_extFunctions import _passAttr
    
    def maskAndMoveAttr(ds,mask,tmp): 
        ds = ds.copy()
        ## also give ds0 and ds1 back masked
        # convert to nii dataset
        if tmp.shape[1] == ds.shape[1]:
            ds2nii = map2nifti(tmp,ds)
        else:
            ds2nii = map2nifti(ds)
        
        ds_masked = fmri_dataset(samples=ds2nii,mask=mask)
        # hand over attribures. Can't take ds_masked.a = ds.a because
        # mapper aren't the same anymore thus generate errors when using new mapper later
        # TODO: make useable and beauty
        attrlist = [
                    'deg2cos_channels', 'subj','degORcosORts','rightROI',
                    'leftROI', 'condition', 'usePh','brainmask', 'leftGM',
                     'rightGM'
                    ]
        ds_masked = _passAttr(ds, ds_masked, attrlist)
        return ds_masked

    def getMask(maskds):
        tmp0 = fmri_dataset(samples=maskds.a.rightROI)
        tmp1 = fmri_dataset(samples=maskds.a.leftROI)
        mask = tmp0.copy(); mask.samples = tmp0.S + tmp1.S
        mask.S[mask.S > 0] = 1
        return mask
    
    mask2niftiPath = "sumMaskROI.nii.gz"
    
    # get all the mask from right and left ROI
    print "using left and rightRoi for masking"
    Masklist = [getMask(ds) for ds in DSlist]

    ## sum them up to see overlap
    summask = Masklist[0].copy()
    # ..get overlap
    summask.samples = summask.samples * 0  # as basis for adding up masks
    for mask in Masklist:
        summask.samples = np.round(summask.S + mask.S)
    #reduce dimension to use it as mask (shape (1,XX) can't be used)
    summask.samples = summask.samples.squeeze() 

    # if intersection mask should be used, just keep voxel which are in both
    # datasets ds0 and ds1, thus have value=2
    if intersection==True: summask.S[summask.S < len(DSlist)] = 0

    mask2nii = map2nifti(summask)
    
    # check if file exist and delete in this case because nb.save can't overwrite
    if os.path.isfile(mask2niftiPath): os.remove(mask2niftiPath) 
    
    nb.save(mask2nii,mask2niftiPath)
    tmp = fmri_dataset(mask2nii, mask=mask2niftiPath)
 
    print "..use generated mask on them"
    outList = [maskAndMoveAttr(ds, mask2niftiPath,tmp) for ds in DSlist]

    return mask2nii, outList, mask2niftiPath

def _pumpUpDS(inputDSlist,getHeaderDSlist,idx_biggestDS):
    """
    eg: _pumpUpDS(getMapperDS,ProjDS)
    => Here give the trainingDS the shape of the biggest getMapperDS
    by giving it a new mask, if it isn't the biggest getMapperDS anyway.
    The "pumped up" DS contains voxel with zero variance but all that matters 
    is its shape ;)
    """
    
    import numpy as np
    from mvpa2.suite import map2nifti, fmri_dataset
    import nibabel as nb
    from HA_extFunctions import _dimreducer, _getSumROIMask
    import os

    trainingDS_toMaxSub =  inputDSlist[:idx_biggestDS]+inputDSlist[(idx_biggestDS+1):] 
    
    # need to make beautiful!
    idx_2ndBiggestDS = \
        [ds.shape == trainingDS_toMaxSub[np.argmax([d.nfeatures for d in trainingDS_toMaxSub])].shape for ds in inputDSlist].index(True)
    
    tds_nifti = map2nifti(getHeaderDSlist[idx_2ndBiggestDS], inputDSlist[idx_2ndBiggestDS])
    del trainingDS_toMaxSub

    # generate sum mask of biggest (not used) ds and biggest ds
    # -> even if one ds has more voxel, it doesn't mean that all 
    # voxel of the smaller ds are within its mask! 
    # Solution: sum mask of biggest and 2nd biggest      
    ttt0 = getHeaderDSlist[idx_biggestDS].copy()
    ttt1 = getHeaderDSlist[idx_2ndBiggestDS].copy()
    
    # take only 1 Vol (no 4D needed)        
    tmp0, tmp1 = ttt0[0].copy(), ttt1[0].copy()
    tmp0.samples[:,:], tmp1.samples[:,:] = 1, 1   
    
    # get sum mask out of the 2 ds
    _,[tmp0_maskDS,tmp1_maskDS],_ = _getSumROIMask([tmp0,tmp1])
    idx = np.array([(tmp0_maskDS.S == 1) | (tmp1_maskDS.S == 1)]).squeeze()
    maskDS = tmp0_maskDS[:,idx == True]
    maskDS.samples[:,:] = 1
    tmp_niiOut = map2nifti(maskDS) # eg: .shape = (64, 86, 45, 1)
    tmp_niiOutname = "tmpMask_NOTsqueezed.nii.gz"
    
    # check if file exist and delete in this case because nb.save can't overwrite
    if os.path.isfile(tmp_niiOutname): os.remove(tmp_niiOutname) 
    
    nb.save(tmp_niiOut, tmp_niiOutname)
    
    # remove 1 dimension to be able to use it as mask
    tmp_outName = _dimreducer(tmp_niiOutname)

    # re-generate biggest trainingDS with mask of biggest inputDSlist
    tmpDS = fmri_dataset(samples=tds_nifti, mask=tmp_outName)

    inputDSlist[idx_2ndBiggestDS] = tmpDS
    
    return inputDSlist
        
# ==================================================================== #
def listsel(inlist, pos):
    """
    just load a list from h5save and give path containing "pos" (ie. subject number)
    """
    
    from mvpa2.base.hdf5 import h5load
    return [l for l in h5load(inlist) if (l.rfind("_"+pos) != -1)][0] 
   
# ==================================================================== #
def parser(in_DSout,in_ROIvoxel):
    """
    Simply takes Output from retmapFlow OR MovieFlow and parse it to Mapper
    or Proj stuff. Depending on Usage as MapperParser or ProjParser
    (by cloning)
    """
    
    # some placeholder
    return [in_DSout,in_ROIvoxel]

def parserSingle(in_file):
    """
    single file version for placeholder for if("do preproc step")
    .. this is now a mess and just a placeholder
    """
    
    import os
    out_file = os.path.abspath(in_file)
    return in_file 

def parserSmoother(in_file):
    """
    .. also for non standardized output functions ..
    And this is also just a placeholder because leftover from old times
    """
    
    smoothed_file = in_file
    return smoothed_file

# ==================================================================== #
def removeGrandMean(infileSave,doit=True):
    """
    takes the mean tSeries across all volumes for each subject and
    substract it from the tSeries of all voxels
    """
    
    from mvpa2.suite import h5save, h5load
    from numpy import unique
    import os
    
    if doit == True:
        # some nodes give infileSaveSave as list, some as string.
        # Thus adopt to it for loading:
        if type(infileSave) == list: infileSave = infileSave[0]
        
        ds = h5load(infileSave)
        print ".. processing single subject"
        ds.samples = (ds.S.T - ds.S.mean(axis=1)).T
        
        DSout = "exclGrandMean.h5"
        h5save(DSout,ds, compression=9)
        DSoutpath = os.path.abspath(DSout)
    # otherwise skip it .. (quick and dirty)
    else:
        DSoutpath = infileSave
    return DSoutpath

# ==================================================================== #
def RetMapBash(delLead, keepVol, infile_all, cond_all,
                SUBJECTS_DIR, openfmriDir,brainmask,prefix="",nscans=4.
                bashscript='studyforrest-data-retinotopy/code/process_retmap'):
    """
    RetMap via AFNI run in bash script
    
    here just call the bash script HA_RetMap_Flow_severalScans.sh 
    Returns a list of all preprocessed files in the retmapdir.
    
    Takes list of cond_all and infiles for all scans and selects the
    number of cond and infiles depending on nscans
    """ 
    
    #TODO: - convert bash to nipype -> write a RetinoProc wrapper
    import subprocess
    import os    
    from HA_extFunctions import _findfile
    
    # for fslmaths: create a gauss kernel of sigma mm and perform mean filtering
    smooth = 4
    
    ## define inputs and bash script path itself ##
    #
    # Takes list of cond_all and infiles for all scans and selects the
    # number of cond and infiles depending on nscans
    infile_all  = infile_all[:nscans]
    cond_all    = cond_all[:nscans]
    # extract the actual subj from the filepath:
    subj = infile_all[0].split("/")[-3]
    
    # transform list of strings of infile_all into one string #
    infiles  = ' '.join(infile_all)
    cond_all = ' '.join(cond_all)
    
    # add quotation marks to the strings to make it readable as a list from the bash script
    infiles  = ''.join(["'",infiles,"'"])
    cond_all = ''.join(["'",cond_all,"'"])
    
    # EPIdir depending on nscans
    EPInewDir      = os.path.abspath("") # to save the stuff in nipype workdir
    
    # merge script and all its inputs #   
    cmd  = ' '.join([bashscript,
                    subj,
                    infile_all[0],
                    infile_all[2],
                    infile_all[3],
                    infile_all[1],
                    str(True),
                    brainmask[0],
                    str(smooth), 
                    EPInewDir,
                    str(True), str(0.0), str(0.0), str(101.25), str(-281.25),
                    str(False), "''","''", "''",
                    str(True)
                    ])

    # run script and get 0 if finished without errors #
    out0 = subprocess.call([cmd],shell=True)
    
    # grab files written by HA_RetMap_Flow_severalScans.sh as output #
    #outdir = os.path.join(openfmriDir,subj,'BOLD',EPIdir)
    here = os.path.abspath("")
    phaseFilesPath = os.path.join(here,subj,'post_processing')
    preprocScansPath = os.path.join(here,subj,'pre_processing')
    
    preprocScans = _findfile(preprocScansPath,'*filt_masked.nii.gz')
    phaseFiles   = _findfile(phaseFilesPath,'combined_*.nii.gz')
    print '.. Finished with node: RetMapBash'
    
    # if bash-script didn't finsh clean. throw Erroer
    if out0 != 0:
        #break
        print ">> ERROR <<"
    else: 
        return preprocScans, os.path.abspath(""), phaseFiles
        
# ==================================================================== #
def collect_RetMapDS(subj,subjlist, nscans, retmappathlist, prefix,
                    tsfilename,ROImaskL,ROImaskR,brainmask,leftGM,rightGM,
                    usePhInsteadTseries=False,channels=10):
    """
    collect data drom all subjects and unites them in a dataset dsAll
    
    phases, coh, tseries, brainmask, and roimask are saved
    
    usePhInsteadTseries : to project and compare phaseMaps instead of tSeries
    """
    
    import os
    from mvpa2.suite import fmri_dataset, vstack, h5save 
    from mvpa2.suite import h5save
    from HA_extFunctions import _dimreducer, _findfile, _deg2cosine
    import numpy as np

    print retmappathlist
    
    DSout = []
    
    ## load data from all subjects
    #
    # load phasemaps separately  
    retmappathlist.sort()
    phaseECC = retmappathlist[1]
    phasePOL = retmappathlist[0]
    phaseECCname = os.path.basename(retmappathlist[1])
    phasePOLname = os.path.basename(retmappathlist[0])

    #squeeze to remove additional dimension
    ### TODO: dont do dimreduction every time
    print('removing 5th dimension of afni images')
    print phaseECC
    print phasePOL
    
    # keep 1st for coherence (?) of phases
    cohECC = _dimreducer(phaseECC, 1)
    cohPOL = _dimreducer(phasePOL, 1)
    # keep 0th for phases
    phaseECC = _dimreducer(phaseECC, 0)
    phasePOL = _dimreducer(phasePOL, 0)
    
    if usePhInsteadTseries: nscans = 2
    
    #take only pairs of scans (1&2, 3&4)
    sub = []
    for scan in range(0,nscans,2):
        # find Scans with nametemplate #
        # (avoid to define name for each scan and condition(clw, ccw..)
        tsfilename.sort()
        retmaptsnameECC = os.path.basename(tsfilename[scan])
        retmaptsnamePOL = os.path.basename(tsfilename[scan+1])
        RetMapTs_ecc = tsfilename[scan]
        RetMapTs_pol = tsfilename[scan+1]

        # save phase values to project
        if usePhInsteadTseries:
            cond = "Phase" 
            samplesE = phaseECC
            samplesP = phasePOL
            targetsE = phaseECCname
            targetsP = phasePOLname
        # save tSeries to project
        else:
            cond = "tSeries" 
            samplesE = RetMapTs_ecc
            samplesP = RetMapTs_pol
            targetsE = retmaptsnameECC
            targetsP = retmaptsnamePOL

        # load ECC dataset                                                                          
        print('loading ECC dataset:', samplesE)
        ecc_tmp = fmri_dataset(
                                samples=samplesE,
                                mask = brainmask,
                                add_fa = {
                                    'rightGM': rightGM[0],
                                    'leftGM': leftGM[0],
                                    'rightROI': ROImaskR[0],
                                    'leftROI': ROImaskL[0],
                                    'phaseECC': phaseECC,
                                    'cohECC': cohECC,
                                    'brainmask': brainmask
                                    },
                                )
        # transform phase to deg to get better result of hypal and also 
        # get easier statistics without circualr statistics
        if usePhInsteadTseries:
            ecc_tmp = _deg2cosine(ecc_tmp, channels=channels)
        sub.append(ecc_tmp)
        
        # load POL dataset
        print('loading POL dataset:', samplesP)
        pol_tmp = fmri_dataset(
                                samples=samplesP,
                                mask = brainmask,
                                targets = targetsP,
                                add_fa = {
                                    'rightGM': rightGM[0],
                                    'leftGM': leftGM[0],
                                    'rightROI': ROImaskR[0],
                                    'leftROI': ROImaskL[0],
                                    'phasePOL': phasePOL,
                                    'cohPOL': cohPOL,
                                    'brainmask': brainmask
                                    },
                                )
        # transform phase to deg to get better result of hypal and also 
        # get easier statistics without circualr statistics
        if usePhInsteadTseries:
            pol_tmp = _deg2cosine(pol_tmp, channels=channels)
        sub.append(pol_tmp)
    
    #stack all subjects and scans into one dataset
    subds = vstack((sub))

    # pass mapper from input to dataset of subj
    subds.a.update(sub[0].a)
    subds.a['brainmask'] = brainmask
    subds.a['leftROI']   = ROImaskL[0]
    subds.a['rightROI']  = ROImaskR[0]
    subds.a['leftGM']    = leftGM[0]
    subds.a['rightGM']   = rightGM[0]
    subds.a['condition'] = cond
    subds.a['subj']      = subj
    subds.a['usePh']     = usePhInsteadTseries
    if usePhInsteadTseries:
        subds.a['deg2cos_channels'] = channels
        subds.a['degORcosORts'] = "cos"   # to check if conversion to deg is needed for analysis
    else:
        subds.a['degORcosORts'] = "ts"   # to check if conversion (to deg) is needed for analysis
    
    ## save dataset for each single subject
    DSout = 'projDS.h5'
    h5save(DSout,subds, compression=9)
    return os.path.abspath(DSout)

# ==================================================================== #
def collect_MovieDS(subj, nscans, inFilesList,
        ROImaskL, ROImaskR, brainmask,leftGM, rightGM, cutAtEdge=4, roimask="", rl=False ):
    
    """
    collect data from all subjects and unites them in a dataset dsMovie
    tseries, and brainmask are saved
    
    cutAtEdge: volumes to cut before and after each neignbouring movie segments
    nVol: number of volumes to keep for each scan/Movie segment
    """
    
    # TODO:
    # - don't load roi here -> avoid preprocessing of if just ROI changed
    import os
    from mvpa2.suite import fmri_dataset, vstack, h5save 
    from mvpa2.suite import h5save
    from HA_extFunctions import _findfile, _roitmpl2roiname
    import numpy as np

    # to delete volume an overlapping edges of scans (of the movie)
    def _volCutter(data, scan, cutAtEdge, nscans):
        # cuts the leading and tailing volumes (n=cutAtEdge), depending
        # if it is a middle or corner (last or first) scan
        nVol = data.shape[0]
        
        ## check for corner scans:
        # 1st scan cut only tailing volumes..
        if scan == 0:   
            data = data[:(nVol-cutAtEdge)]
        # last scan cut only leaeding volumes..
        elif scan == (nscans-1): 
            data = data[cutAtEdge:]
        # others cut leading AND tailing volumes due to more overlap to
        # neighbor volumes
        else:
            toKeep = range(cutAtEdge,(nVol-cutAtEdge))
            data = data[toKeep]
        return data

    # sort folderlist to keep order of subj for reading
    inFilesList.sort() 

    ## load data for a single subject
    sub, mapperDS = [],[]
    # older version hand over lists
    if type(ROImaskR) == list: ROImaskR = ROImaskR[0]
    if type(ROImaskL) == list: ROImaskL = ROImaskL[0]
    if type(leftGM) == list: leftGM = leftGM[0]
    if type(rightGM) == list: rightGM = rightGM[0]

    for scan,inFile in zip(range(len(inFilesList)),inFilesList):
        ## use brainmask from ProjDS
        # load scanspecific and subjspecific brainmask 
        print('loading dataset:', inFile)
        if not roimask == "":
            data_tmp = fmri_dataset(
                                    samples=inFile,
                                    mask = brainmask,
                                    chunks = 'subj' + str(subj),
                                    targets = inFile,
                                    add_fa = {
                                        'leftGM': leftGM,
                                        'rightGM': rightGM,
                                        'leftROI': ROImaskL,
                                        'rightROI': ROImaskR,
                                        'brainmask': brainmask,
                                        'roimask': roimask},
                                    )
        else:
            data_tmp = fmri_dataset(
                                    samples=inFile,
                                    mask = brainmask,
                                    chunks = 'subj' + str(subj),
                                    targets = inFile,
                                    add_fa = {
                                        'leftGM': leftGM,
                                        'rightGM': rightGM,
                                        'leftROI': ROImaskL,
                                        'rightROI': ROImaskR,
                                        'brainmask': brainmask},
                                    )
        # cut the volumes at the edges of movieSegments (scanwise!)
        data_tmp = _volCutter(data_tmp, scan, cutAtEdge, nscans)

        # if problem with missing volumes for sub004 (set rl=True)
        if rl: data_tmp = data_tmp[:3468]
        
        # restrict to GM
        data_tmp = data_tmp[:,np.logical_or(data_tmp.fa.leftGM,data_tmp.fa.rightGM)]

        print "DEBUG .. data_tmp.shape: ", data_tmp.shape

        sub.append(data_tmp)
        data_tmp = []
    #stack all subjects and scans into one dataset
    mapperDS = vstack((sub))
    # pass mapper from input to dataset of subj
    mapperDS.a.update(sub[0].a)
    mapperDS.a['brainmask'] = brainmask
    mapperDS.a['leftROI'] = ROImaskL
    mapperDS.a['rightROI'] = ROImaskR
    mapperDS.a['leftGM'] = leftGM
    mapperDS.a['rightGM'] = rightGM
   
    # if roimask is provided cause used
    if not roimask == "": 
        print "DEBUG .. attaching roimask"
        roimaskDS = fmri_dataset(samples=roimask)
        mapperDS.a['roimask'] = roimask
    
    DSout  = 'mapperDS.h5'
    h5save(DSout,mapperDS, compression=9)
    return os.path.abspath(DSout)

# ==================================================================== #
def DSjoiner(subjMapperDS):
    
    """
    takes the datasets of all subjects and merge them into a List of all files
    .. also could use JoinNode, but this is the easy solution FOR NOW.
    joining multiple sources also needs to be tested!
    """
    
    from mvpa2.suite import h5save
    import os
    subjMapperDS.sort()
    DSout  = 'DSList.h5'
    h5save(DSout,subjMapperDS, compression=9)
    return os.path.abspath(DSout)

# ==================================================================== #
def restrictDs2ROIvoxel(dsAll, noROI=False, useCohMask=False, cohThres=0.00, makeROIdsNifti=False):
    """
    mask dsAll with ROIvoxel-mask and save them as ROIvoxel.
    
    short: just mask dsAll with a fa 
    ROIvoxel contains the dsAll restricted to only its ROIvoxel.
    - check if ROIvoxel already generated and skip generating in that case
    - Can handle singleSubject DS or a list of DS
    - !! cohMask only generated for single Subject input!
    """
    
    import numpy as np
    import os
    from mvpa2.suite import h5save,h5load, fmri_dataset, map2nifti
    import nibabel as nb
    from numpy import unique
    from HA_extFunctions import _back2nifti, _dimreducer
    
    if type(dsAll) == list: DS = h5load(dsAll[0])
    else: DS = h5load(dsAll)

    # first restrict to GM, cause same for both further maskings
    if type(DS) == list:
        DS = [ds[:,np.logical_or(ds.fa.leftGM,ds.fa.rightGM)] for ds in DS]
    else: DS = DS[:,np.logical_or(DS.fa.leftGM,DS.fa.rightGM)]

    def get_roivoxel(DS):
        print "DEBUG: .. restirct to GM and ROI"
        ROIvoxel = DS[:,np.logical_or(
            np.logical_and(DS.fa.rightGM,DS.fa.rightROI > 0),
            np.logical_and(DS.fa.leftGM,DS.fa.leftROI > 0),
                )]
        return ROIvoxel

    # Can handle singleSubject DS or a list of DS (if type...):
    ## Mask to voxel where cohThres for POL and ECC is >= cohThres    
    def get_cohVoxel(DS,cohThres):
        print "DEBUG: .. generating cohmask from cohECC and cohPOL"
        CohVoxel = DS[:,np.logical_and(DS.fa.cohPOL >= cohThres, DS.fa.cohECC >= cohThres)]
        
        ## get a cohThresMask for later usage (i.e. masking Movie or MapperDS).
        # only works on single subjects        
        dsback2niftiPath = 'cohMask_thres'+str(cohThres)+'.nii.gz'
        
        cohmask = CohVoxel.copy()
        cohmask.S[:,:] = 1  # set all values that surpassed cohThres to 1 (binarize)

        dsback2nifti = map2nifti(cohmask[0])
        #reduce dimension to use it as mask (shape (1,XX) can't be used)
        dsback2nifti= dsback2nifti.get_data().squeeze()
        headertmp = map2nifti(CohVoxel[0])
        outnifti = nb.Nifti1Image(dsback2nifti, headertmp.get_affine())
        
        # check if file exist and delete in this case because nb.save can't overwrite
        if os.path.isfile(dsback2niftiPath): os.remove(dsback2niftiPath) 
        
        nb.save(outnifti, dsback2niftiPath)
        
        #print "DEBUG: .. new maskfile: ",  dsback2niftiPath
        return CohVoxel, dsback2niftiPath
    
    print "DEBUG: ..  shape before all masking : ", DS.shape
    # if roimask is provided (as for movie data that got cohmask
    # from retmap data within getMovieDS: just apply the cohmask
    if hasattr(DS.a,"roimask"):
        print "DEBUG: .. using roimask"
        ROIvoxel = DS[:,DS.fa.roimask == 1]
        roimaskpath = ""
        print "DEBUG: ..  shape after all masking : ", ROIvoxel.shape
        roiNiftiPaths = [] # placeholder
    # if retmap or other ProjDS is input, thus need to generate a 
    # roimask, rather applying it .. generate mask depending on
    # input params
    else:
        ## decide if ..
        # .. using coh threshold to mask
        print "DEBUG: useCohMask, noROI: ", useCohMask, noROI
        if useCohMask:
            ROIvoxel, _ = get_cohVoxel(DS,cohThres)
            print "DEBUG: ..  shape after cohmasking : ", ROIvoxel.shape
        else:
            print "DEBUG: .. no cohmasking"
            ROIvoxel = DS
        
        # .. using ROI for masking
        if noROI: # just mask with graymask as done in first step, not more
            ROIvoxel = ROIvoxel
            print "DEBUG: ..no ROI masking"
        else: 
            ROIvoxel = get_roivoxel(ROIvoxel)
            print "DEBUG: ..  shape after ROImasking : ", ROIvoxel.shape
            
        # generate ROImask to just use for mapperDS (inkl. ROI and cohMask)
        print "DEBUG: .. generating ROImask from voxel left over after ROI- and coh-Masking"
        roiVoxelmask = ROIvoxel[0,ROIvoxel.S[0,:] != 0].copy()
        roiVoxelmask.S[:,:] = 1 # binarize to mask

        #reduce dimension to use it as mask (shape (1,XX) can't be used)
        roimaskpath = 'roiVoxelMask_coh+ROI.nii.gz'
        dsback2nifti = map2nifti(roiVoxelmask[0]).to_filename(roimaskpath)
        # reduce dim: i.e. (1, 64, 86, 45) -> (64, 86, 45). Otherwise it can't be loaded as fa
        roimaskpath = _dimreducer(roimaskpath) 
        
        # export masked DS as nifti
        if makeROIdsNifti: 
            print "DEBUG: .. saving masked DS as nifti"
            # just use the original phase values as saved in the fa
            tmp = ROIvoxel.copy()
            tmp.samples = np.array([ROIvoxel.fa.phaseECC, ROIvoxel.fa.phasePOL])
            roiNiftiPaths = _back2nifti(ROIvoxel, 1, tmp, 2, -1, -1, 'Sub'+DS.a.subj, nameadd="")
        else: roiNiftiPaths = []
        
    DSout  = 'ROIvoxel.h5'
    h5save(DSout,ROIvoxel, compression=9)
    return os.path.abspath(DSout), os.path.abspath(roimaskpath), roiNiftiPaths

# ==================================================================== #
def get_refMapperDS(dsAll, nsubj):
    """
    Generating reference dataset with voxel outside (V1)roi 
    """

    from mvpa2.suite import h5save, h5load
    import os
    from HA_extFunctions import _roitmpl2roiname
    
    dsAll = [h5load(d) for d in h5load(dsAll)]
    
    # reference space for spheres is outside the ROImasks (where it is 0)
    # thus mask = brainmask - V1
    refdsAll = [ds[:,(ds.fa.rightROI == 0) & (ds.fa.leftROI == 0)] for ds in dsAll]

    DSout  = 'refdsAll.h5'
    h5save(DSout,refdsAll, compression=9)
    return os.path.abspath(DSout)

# ==================================================================== #
def getMapperANDCommonspace(getMapperDSsave, ProjDS, toSubj, usePh=False, preConn=False):
    """
    generate Commonspace and according Mapper. If the subject in which
    space it will be mapped, shouldn't be left out: just get Mapper and 
    CS once for all subjects (saves a lot a space and time..)
    """
    
    from mvpa2.algorithms.hyperalignment import Hyperalignment
    from mvpa2.suite import h5save, h5load, zscore    
    from mvpa2.mappers.procrustean import ProcrusteanMapper
    from HA_extFunctions import _pumpUpDS 
    import os
    import numpy as np 
    import nibabel as nb
    import copy

    print "DEBUG:"
    print "getMapperDSsave = ",getMapperDSsave
    print "ProjDS = ",ProjDS
    print "toSubj = ",toSubj

    getMapperDSlist = h5load(getMapperDSsave) # load filelist
    ProjDSlist      = h5load(ProjDS)
    nsubj           = len(ProjDS)

    # load all from thier list
    if type(getMapperDSlist[0]) == str:
        getMapperDS = [h5load(filename) for filename in getMapperDSlist]
    else:
        getMapperDS = copy.copy(getMapperDSlist)

    if type(ProjDSlist[0]) == str:
        ProjDS = [h5load(filename) for filename in ProjDSlist]
    else:
        ProjDS = [filename for filename in ProjDSlist]

    """ hyperalign 
    

    biggest trainingDS mustn't be smaller then biggest getMapperDS
    
    If biggest trainingDS is bigger the biggest getMapperDS, this DS
    can't be projected into the CS, because the mapper-size depends
    on the biggest trainingDS.shape.
    
    => Here give the trainingDS the shape of the biggest getMapperDS
    by giving it a new mask, if it isn't the biggest getMapperDS anyway.
    The "pumped up" DS contains voxel with zero variance but al that matters 
    is its shape ;)
    """
    
    # take ProjDS to according getMapperDS to use its mapper to get the 
    # mask-nifti of the biggest getMapperDS
    idx_biggestDS = np.argmax([d.nfeatures for d in getMapperDS])

    # only if biggestDS is excluded from training, "pump up" 2nd biggest ds
    if idx_biggestDS == toSubj:
        # take header for nifti convert from its header
        if preConn: 
            # the ProjDS isn't projected, thus its attrivutes aren't needed
            refDS = getMapperDS  
        else:
            refDS = ProjDS 
        getMapperDS = _pumpUpDS(getMapperDS,refDS,idx_biggestDS)
        # also pump up ProjDS. Otherwise conflict if dsfrom.shape 
        # doesn't fit to m.proj.shape
        ProjDS      = _pumpUpDS(ProjDS,ProjDS,idx_biggestDS)
    ##

    print 'doing hyperalignment for leave-one-out subj number',toSubj
    # get Mapper for later forward and reverse projection
    # can be subset of DS to generate CS from. This builds up the CS.
    # can be used to generate a CS out of SubA&B but to get a mapper 
    # of subC into that CS, without changing th CS by subC
    trainingDS = copy.copy(getMapperDS)
    if (toSubj == -1) or (toSubj == [-1]):
        # take all subj data to get Mapper and Commonspace
        leaveSubjStr = "NO"
    else:
        # exclude "toSubj" from CS-generation
        del trainingDS[toSubj]
        leaveSubjStr = str(toSubj+1)

    ha = []
    # (as proven good in hypaltest.py)
    ha = Hyperalignment(alpha=1,ref_ds=np.argmax([d.nfeatures for d in trainingDS]),
            level2_niter=5,
            alignment=ProcrusteanMapper(demean=False,space='commonspace',svd='dgesvd', oblique_rcond=-1.0)) 
    ha.train(trainingDS)
    print "..training done!"

    print "..getting mapper"
    mapperAll = ha(getMapperDS)
    
    #print "DEBUG: .."
    #[a.shape for a in ProjDS]
    #[a.shape for a in getMapperDS]
    #[a.proj.shape for a in mapperAll]
    #[a.shape for a in trainingDS]

    # bring all subjects into Commonspace
    # They have to have same number of subjects for ProjDS and mapperAll.  
    CS = []
    for m,dsfrom in zip(mapperAll,ProjDS):
        print "..mapping into CS"
        cs = m.forward(dsfrom)
        if usePh == False:
            print "..doing zscoring of tSeries CS"
            _ = zscore(cs)
        CS.append(cs)
    
    savekey = '_toSu'+str(toSubj)
    MapperOut      = 'MapperAllminus1'+savekey+'.h5'
    CommonspaceOut = 'Commonspace'+savekey+'.h5'
    CSout          = 'CSsubj'+savekey+'.h5'   
    h5save(MapperOut,mapperAll, compression=9)
    h5save(CommonspaceOut,ha.commonspace, compression=9)
    h5save(CSout,CS, compression=9)
    return os.path.abspath(MapperOut), \
            os.path.abspath(CommonspaceOut), \
            os.path.abspath(CSout)

# ==================================================================== #
def getSphereMeanDS(refDSsave,sparse_radius,nproc=10):
    """
    get spheres within the brainmask of each subject, calculates the 
    mean tSeries for each and gives this as output (..for later hyperhyper .. yeah!)
    
    Inputs ::
    refDSsave     : h5save of list of datasets to calculate sphere means from
    sparse_radius : radius of sphreres (error if to big)
    nproc         : number of cpu to use to calculate seed_mean
    
    Outputs ::
    path to saved h5file with sphere mean tSeries of every subject
      (may vary in number of spheres across subjects due to different 
       size of brains)
    """

    from HA_extFunctions import _get_seed_mean
    from mvpa2.suite import h5load, h5save, fmri_dataset, map2nifti
    import os
    
    refDSlist = h5load(refDSsave)    
    nsubj = len(refDSlist)
    
    Sm_suSpace = []
    for su in range(nsubj):
        print ". processing subject: ",su
        refds = h5load(refDSlist[su])
        # subjectwise spheres this subjectspecific brainmasks
        print "..get brainmask of subject"
        brainmask = refds.a.brainmask[0]

        # prepare masking (adopt shape from 4D to 3D, otherwise is not working)
        bm = refds.copy()
        bm = bm[0]
        bm.samples = bm.S[0,bm.samples[0,:] != 0]

        print "..get spheres"
        seed_mean, _ , _ = _get_seed_mean(bm,sparse_radius, nproc)
        
        print ".. get seed_mean of subject: ", su
        tmp = seed_mean(refds)
        # maybe useful for later steps ..
        tmp.a['sparse_radius'] = sparse_radius
        tmp.a['leftGM'] = refds.a.leftGM
        tmp.a['rightROI'] = refds.a.rightROI
        tmp.a['brainmask'] = refds.a.brainmask
        tmp.a['rightGM'] = refds.a.rightGM
        tmp.a['leftROI'] = refds.a.leftROI
        if hasattr(refds.a,"deg2cos_channels"):
            tmp.a['deg2cos_channels'] = refds.a.deg2cos_channels
        if hasattr(refds.a,"degORcosORts"):
            tmp.a['degORcosORts'] = refds.a.degORcosORts
        if hasattr(refds.a,"subj"):
            tmp.a['subj'] = refds.a.subj
        if hasattr(refds.a,"usePh"):
            tmp.a['usePh'] = refds.a.usePh
        if hasattr(refds.a,"condition"):
            tmp.a['condition'] = refds.a.condition

        Sm_suSpace.append(tmp)

    # save it as h5
    SeedmeanOut = 'SeedmeanSubjectwise.h5'
    h5save(SeedmeanOut,Sm_suSpace,compression=9)
    return os.path.abspath(SeedmeanOut)

# ==================================================================== #
def make_Connectome(sparse_radius, ROIvoxel, refDS, brainmask="needFix", 
                    refDSisMean=False, nproc=10,connectomeMethod="seed_meanAND_SubjTSmean"):
    """
    generate Connectome with spheres outside roi (i.e. wholebrain - ROImask)
    
    checks if connectome is already finished and load in that case
    - method: seed_mean, pca_ica, atlasROI
    
    output :: 
       - Connectome_all: Connectomes orderd by subjectorder 
    """
    
    # TODO:
    # - unify all the different cases by outsourcing the differneces
    # - try to implement subjectwise brainmasks
    # - get seedmask and all other stuff from one fixed volume (as for atlas),
    #   NOT subjectwise (as now??)
    
    import os
    import numpy as np
    import scipy.spatial.distance as sd
    from mvpa2.datasets import Dataset, vstack
    from mvpa2.suite import h5save, h5load, fmri_dataset, zscore
    from HA_extFunctions import _get_seed_mean, _getIntersectionBrainMask
    from sklearn.decomposition import FastICA, PCA

    if type(refDS) == list: refDSlist = h5load(refDS[0])
    else:                   refDSlist = h5load(refDS)
    
    if type(ROIvoxel) == list: ROIvoxelList = h5load(ROIvoxel[0])
    else:                      ROIvoxelList = h5load(ROIvoxel)

    nsubj  = len(ROIvoxelList)    
    print ">>> make_Connectome"   
    print "DEBUG nsubj =",nsubj   
    print "DEBUG refDS = ",refDS       
    print "DEBUG ROIvoxel = ",ROIvoxel
    
    # if no Seed_mean is generated, cause the Sphere means are
    # the input (alias preConenctome) skip the generation of intersectionMask,
    # otherwise (if refDSisMean is False), generate intersectionMask
    if refDSisMean is False: # otherwise spare Seed_mean stuff cause not needed
        ### get intersectionMask ##
        ## a brainmask to use on all subj datasets. It is the essenece of Voxel
        ## which are present in all subj, thus voxel of subjBrains that are 
        ## only present in this or also other subj, will be excluded.
        # => intersection Mask
        print "..getting IntersectionMask"
        intersectionMaskPath, intersectionMask = _getIntersectionBrainMask(refDSlist)
        bm = fmri_dataset(samples=intersectionMaskPath,mask=intersectionMaskPath)
        seed_mean, _ , _ = _get_seed_mean(bm,sparse_radius, nproc)  # a.voxel_dim bigger then voxel indices reach
        # ==================================================================== #
        # get mean in each seed, subjectwise
        tt, sm_su = [], []
        for su in range(nsubj):
            if type(refDSlist[su]) == str: refds = h5load(refDSlist[su])
            else: refds = refDSlist[su]
            
            print ".. getting seed_mean of subject: ", su
            print "DEBUG: refds.a.leftGM: ", refds.a.leftGM
            print "DEBUG: refds.shape:", refds.shape
            tmp = seed_mean(refds)            
            sm_su.append(tmp.S)
            tt.append(tmp)

        # get a mean tSeries across all subj
        if connectomeMethod == "seed_meanAND_SubjTSmean":
            sm_su = np.array(sm_su)
            # get mean tSeries across all subj
            m = np.mean(sm_su,axis=0)
            print m.shape
            meanDS = tt[0]
            meanDS.S[:,:] = m
            ss = meanDS
            ss.fa['voxel_indices_spheres'] = ss.fa.voxel_indices
            dsmIn = ss
        if connectomeMethod == "sphere_mainPCAcomp":
            print "..reducing number of elements by PCA and ICA"
            # reduce data by 1st doing PCA and then ICA on these
            # components
            pca = []
            pca = PCA(n_components=len(tmp)) # define
            pcaComp = pca.fit_transform(tmp) # reduce dimensions
            scands = h5load(ROIvoxelList[su])

            ica=[]
            ica = FastICA(n_components=len(tmp),random_state=0)
            icaComp = ica.fit_transform(pcaComp)
            ss = icaComp
    else:
        intersectionMaskPath = "none"
                
    Connectome, DSM = [],[]
    for su in range(nsubj):
        print 'debug: su, nsubj', su, nsubj
        # if no seed_means are needed ot the refDS is already any kind,
        # of Connectome: Load the refDS directly as dsmIn to correlate
        # the ROIvoxel tSeries with to get a Connectome in the end (ConenctomeOut)
        if refDSisMean is True:
            if type(refDSlist[0]) == str: dsmIn = h5load(refDSlist[su])
            else: dsmIn = refDSlist[su]
        
        # ds inside ROI. Correlate MapperDS with 
        scands = h5load(ROIvoxelList[su])

        if (connectomeMethod == "seed_mean_Only") and (refDSisMean is True):
            # get the Connectome subjectwise, 
            # thus dsmIn is NOT meanTS across all subj
            # Instead, average the Connectome across subjects
            # 1.) load seed mean to get subjectspecific Connectome
            dsmIn = tt[su]
        
        print "shape of ds to compare: ", scands.shape, dsmIn.shape   
        dsm = sd.cdist(dsmIn.S.T,scands.S.T, 'euclidean')
        
        # save it for eventual later averaging across all subj
        DSM.append(dsm) 
        
        # >> DEBUG
        print ">>>> mean(dsm): ", np.mean(np.diag(dsm))
        print ">>>> max(scands): ", np.max(scands.S)
        print ">>>> var(dsm): ", dsm.var()

        print "connectome.shape: ", dsm.shape
        connTmp = Dataset(dsm, sa={
                            'chunks': np.ones(dsm.shape[0]),
                            },
                            fa={
                            'voxel_indices': scands.fa.voxel_indices,
                            })
        con = connTmp
        # keep brainmask for later Connectome
        con.a['brainmask'] = scands.a.brainmask
        Connectome.append(con)  

    # save it as h5
    ConnectomeName  = 'Connectome11.h5'
    h5save(ConnectomeName,Connectome,compression=9)
    ConnectomeOut = os.path.abspath(ConnectomeName)
    
    print "DEBUG: "
    print "ConnectomeOut = ", ConnectomeOut
    print " ==== end make_Connectome ==== "
    return ConnectomeOut, intersectionMaskPath
            
# ==================================================================== #
def leave1out_proj(Mapper,CS,ProjDS, toSubj, nVol, makenifti,
                    nscansMapper,makeMeanCS=False):
    """
    leave one (of the subjects) out and project into left out subject space
    
    build commonspace 1st time without 1 subj, then extend it with that 
    subj and project back into its original space
    - use 'ProjDS' data to project this subjects' data 
      into the leave-subj space
    - The getMapperDS can be anything (RetMap, Connectome etc.).
      It is just used to generate the mapper for forward and reverse projection 
      of the ProjDS-Datasets. Its Commonspace is saved for further
      (yet unknown) purposes.
    - plot color.. subA vs SubB voxel cdist
    """

    #TODO: 
    # - maybe take mean of all keep-subj in CS to project to leave-subj space?
    # - only one comparison for meanCS (all fromSubj-CS are the same due to meanCS

    from HA_extFunctions import (
                    _back2nifti, _shifted_cosine, _ds2deg, _getCosCdist,
                    _getMeanCS
                    )
    import numpy as np
    import scipy.spatial.distance as sd
    from mvpa2.suite import h5save, h5load, vstack
    import os
        
    CS = h5load(CS[toSubj])
    Mapper = h5load(Mapper[toSubj])
    ProjDSlist = h5load(ProjDS)
    
    # load all from thier list
    ProjDS = [h5load(filename) for filename in ProjDSlist]
    
    # use header of the subject(space) it will be projected in
    toCs = CS[toSubj]
    
    ## get number of Volumes per scan ##
    nsubj           = len(CS)
    subjList        = range(nsubj)

    # get mean of CS excluding the toSubj
    if makeMeanCS: CSexlToSubj = CS[:toSubj] + CS[(toSubj+1):]

    # project every Subj into "toSubj" andCorrelate data inside CS
    Backproj2toSubj = [['backproj','fromSubj','toSubj']]
    CScorrToSubj = [['CScorrToSubj','fromSubj','toSubj']]
    CompToSubj = [['corrToSubj2fromSubj','fromSubj','toSubj']]
    outnifti = [['savedniftis','toSubj','fromSubj','nscansMapper']]    
    CompToSubj_data = [['dsback','dsorg','toSu']]
    # save just filenames without further information. More complex lists
    # are more complicated to give into further processing
    outnifti_raw, nameAdd = [], ""    
    for fromSubj in subjList:
        cs = CS[fromSubj]
        print "..processing subj ",fromSubj, " to Subj ",toSubj

        ## voxelwise cdist in CS for all subj -> compare to cdist in ProjDS-Dataset
        print "DEBUG: tosubj.shape, fromshub.shape: ", toCs.shape, cs.shape 
        
        # reverse projection from CS back into subjectspace
        #-->> THE FOLLOWING IS WRONG!!!: mapper = Mapper[fromSubj], because I want to project TOSUBJ!!!
        mapper = Mapper[toSubj] # this is right
        
        # if using meanCS, overwrite subjectspecific cs
        if makeMeanCS:
            cs = _getMeanCS(CSexlToSubj, cs)
            nameAdd = "CSmean"
        
        print "DEBUG: .. mapperAll[toSubj].proj.shape: ",mapper.proj.shape
        backproj = mapper.reverse(cs)
                   
        ## update attributes with data from toCS dataset
        backproj.a.update(ProjDS[toSubj].a)
        print "DEBUG: .. backproj.fa.keys", backproj.fa.keys()
        
        Backproj2toSubj.append([backproj,fromSubj,toSubj])
        print ">> makenifti = ",makenifti
        # export dataset, split into it scans pieces and convert2nifti #
        if makenifti == True:
            
            ## if phasevalues (thus converted): transform phase cosine values back to deg        
            if (hasattr(backproj.a,"usePh") == True):
                if  backproj.a.degORcosORts == "cos":
                    # if cos values, transform to deg
                    # if phasevalues projected, nscansProj is probably smaller
                    nscansProj = backproj.shape[0]/backproj.a.deg2cos_channels
                    print "DEBUG: .. convert cos2deg"
                    niftiDS = _ds2deg(backproj)
                # if usePh = True, nChannelsnVol are used, but just 2 phasevals,
                # thus set nVol = 1 for back2nifti, AFTER cos2deg!
                nVol = 1
            else: nscansProj = ProjDS[0].shape[0] / nVol
            
            print '..splitting dsback into its scanparts and convert them into nifti'
            print "DEBUG: backproj.shape: ",niftiDS.shape
            print "DEBUG: fromSubj, toSubj: ", fromSubj, toSubj
            savedniftis = _back2nifti(toCs, nVol, niftiDS, nscansProj,
                                        nscansMapper,fromSubj, toSubj,nameAdd)
            outnifti.append(zip(
                savedniftis,
                len(savedniftis)*[toSubj],
                len(savedniftis)*[fromSubj],
                len(savedniftis)*[nscansProj]
                ))
            outnifti_raw.append(savedniftis)
        # compare the data of subject "fromSubj" projected to "toSubj"-space
        # with the original "toSubj" data
        ##
        # for working comparison, ignore Voxel outside brain,
        # thus voxel with zero variance. They also generate NaN in cdist('co')                        
        projDS = ProjDS[toSubj]
        
        ## get intersection mask
        #print "..getting IntersectionMask"
        #intersectionMaskPath, intersectionMask = _getIntersectionBrainMask([backproj,projDS])
        
        #from mvpa2.suite import fmri_dataset
        #bm = fmri_dataset(samples=intersectionMaskPath,mask=intersectionMaskPath)
        #backproj_masked = fmri_dataset(samples=backproj,mask=intersectionMaskPath)
        #refdsAll = backproj[:,(bm.samples == 1)]
        
        #_pumpUpDS([backproj,projDS],[projDS],0)
        
        idx = np.array([(backproj.S.T.var(1) == 0) | (projDS.S.T.var(1) == 0)]).squeeze()
        dsback = backproj[:,idx == False]
        dsorg = projDS[:,idx == False]
        if hasattr(backproj.a,"usePh") == True:
            print "DEBUG: .. dsback.shape, dsorg.shape", dsback.shape, dsorg.shape
            compToSubj = _getCosCdist(dsback,dsorg)
            csCorrToSubj = _getCosCdist(toCs,cs)
        else:
            print '>>> Just correlating across whole tSeries; NOT scanwise!!'
            compToSubj = np.mean(np.diag(sd.cdist(dsback,dsorg,'co')))
            csCorrToSubj = np.mean(np.diag(sd.cdist(toCs.S.T,cs.S.T,'co')))
        ##
        
        print "compToSubj: ",compToSubj
        print "csCorrToSubj: ",csCorrToSubj
        CompToSubj.append([compToSubj,toSubj])
        CScorrToSubj.append([csCorrToSubj,fromSubj,toSubj])
    
    # all the same for all fromSu, because meanOthers used ->
    # other are all subj except toSubj
    CompToSubj_data.append([[dsback,dsorg,toSubj]])
        
    savekey = \
        '_'+str(nscansMapper)+'scans'+ \
        '_leavesubjnr'+str(toSubj+1)+'of'+str(nsubj) +'subj.h5'
    print ".. saving with key: ", savekey
    Backproj2toSubjOut = os.path.join('Backproj2toSubj'+savekey)
    CScorrToSubjOut    = os.path.join('CScorrToSubj'+savekey)
    CompToSubjOut      = os.path.join('CorrToSubj2fromSubj'+savekey)
    CompToSubj_dataOut = os.path.join('CorrToSubj2fromSubj_data'+savekey)

    # save outfiles and return it
    Backproj2toSubjSave = Backproj2toSubj
    CScorrToSubjSave    = CScorrToSubj
    CompToSubjSave      = CompToSubj
    CompToSubj_dataSave = CompToSubj_data
    h5save(Backproj2toSubjOut,Backproj2toSubjSave, compression=9)
    h5save(CScorrToSubjOut,CScorrToSubjSave, compression=9)
    h5save(CompToSubjOut,CompToSubjSave, compression=9)
    h5save(CompToSubj_dataOut,CompToSubj_dataSave, compression=9)
    return os.path.abspath(Backproj2toSubjOut), \
            os.path.abspath(CScorrToSubjOut), \
            os.path.abspath(CompToSubjOut), \
            outnifti,outnifti_raw, \
            os.path.abspath(CompToSubj_dataOut)

# ==================================================================== #
def CompareAnatAlign2HyperAlign(ProjDSsave,Backproj2toSubj,usePhInsteadTseries):
    """
    get if there is improvement or not by hyperalignment, compared to anatomical alignment
    """

    from mvpa2.suite import h5load, h5save
    import scipy.spatial.distance as sd
    import numpy as np
    from HA_extFunctions import (
        _roitmpl2roiname, _getSumROIMask, _ds2deg, _get_avg_euc_dist, _getCosCdist,
        )
    import copy
    import os
    
    print "DEBUG: "
    print "ProjDSsave = ", ProjDSsave
    print "Backproj2toSubj = ", Backproj2toSubj
    
    ProjDSlist = h5load(ProjDSsave)
    Backproj2toSubj.sort()
    nsubj = len(ProjDSlist)
    
    # compare the "toSubj", thus the aims with the postHPYAL and the other 
    # preHypal datasets to get the anatomical aligned correlation of the ds
    # and the correlation after HYPAL
    hypalcdistMean = []
    Anat2Hypal = [['anat-hypal[ECC,POL]','anatcdist[ECC,POL]','hypalcdist[ECC,POL]','fromSubj','toSubj', 'nVoxel[hypal,anatal]']]
    print "anat2hypal (positve=improvement by hypal), fromSubj, toSubj: "
    for toSubj in range(nsubj):
        refds = h5load(ProjDSlist[toSubj])
        # load dataset from a list, thus "toSub" defines index of toSubj in list
        tmpPostHypal = h5load(Backproj2toSubj[toSubj])[1:]        
        
        for fromSubj in range(nsubj):
            preHypal = h5load(ProjDSlist[fromSubj]) # subj in its own (anatomical) space
            postHypal = tmpPostHypal[fromSubj][0]          
            
            ### uncomment to use a mask of intersection between pre- and postHypal (very restrictive!)
            ## select voxel which are in all prejections  ##
            _,[refds_masked,postHypal_masked],_ = _getSumROIMask([refds,postHypal])
            idxPP = np.array([
                (refds_masked.S.T.var(1) == 0) | (postHypal_masked.S.T.var(1) == 0)
                ]).squeeze()
            fromSuDsHypal = postHypal_masked[:,idxPP == False].copy()
            toSuDs        = refds_masked[:,idxPP == False].copy()
            
            ## voxelwise cdist for all subj: 
            # compare "pre vs pre" and "pre vs post" to get the qualtiy of
            # alignment of anatomical and hyperaligment
            if usePhInsteadTseries:
                print ">>>>>!!! get ECC&POL with cosine Transformed values for comparablility "
                channels = refds_masked.a.deg2cos_channels
                "DEBUG: ..processing cdist hypalcdist"
                print "DEBUG: ..  toSuDs.shape, fromSuDsHypal.shape :", toSuDs.shape, fromSuDsHypal.shape
                hypalcdist = _getCosCdist(toSuDs, fromSuDsHypal, channels)
                #print "DEBUG: ..hypalcdist, fromSubj, toSubj: ", hypalcdist, fromSubj, toSubj
            else: 
                print "DEBUG: ..processing cdist hypalcdist"
                hypalcdist = np.mean(np.diag(sd.cdist(toSuDs,fromSuDsHypal,'co')))  #shape: vol,nROIvoxel
            print "hypalcdist = ", hypalcdist
            Anat2Hypal.append([[], [], hypalcdist ,fromSubj,toSubj, [preHypal.shape[1],postHypal.shape[1]]])
    
    savekey = '_'+str(nsubj) +'subj.h5'
    print ".. saving"
    DSout = os.path.join('Anat2Hypal'+savekey)
    h5save(DSout,Anat2Hypal, compression=9)
    print "DSout = ", os.path.abspath(DSout)
    return os.path.abspath(DSout)

# ==================================================================== #
def niftiselector(subj,outnifti,subjlist):
    """
    selct the niftis from the outnifti list which should be projected on
    surface. Need this function to extract just the filenames and paths.
    ToSubj and fromSubj information is not nedded.
    """
    
    import numpy as np
    toSubj = subjlist.index(subj)          
    toSubjTmp = outnifti[toSubj]
    toSubjNiftilist = list(np.hstack(toSubjTmp))

    return toSubjNiftilist

# ==================================================================== #
def getAnatCdist(DSlist, makenifti=False):
    """
    Take the phasemaps (in grptmpl-space) and correlate them (in grptmpl-space)
    """

    from mvpa2.suite import h5save, h5load
    from HA_extFunctions import (_getCosCdist, _getMeanCS, _getSumROIMask,
            _ds2deg, _back2nifti)
    import os
    import numpy as np
    
    DS = [h5load(filename) for filename in h5load(DSlist)]

    ## get sum mask from all subj ROI to be able to mean them ##
    # before different shaped ds can't be averaged into one ds
    maskfile,DS_sumMasked,_ = _getSumROIMask(DS)
    
    ## get number of Volumes per scan ##
    nsubj           = len(DS)
    subjList        = range(nsubj)
    
    AnatCdistSummary = [['anatCdist_meanOther','toSubj']]
    DSmeanOther = [['DSmeanOther','toSubj']]
    CompToSubjOut_data = [['dsback_masked','dsorg','toSubj']]
    outnifti = [['savedniftis','toSubj','fromSubj','nscansMapper']]    
    outnifti_raw = []
    for toSubj in subjList:

        DSOther = DS_sumMasked[:toSubj] + DS_sumMasked[(toSubj+1):]
        ds = DS_sumMasked[toSubj] # reference Subject to compare with
        dsMean_other = _getMeanCS(DSOther, ds)
        
        # give nifti of meanDS_other for each toSubj
        # -> toSubj is not included in dsMean_other!
        if makenifti == True:
            ## if phasevalues (thus converted): transform phase cosine values back to deg        
            if (hasattr(dsMean_other.a,"usePh") == True):
                if  dsMean_other.a.degORcosORts == "cos":
                    # if cos values, transform to deg
                    # if phasevalues projected, nscansProj is probably smaller
                    nscansProj = dsMean_other.shape[0]/dsMean_other.a.deg2cos_channels
                    print "DEBUG: .. convert cos2deg"
                    ds2nii = _ds2deg(dsMean_other)
                else: ds2nii = dsMean_other
                
                # if usePh = True, nChannelsnVol are used, but just 2 phasevals,
                # thus set nVol = 1 for back2nifti, AFTER cos2deg!
                nVol = 1
            else: 
                nscansProj = DS[0].shape[0] / nVol
                ds2nii = dsMean_other
            
            print '..splitting dsback into its scanparts and convert them into nifti'
            #print "DEBUG: backproj.shape: ",niftiDS.shape
            #print "DEBUG: fromSubj, toSubj: ", fromSubj, toSubj
            savedniftis = _back2nifti(ds2nii, nVol, ds2nii, nscansProj,
                                        -1,-1, toSubj,"anatMean")
            outnifti.append(zip(
                savedniftis,
                len(savedniftis)*[toSubj],
                len(savedniftis)*[-1],
                len(savedniftis)*[nscansProj]
                ))
            outnifti_raw.append(savedniftis) 
    
        # for working comparison, ignore Voxel outside brain,
        # thus voxel with zero variance. They also generate NaN in cdist('co')                        
        idx = np.array([(ds.S.T.var(1) == 0) | (dsMean_other.S.T.var(1) == 0)]).squeeze()
        ds_masked = ds[:,idx == False]
        dsMean_other_masked = dsMean_other[:,idx == False]
        if hasattr(ds_masked.a,"usePh") == True:
            print "DEBUG: .. dsback.shape, dsorg.shape", ds_masked.shape, dsMean_other_masked.shape
            compToSubj = _getCosCdist(ds_masked,dsMean_other_masked)
        ##
        AnatCdistSummary.append([compToSubj, toSubj, dsMean_other_masked.shape])
        DSmeanOther.append([dsMean_other,toSubj])
        CompToSubjOut_data.append([[dsMean_other_masked, ds_masked,toSubj]])
    savekey = \
        '_nSubj'+str(nsubj)+'.h5'
    print ".. saving with key: ", savekey
    AnatCdistSummarySave = os.path.join('AnatCdistSummary'+savekey)
    DSmeanOtherSave = os.path.join('DSmeanOther'+savekey)
    CompToSubjOut_dataSave = os.path.join('CorrToSubj2fromSubj_data'+savekey)
    
    # save outfiles and return it
    h5save(AnatCdistSummarySave,AnatCdistSummary, compression=9)
    h5save(DSmeanOtherSave, DSmeanOther, compression=9)
    h5save(CompToSubjOut_dataSave, CompToSubjOut_data, compression=9)

    return os.path.abspath(AnatCdistSummarySave), \
           os.path.abspath(DSmeanOtherSave), \
           outnifti,outnifti_raw, \
           os.path.abspath(CompToSubjOut_dataSave)
