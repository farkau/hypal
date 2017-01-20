## FK2014: all teh different workflows used for the anat_cdist or
# hyperalignment meta workflow

#############################################
## See LICENSE for licensing and copyright ##
#############################################

import nipype.pipeline.engine as pe # pypeline engine
import nipype.interfaces.utility as util
from HA_configLoader_v2 import configLoader
from nipype import Node, Workflow, SelectFiles
from nipype.interfaces.utility import Function, Select

########################################################################
# RetMap flow (PreProc, RetMapping, Datacollection to DS)
########################################################################
def create_rmFlow(cfgFile,name='retMapFlow',useGrpTmpl=False):
    """ Do Retinotopic mapping with AFNI (via bash..still), collect the data into Datasets
        and give also Dataset restricetd just to ROIvoxel.
        The loaded cfg tells which processing steps should be made 

    Parameters
    ----------

    name : name of workflow (default: retMapFlow)
    cfgFile : path to used config file (i.e. "/AbsPathToFile/hypal.cfg"

    Inputs::
    inputspec.subj : subj to process

    Outputs::
    outputspec.DSout           : (Phasemap &) tSeries dataset, inkl. ROI and Brainmask
    outputspec.ROIvoxel        : DS dataset restriced to only ROIvoxel
    outputspec.preprocScans    : path to phase files from retmapping (eg. ccw, con, clw, exp)
    outputspec.outdir          : path to parent folder were preprocScans are saved
    outputspec.phaseFiles      : merged phasemaps for eccentricity and polar angle (2 files)
    outputspec.surfacefiles_R  : surface projection of phase files (right hemis.)
    outputspec.surfacefiles_L  : surface projection of phase files (left hemis.)
    outputspec.coreg2anat      : coregistered phaseFiles (nifti)
    outputspec.DSout           : Phasemaps and tSeries of Retmapping saved in Dataset
    outputspec.ROIvoxel        : masked DSout. Masking used roimasks
    outputspec.roimasks        : mask to mask DSout (nifti)
    """

    retMapFlow = pe.Workflow(name=name)

    from HA_extFunctions import (
        RetMapBash, collect_RetMapDS,removeGrandMean, restrictDs2ROIvoxel,
        DSjoiner, parserSingle
        )
    import nipype.interfaces.freesurfer as fs
    from nipype.interfaces import fsl
    def listsel(inlist): return inlist[0]
    def sub2subj(sub): return "sub-"+sub
    
    # ================================================================ #
    ## define workflow specific input
    
    # Infosource
    inputspec = pe.Node(interface=util.IdentityInterface(
        fields=['subj']
        ), name='inputspec')   

    # ---------------------------------------------------------------- #
    ## load cfgFile to get all the settings ##
    # (better then passing them via inputspec(!?))
    #
    print "..loading config from: ", cfgFile
    cfg = configLoader(cfgFile)
    print "..done"
    # ================================================================ #

    # generate list depedning on nscans to use datagrabber depending on nscans
    # without connection to nscans it gets uncoupled from the nscans dependency,
    # thus it generates everything for all subjects
    def getscanlist(nscans):
        scanlist = []
        [(scanlist.append("%03d" % scan)) for scan in range(1,nscans+1)]
        return scanlist
    if useGrpTmpl == False: datagrabber = Node(SelectFiles(cfg.RetMap.templates), name="datagrabber"+name)
    else: datagrabber = Node(SelectFiles(cfg.RetMap_grpTmpl.templates), name="datagrabber"+name)
    datagrabber.inputs.force_lists    = True
    datagrabber.inputs.sort_filelist  = True
    datagrabber.inputs.base_directory = cfg.datadir
    datagrabber.inputs.scanlist       = getscanlist(cfg.RetMap.nscans)
    retMapFlow.connect(inputspec,'subj',datagrabber,'subj')

    ## take data from grptmpl-script and generate Phasemaps and PreProcData ##
    # > SUBJECTWISE < #
    getDSandPhase = pe.Node(name='getDSandPhase'+name, interface=Function(
            function     = RetMapBash,
            input_names  = ['delLead','keepVol','infile_all','cond_all',
                            'SUBJECTS_DIR','openfmriDir','brainmask',
                            'prefix','nscans','retmapPath'],
            output_names = ['preprocScans','outdir','phaseFiles']))
            
    getDSandPhase.inputs.delLead      = cfg.RetMap.delLead
    getDSandPhase.inputs.keepVol      = cfg.RetMap.nVol
    getDSandPhase.inputs.cond_all     = cfg.RetMap.cond_list
    getDSandPhase.inputs.SUBJECTS_DIR = cfg.FSdir
    getDSandPhase.inputs.openfmriDir  = cfg.datadir
    getDSandPhase.inputs.retmapPath   = cfg.retmapPath
    retMapFlow.connect(datagrabber,'func',getDSandPhase,'infile_all')
    retMapFlow.connect(datagrabber,'brainmask',getDSandPhase,'brainmask')

    # ================================================================ #
    
    ## get dataset to project ##
    # stays the same all the time of this pipeline #
    getDS = pe.Node(name = 'getDS'+name, interface=Function(
            function     = collect_RetMapDS,
            input_names  = ['subj','subjlist','nscans','retmappathlist',
                            'prefix', 'tsfilename','ROImaskL',
                            'ROImaskR','brainmask','leftGM','rightGM','usePhInsteadTseries'],
            output_names = ['DSout']),
            unique       = True)
    getDS.inputs.subjlist   = cfg.subjlist
    getDS.inputs.nscans     = cfg.RetMap.nscans
    getDS.inputs.prefix     = cfg.RetMap.prefix
    getDS.inputs.tsfilename = cfg.RetMap.tsfilename
    getDS.inputs.usePhInsteadTseries = cfg.RetMap.usePhInsteadTseries
    retMapFlow.connect(datagrabber,'brainmask',getDS,'brainmask')
    retMapFlow.connect(getDSandPhase,'phaseFiles', getDS,'retmappathlist')
    retMapFlow.connect(getDSandPhase,'preprocScans', getDS,'tsfilename')
    retMapFlow.connect(inputspec,'subj',getDS,'subj')
    retMapFlow.connect(datagrabber,'rightROI',getDS,'ROImaskR')
    retMapFlow.connect(datagrabber,'leftROI',getDS,'ROImaskL')
    retMapFlow.connect(datagrabber,'rightGM',getDS,'rightGM')
    retMapFlow.connect(datagrabber,'leftGM',getDS,'leftGM')

    # parser just has input=output but this way the tow different
    # if-conditions lead to the same output node. This is easier
    # to connect without further changes in the workflow
    parser1 = pe.Node(name = 'parser1'+name,
         interface=Function(function=parserSingle,
         input_names  = ['in_file'], output_names = ['in_file']))
    ## if phasevalues projeced, don't remove mean
    if cfg.RetMap.usePhInsteadTseries:
        retMapFlow.connect(getDS,'DSout',parser1,'in_file')
    ## remove average of all voxels tSeries, else
    else:
        removeMean = pe.Node(name='removeMean'+name,
                            interface=Function(
                                function     = removeGrandMean,
                                input_names  = ['infileSave'],
                                output_names = ['DSout']),
                            )
        retMapFlow.connect(getDS,'DSout',removeMean,'infileSave')
        retMapFlow.connect(removeMean,'DSout',parser1,'in_file')

    ### mask DS with fa from DS and save them as ROIvoxel. ##
    # - skipps generation of ROIvoxel.h5 if it is already generated
    getDsROIvoxel = pe.Node(name='getDsROIvoxel'+name,interface=Function(
            function     = restrictDs2ROIvoxel,
            input_names  = ['dsAll','noROI','useCohMask','cohThres','makeROIdsNifti'],
            output_names = ['DSout','roimask','roiNifti']),
            unique       = True)
    getDsROIvoxel.inputs.cohThres   = cfg.RetMap.cohThres
    getDsROIvoxel.inputs.useCohMask = cfg.RetMap.useCohMask
    getDsROIvoxel.inputs.noROI      = cfg.noROI
    getDsROIvoxel.inputs.makeROIdsNifti = True
    retMapFlow.connect(parser1,'in_file',getDsROIvoxel,'dsAll')
    
    # ================================================================ #
    ### project the retmaps onto surface ###
    #
    coregPh2FS = pe.MapNode(name = 'coregPh2FS'+name,
                            interface=fsl.FLIRT(apply_xfm=True),
                            iterfield=['in_file'])    
    retMapFlow.connect(datagrabber,('anatref',listsel),coregPh2FS,'reference')
    retMapFlow.connect(datagrabber,('regfile',listsel),coregPh2FS,'in_matrix_file')
    
    # for grptmpl, first warp image back from grpspace to subjspace,
    # then applyxfm (in postmat). Als mask with mask from reference space
    if useGrpTmpl:
        warpGrp2subj = pe.MapNode(name = 'warpGrp2subj'+name,
                                interface=fsl.ApplyWarp(),
                                iterfield=['in_file'])
        retMapFlow.connect(getDsROIvoxel,'roiNifti',warpGrp2subj,'in_file')
        retMapFlow.connect(datagrabber,('warpGrp2subj',listsel),warpGrp2subj,'field_file')
        retMapFlow.connect(datagrabber,('refbold3Tp2',listsel),warpGrp2subj,'ref_file')
        
        retMapFlow.connect(warpGrp2subj,'out_file',coregPh2FS,'in_file')
    else:
        # just use applyxfm
        retMapFlow.connect(getDsROIvoxel,'roiNifti',coregPh2FS,'in_file')
    #
    ####

    ## project on surface ##
    # only project, no further coreg
    proj_on_surf_lh = pe.MapNode(name = 'proj_on_surf_lh_'+name,
                        interface=fs.SampleToSurface(
                            hemi='lh',
                            subjects_dir=cfg.FSdir,
                            sampling_method="average",
                            sampling_range=(0.0,0.75,0.001),
                            sampling_units="frac",
                            surface = "white",
                            reg_header=True
                            ),
                            iterfield=['source_file'])
    retMapFlow.connect(coregPh2FS,'out_file',proj_on_surf_lh,'source_file')
    retMapFlow.connect(inputspec,('subj',sub2subj),proj_on_surf_lh,'subject_id')
    proj_on_surf_rh = pe.MapNode(name = 'proj_on_surf_rh_'+name,
                        interface=fs.SampleToSurface(
                            hemi='rh',
                            subjects_dir=cfg.FSdir,
                            sampling_method="average",
                            sampling_range=(0.0,0.75,0.001),
                            sampling_units="frac",
                            surface = "white",
                            #projection_stem="",
                            reg_header=True
                            ),
                            iterfield=['source_file'])
    retMapFlow.connect(coregPh2FS,'out_file',proj_on_surf_rh,'source_file')
    retMapFlow.connect(inputspec,('subj',sub2subj),proj_on_surf_rh,'subject_id')

    # ================================================================ #
    ### Collect all single subject datasets into one ##
    joiner = pe.JoinNode(name = 'joiner'+name, interface=Function(
            function     = DSjoiner,
            input_names  = ['subjMapperDS'],
            output_names = ['DSout']),
            joinsource   = 'info_subjectsMain',
            joinfield    = 'subjMapperDS',
            unique       = True)
    joinAllDS = joiner.clone('joinAllDS'+name)
    joinROIvoxelDS = joiner.clone('joinROIvoxelDS'+name)
    joinCohMask = joiner.clone('joinCohMask'+name)

    # collect all single subject datasets
    # somehow the getDS gives a list, instead of a string for the fileadress
    retMapFlow.connect(parser1,'in_file',joinAllDS,'subjMapperDS')
    # collect all single subject ROIvoxel datasets
    retMapFlow.connect(getDsROIvoxel,'DSout',joinROIvoxelDS,'subjMapperDS')
    # collect all roimasks (will be empty if cfg.RetMap.useCohMask is False)
    retMapFlow.connect(getDsROIvoxel,'roimask',joinCohMask,'subjMapperDS')

    # ================================================================ #
    ## give all the output ##
    outputspec = pe.Node(interface=util.IdentityInterface(
        fields=['DSout','ROIvoxel',
        'preprocScans','outdir','phaseFiles','surfacefiles_R','surfacefiles_L',
        'coreg2anat','roimasks']
        ), name='outputspec')
    retMapFlow.connect(joinAllDS,'DSout',outputspec,'DSout')
    retMapFlow.connect(joinROIvoxelDS,'DSout', outputspec, 'ROIvoxel')
    retMapFlow.connect(joinCohMask,'DSout', outputspec, 'roimasks')
    # the projection stuff
    retMapFlow.connect(getDSandPhase,'preprocScans',outputspec,'preprocScans')
    retMapFlow.connect(getDSandPhase,'outdir',outputspec,'outdir')
    retMapFlow.connect(getDSandPhase,'phaseFiles',outputspec,'phaseFiles')
    retMapFlow.connect(proj_on_surf_rh,'out_file',outputspec,'surfacefiles_R')
    retMapFlow.connect(proj_on_surf_lh,'out_file',outputspec,'surfacefiles_L')
    retMapFlow.connect(coregPh2FS,'out_file',outputspec,'coreg2anat')
    return retMapFlow

########################################################################
# Movie flow: SingleSubject - PreProc
########################################################################
def create_moviePPflow(
                        cfgFile,
                        condition="Movie",
                        name='moviePreProcFlow'):
    """ Preprocess and collect all movie data. The loaded cfg-file tells
    which processing steps should be done or not
    
    Parameters
    ----------

    name    : name of workflow (default: moviePreProcFlow)
    cfgFile : path to used config file (i.e. "/AbsPathToFile/hypal.cfg"
    
    Inputs::
    inputspec.nscans : number of scans to use/preproc from mapperDS
    inputspec.subj   : subj to process
    roimasks : nifti to mask the DS, thus to generate ROIvoxel
    Outputs::
    outputspec.DSout    : complete preprocessed and collected movie dataset
    outputspec.ROIvoxel : DS dataset restriced to only ROIvoxel

    """
    moviePreProcFlow = pe.Workflow(name=name)

    from nipype.interfaces.freesurfer import MRIConvert
    from nipype.interfaces import fsl
    from nipype.interfaces import afni
    from HA_extFunctions import (
        collect_MovieDS, removeGrandMean,
        restrictDs2ROIvoxel, DSjoiner, parserSingle, parserSmoother,
        listsel,
        )    
    # ================================================================ #
    # Infosource
    inputspec = pe.Node(interface=util.IdentityInterface(
        fields=['nscans','subj','roimasks']), name='inputspec')   

    # ---------------------------------------------------------------- #
    ## load cfgFile to get all the settings ##
    # (better then passing them via inputspec(!?))
    #
    print "..loading config from: ", cfgFile
    cfg = configLoader(cfgFile)
    if condition == "Movie": cond = cfg.Movie
    elif condition == "Other": cond = cfg.Other
    print "..done"
    
    ## define workflow specific input
    print "> ",name," -- \n temporalSmooth: \t", cond.tempSmooth, \
            "\n condition: \t\t",condition, \
            "\n fixTR: \t\t",cond.fixTR,  \
            "\n spatialSmooth: \t",cond.spatialSm, \
            "\n despike: \t\t",cond.despike
    
    # ================================================================ #
    # get data to train mapper (hypal) #
    # generate list depedning on nscans to use datagrabber depending on nscans
    # without connection to nscans it gets uncoupled from the nscans dependency,
    # thus it generates everything for all subjects
    def getscanlist(nscans):
        scanlist = []
        [(scanlist.append("%03d" % scan)) for scan in range(1,nscans+1)]
        return scanlist

    datagrabber = Node(SelectFiles(cond.templates), name="datagrabber"+name)
    datagrabber.inputs.base_directory = cfg.datadir
    datagrabber.inputs.force_lists    = True
    datagrabber.inputs.sort_filelist  = True

    moviePreProcFlow.connect(inputspec,'subj',datagrabber,'subj')
    moviePreProcFlow.connect(inputspec,('nscans',getscanlist),datagrabber,'scanlist')
    
    ## placeholder for skipping preproc steps 
    skipper = pe.MapNode(name='skipper'+name,interface=Function(
            function = parserSingle,
            input_names  = ['in_file'],
            output_names = ['out_file']),
            iterfield=['in_file'])
    # ---------------------------------------------------------------- #
    ## preproc Movie data ##
    # 
    if cond.fixTR is True:
        fixMovieTR = pe.MapNode(name="fixMovieTR"+name, interface=MRIConvert(
                                args = "-tr " + str(cond.TR*1000.),
                                out_type = "niigz",
                                subjects_dir = cfg.FSdir),
                                iterfield=['in_file'])
    else: 
        print "..skip fix TR"
        fixMovieTR = skipper.clone("fixMovieTR"+name)
    moviePreProcFlow.connect(datagrabber,'func',fixMovieTR,'in_file')

    if cond.spatialSm is True:
        spatialSmooth = pe.MapNode(name="spatialSmooth"+name,interface=fsl.SUSAN(
                                brightness_threshold = cond.brightness_threshold,
                                fwhm                 = cond.fwhm),
                                iterfield=['in_file'])
    else:
        print "..skip spatial Smooth"
        spatialSmooth = pe.MapNode(name="spatialSmooth"+name,
                interface=Function(function = parserSmoother,
                input_names  = ['in_file'], output_names = ['smoothed_file']),
                iterfield=['in_file'])
    moviePreProcFlow.connect(fixMovieTR,'out_file',spatialSmooth,'in_file')

    if cond.tempSmooth is True:
        temporalSmooth = pe.MapNode(name="temporalSmooth"+name,interface=afni.Bandpass(
                                lowpass = 1./cond.lowpass,
                                highpass = 1./cond.highpass,
                                outputtype = 'NIFTI_GZ'),
                                iterfield=['in_file'])
    else: 
        print "..skip temporal Smooth"        
        temporalSmooth = skipper.clone("temporalSmooth"+name)
    moviePreProcFlow.connect(spatialSmooth,'smoothed_file',temporalSmooth,'in_file')
    
    if cond.despike is True: 
        despike = pe.MapNode(name="despike"+name,interface=afni.Despike(
                                outputtype = 'NIFTI_GZ'),
                                iterfield=['in_file'])
    else: despike = skipper.clone("despike"+name)
    moviePreProcFlow.connect(temporalSmooth,'out_file',despike,'in_file') 

    # ---------------------------------------------------------------- #
    # 1st get data for each single subject
    getDS = pe.Node(name = 'getDS'+name, interface=Function(
            function     = collect_MovieDS,
            input_names  = ['subj','nscans','inFilesList','ROImaskL',
                            'ROImaskR','brainmask','leftGM','rightGM','cutAtEdge','roimask','rl'],
            output_names = ['DSout']),
            )
    getDS.inputs.cutAtEdge  = cond.cutAtEdge
    getDS.inputs.rl         = cond.rl
    moviePreProcFlow.connect(despike,'out_file', getDS,'inFilesList')
    moviePreProcFlow.connect(inputspec,'subj',getDS,'subj')
    moviePreProcFlow.connect(inputspec,'nscans',getDS,'nscans')
    moviePreProcFlow.connect(datagrabber,'brainmask',getDS,'brainmask')
    moviePreProcFlow.connect(datagrabber,'rightROI',getDS,'ROImaskR')
    moviePreProcFlow.connect(datagrabber,'leftROI',getDS,'ROImaskL')    
    moviePreProcFlow.connect(datagrabber,'rightGM',getDS,'rightGM')
    moviePreProcFlow.connect(datagrabber,'leftGM',getDS,'leftGM')
    # if (coh of phasemaps should be used for masking the toProject dataset)
    # AND (is provided by a give retmap dataset):
    # (is empty anyway if useCohMask is False)
    if cfg.RetMap.useCohMask and ((cfg.ProjMapper.MapperDS == "RetMap") or (cfg.ProjMapper.ProjDS == "RetMap")):
        roimaskSel = pe.Node(name = 'roimaskSel'+name, interface=Function(
            function     = listsel,
            input_names  = ['inlist','pos'],
            output_names = ['out']),
            )
        moviePreProcFlow.connect(inputspec,'subj',roimaskSel,'pos')
        moviePreProcFlow.connect(inputspec,'roimasks',roimaskSel,'inlist')
        moviePreProcFlow.connect(roimaskSel,'out',getDS,'roimask')

    ## remove average of all voxels tSeries
    removeMean = pe.Node(name='removeMean'+name,
                        interface=Function(
                            function     = removeGrandMean,
                            input_names  = ['infileSave','doit'],
                            output_names = ['DSout']))
    removeMean.inputs.doit = cond.removeGrandMean
    moviePreProcFlow.connect(getDS,'DSout',removeMean,'infileSave')

    ## mask DS with fa from DS and save them as ROIvoxel. ##
    getDsROIvoxel = pe.Node(name='getDsROIvoxel'+name,interface=Function(
            function     = restrictDs2ROIvoxel,
            input_names  = ['dsAll','noROI','useCohMask','cohThres'],
            output_names = ['DSout','roimask','makeROIdsNifti']),
            )
    getDsROIvoxel.inputs.noROI = cfg.noROI
    getDsROIvoxel.inputs.makeROIdsNifti = False
    # if retmap is used with a coh threshold, this also needs to be 
    # used for the movieDS-roiVoxel selection. Otherwise the Hypal doesn't work;
    # For Hypel the mapperTrainingDS (aka movieDS) and the projDS (aka retmapROIvoxelDS)
    # need to have a dimension in common, here the nROIvoxel.
    if (cfg.ProjMapper.MapperDS == "RetMap") or (cfg.ProjMapper.ProjDS == "RetMap"):
        getDsROIvoxel.inputs.cohThres   = cfg.RetMap.cohThres
        getDsROIvoxel.inputs.useCohMask = cfg.RetMap.useCohMask
        
    moviePreProcFlow.connect(removeMean,'DSout',getDsROIvoxel,'dsAll')
    # ---------------------------------------------------------------- #
    ## Collect all single subject datasets into one ##
    joiner = pe.JoinNode(name = 'joiner'+name, interface=Function(
            function     = DSjoiner,
            input_names  = ['subjMapperDS'],
            output_names = ['DSout']),
            joinsource   = 'info_subjectsMain',
            joinfield    = 'subjMapperDS',
            unique       = True)
   
    # collect all single subject datasets
    joinAllDS = joiner.clone('joinAllDS'+name)
    moviePreProcFlow.connect(removeMean,'DSout',joinAllDS,'subjMapperDS')
    
    # collect all single subject ROIvoxel datasets
    joinROIvoxelDS = joiner.clone('joinROIvoxelDS'+name)
    moviePreProcFlow.connect(getDsROIvoxel,'DSout',joinROIvoxelDS,'subjMapperDS')

    # ================================================================ #
    ## give all the output ##
    outputspec = pe.Node(interface=util.IdentityInterface(
            fields=['DSout','ROIvoxel']
            ), name='outputspec')
    moviePreProcFlow.connect(joinAllDS,'DSout',outputspec,'DSout')
    moviePreProcFlow.connect(joinROIvoxelDS,'DSout',outputspec,'ROIvoxel')
    return moviePreProcFlow

#########################################################################
## Connectome flow: do Hyper,Hyper! ##
# [Scooter et al., 1994; https://www.youtube.com/watch?v=RHVSshgPlQs]
#########################################################################
def create_doHyperHyper(cfgFile,name='hyperhyperflow'):
    """ Get PreConnectome (Connectome from Hyperaligned Spheres). 
    
    1.) hyperalign spheres of whole brain,
       'getMapperDSsave','ProjDS' need to be the sphere tSeries 
       
    2.) get CS and Mappers to project sphere tSeries into thier 
        subjectspecific own Commonspace (CSsubjOut)
   
    3.) get Connectome with spheres in subjectspecific CS and ROIvoxel ##

    4.) use the Conenctome to get Mapper and project ProjDS
     --> done in postHypal-workflow via main function
    
    Parameters
    ----------

    name    : name of workflow (default: hyperhyperflow)  [Scooter et al., 1994, ... and the Beat Goes On!]
    cfgFile : path to used config file (i.e. "/AbsPathToFile/hypal.cfg"

    Inputs::
    inputspec.sparse_radius :
        radius of spheres. Hyperalign these spheres into thier 
        subjectspecific Commonspace (CS). Use this CS get Connectome
        with ProjDS  for Hyperalignment.
    inputspec.Mapper_refDS  : 
        Dataset to train Mapper with

    Outputs::
    outputspec.sphereCommonspace              :
        Spheres projected into thier subject specific Commonspace.
        Might be useful to Compare quality of first hyperalginmetn step.        
    outputspec.Connectome_ProjDSvsHypalSphere : 
        Connectome out of "Hyperaligned spheres in thier CS" AND
        "ProjDS tSeries"
           
    
    """
    hyperhyperflow = pe.Workflow(name=name)
    
    from HA_extFunctions import (
        getMapperANDCommonspace,
        make_Connectome,
        getSphereMeanDS,
        )
    
    # ================================================================ #
    ## define workflow specific input
    
    # Infosource
    inputspec = pe.Node(interface=util.IdentityInterface(
        fields=['sparse_radius','Mapper_refDS','Mapper_ROIvoxel'
                ]), name='inputspec')   
    # ---------------------------------------------------------------- #
    ## load cfgFile to get all the settings ##
    # (better then passing them via inputspec(!?))
    #
    print "..loading config from: ", cfgFile
    cfg = configLoader(cfgFile)
    print "..done"
    # ================================================================ #

    
    ### hyperalign spheres of whole brain,
    #
    ## 'getMapperDSsave','ProjDS' need to be the sphere tSeries ##
    getSphereMeanTS = pe.Node(name='getSphereMeanTS',interface=Function(
        function     = getSphereMeanDS,
        input_names  = ['refDSsave','sparse_radius','nproc'],
        output_names = ['SeedmeanOut']))
    getSphereMeanTS.inputs.nproc = cfg.nproc
    hyperhyperflow.connect(inputspec,'sparse_radius',getSphereMeanTS,'sparse_radius')
    hyperhyperflow.connect(inputspec,'Mapper_refDS',getSphereMeanTS,'refDSsave')
    
    ## get CS and Mappers to project sphere tSeries into thier ##
    # subjectspecific own Commonspace (CSsubjOut) #
    preHypergetMapperCS = pe.MapNode(name='preHypergetMapperCS',interface=Function(
            function     = getMapperANDCommonspace,
            input_names  = ['getMapperDSsave','ProjDS','toSubj','usePh','preConn'],
            output_names = ['MapperOut','CommonspaceOut','CSsubjOut']),
            iterfield    = ['toSubj'])
    preHypergetMapperCS.inputs.toSubj  = cfg.toSubj4Mapper_list
    preHypergetMapperCS.inputs.usePh   = cfg.ProjDS.usePhInsteadTseries
    preHypergetMapperCS.inputs.preConn = cfg.usePreConnectome
    hyperhyperflow.connect(getSphereMeanTS,'SeedmeanOut',preHypergetMapperCS,'ProjDS')
    hyperhyperflow.connect(getSphereMeanTS,'SeedmeanOut',preHypergetMapperCS,'getMapperDSsave')
    
    ## get Connectome with spheres in subjectspecific CS and ROIvoxel ##
    getConnectome = pe.Node(name='getConnectome',interface=Function(
        function     = make_Connectome,
        input_names  = ['sparse_radius','ROIvoxel','refDS','brainmask',
                        'refDSisMean','nproc','connectomeMethod'],
        output_names = ['ConnectomeOut','intersectionMaskPath']))
    # if using ConnectomeCStSeries, no further spheres needed, thus
    getConnectome.inputs.sparse_radius = -1
    getConnectome.inputs.refDSisMean   = True
    getConnectome.inputs.nproc         = cfg.nproc
    hyperhyperflow.connect(inputspec,'Mapper_ROIvoxel',getConnectome,'ROIvoxel')
    hyperhyperflow.connect(preHypergetMapperCS,'CSsubjOut',getConnectome,'refDS')
    #
    ###

    ###  use the Conenctome to get Mapper and project ProjDS
    # --> done in postHypal-workflow via main function

    # ================================================================ #
    ## give all the output ##
    outputspec = pe.Node(interface=util.IdentityInterface(
        fields=['sphereCommonspace','Connectome_ProjDSvsHypalSphere']
        ), name='outputspec')
    hyperhyperflow.connect(preHypergetMapperCS,'CSsubjOut',outputspec,'sphereCommonspace')
    hyperhyperflow.connect(getConnectome,'ConnectomeOut',outputspec,'Connectome_ProjDSvsHypalSphere')
    return hyperhyperflow

########################################################################
# get Commonspace , Backprojection and Comparisons
########################################################################
def create_postHypalAnalysis(cfgFile,name='postHypalAnalysis'):
    """ Get the backprojected data, do correlation inside CS to get
        quality of Hypal and compare with anatalignment. 
    1.) getMapperCS: Get Mapper by using the InputConnectome 
        
    2.) getbackproj: Get backprjection into subjectspace by Hypal and
        and also correlate tSeries in CS
    3.) anatVShyperCorr: compare Hypal and AnatAlignment
    
    Parameters
    ----------

    name    : name of workflow (default: preConnFlow)
    cfgFile : path to used config file (i.e. "/AbsPathToFile/hypal.cfg"

    Inputs::
    inputspec.getMapperDSsave  : Connectome or MapperDS to get Mapper to project to CS.
                              If Connectome or Mapper is given, depends on useConnectome=True/False
    inputspec.ProjROIvoxel  : ProjDS restricted to ROIvoxel. Project these voxel 
                              by Hypal and compare with AnatAl
    inputspec.nscans        : number scans to use to get Connectome (nscans of MapperDS)
    inputspec.ProjDSout     : Whole ProjDS. It will be restriced to a allSubj-anatomical
                              intersection mask to generate a "new" ProjROIvoxelDS
                              to take as basis for the comparison of anatomical 
                              alignment with Hyperalignment. 
                              Because anatomical Alignment (AnatAl) needs
                              not only intersubject, interscan intersection
                              the new anatomical Intersectionmask is necessary :(

    Outputs::
    outputspec.AnatAlign2HyperOut : The output is not really needed. All stuff is already
                                    saved in the h5dir. 
    outputspec.Backproj2toSubjOut : Backprojected into subjects space
    outputspec.CScorrToSubjOut    : Correlation of the Commonspace data to the target subject
    outputspec.CompToSubjOut      : Correlation of the target subj data with the projected data output
    outputspec.CompToSubjOut_data : data used to get CompToSubjOut
    outputspec.Mapper             : the matrix used to project/transform the data
    outputspec.CSsubj             : the subjects commonspace
    outputspec.outnifti           : the projected data converted as nifti,
                including some info which scan from which to which subject
    outputspec.outnifti_raw       : as outnifti but without further
                information (I know, this is a mess)
    outputspec.surfacefiles_R     : outnifti as surface files for right hemis
    outputspec.surfacefiles_L     : .. and for left hemis
    outputspec.coreg2anat         : coregistered outfiles
    outputspec.meanCS_Backproj2toSubjOut : FOR ALL FOLLOWING:
            same as the others without meanCS, except here it is not single subject
            to subject transformation. Here teh data is averaged in the 
            commonspace and then transfered/projected.
            Doing this in style of Gubtupalli et al., 2016.
            It significantly invreases the quality of the hyperalignment 
            output.
    outputspec.meanCS_CScorrToSubjOut    : see above
    outputspec.meanCS_CompToSubjOut      : see above
    outputspec.meanCS_outnifti           : see above
    outputspec.meanCS_outnifti_raw       : see above
    outputspec.meanCS_CompToSubjOut_data : see above

   
    """
    postHypalAnalysis = pe.Workflow(name=name)
    
    from HA_extFunctions import (
        getMapperANDCommonspace,
        leave1out_proj,
        CompareAnatAlign2HyperAlign,
        )
    import nipype.interfaces.freesurfer as fs
    from nipype.interfaces import fsl
    def listsel(inlist): return inlist[0]

    # ================================================================ #
    ## define workflow specific input
    
    # Infosource
    inputspec = pe.Node(interface=util.IdentityInterface(
        fields=['cfgFile','ProjROIvoxel',
                'getMapperDSsave','nscans'
                ]), name='inputspec')   
    # ---------------------------------------------------------------- #
    ## load cfgFile to get all the settings ##
    # (better then passing them via inputspec(!?))
    #
    print "..loading config from: ", cfgFile
    cfg = configLoader(cfgFile)
    print "..done"
    # ================================================================ #
    ## get CS and Mappers ##
    getMapperCS = pe.MapNode(name='getMapperCS',interface=Function(
            function     = getMapperANDCommonspace,
            input_names  = ['getMapperDSsave','ProjDS','toSubj','usePh','preConn'],
            output_names = ['MapperOut','CommonspaceOut','CSsubjOut']),
            iterfield    = ['toSubj'])
    getMapperCS.inputs.toSubj = cfg.toSubj4Mapper_list
    getMapperCS.inputs.usePh  = cfg.ProjDS.usePhInsteadTseries
    # if noConn, use reshaping mechanism as for preConn usage
    if cfg.connectomeKey is 'noConn': getMapperCS.inputs.preConn = True
    postHypalAnalysis.connect(inputspec,'ProjROIvoxel',getMapperCS,'ProjDS')
    postHypalAnalysis.connect(inputspec,'getMapperDSsave',getMapperCS,'getMapperDSsave')

    ## project data back to according subjectspace ##
    # run it with leaveoutscans depending on nscans (here 4 to 8 nscans)
    # ---------------------------------------------------------------- #
    # if RetMap is projected AND only the phase values (NOT tSeries),
    # only the phaseValues need to be backprojected. Thus just "nChannels" Volumes
    # per condition (efd,pfd). nChannels, cause each phase deg value gets
    # converted in nChannels
    if cfg.ProjDS.usePhInsteadTseries: nVol = 2
    else: nVol = cfg.ProjDS.nVol
    # ---------------------------------------------------------------- #
    getbackproj = pe.MapNode(name='getbackproj',interface=Function(
            function     = leave1out_proj,
            input_names  = ['Mapper','CS','ProjDS', 'toSubj',
                            'nVol','makenifti','nscansMapper',
                            ],
            output_names = ['Backproj2toSubjOut','CScorrToSubjOut',
                            'CompToSubjOut','outnifti','outnifti_raw',
                            'CompToSubjOut_data']),
            iterfield    = ['toSubj'])                                                                                              
    getbackproj.inputs.toSubj    = cfg.toSubj_list
    getbackproj.inputs.makenifti = cfg.makenifti
    getbackproj.inputs.nVol      = nVol
    postHypalAnalysis.connect(inputspec,'nscans',getbackproj,'nscansMapper')
    postHypalAnalysis.connect(inputspec,'ProjROIvoxel',getbackproj,'ProjDS')
    postHypalAnalysis.connect(getMapperCS,'MapperOut',getbackproj,'Mapper')
    postHypalAnalysis.connect(getMapperCS,'CSsubjOut',getbackproj,'CS')
    
    getbackproj_meanCS = pe.MapNode(name='getbackproj_meanCS',interface=Function(
            function     = leave1out_proj,
            input_names  = ['Mapper','CS','ProjDS', 'toSubj',
                            'nVol','makenifti','nscansMapper','makeMeanCS'
                            ],
            output_names = ['Backproj2toSubjOut','CScorrToSubjOut',
                            'CompToSubjOut','outnifti','outnifti_raw',
                            'CompToSubjOut_data']),
            iterfield    = ['toSubj'])                                                                                              
    getbackproj_meanCS.inputs.toSubj     = cfg.toSubj_list
    getbackproj_meanCS.inputs.makenifti  = cfg.makenifti
    getbackproj_meanCS.inputs.nVol       = nVol
    getbackproj_meanCS.inputs.makeMeanCS = True
    postHypalAnalysis.connect(inputspec,'nscans',getbackproj_meanCS,'nscansMapper')
    postHypalAnalysis.connect(inputspec,'ProjROIvoxel',getbackproj_meanCS,'ProjDS')
    postHypalAnalysis.connect(getMapperCS,'MapperOut',getbackproj_meanCS,'Mapper')
    postHypalAnalysis.connect(getMapperCS,'CSsubjOut',getbackproj_meanCS,'CS')

    ## do correlation of before and after HYPAL datasets to get baseline 
    # for alignment quality etc.
    # get the comparison AnatAlignVSHyperal
    anatVShyperCorr = pe.Node(name='anatVShyperCorr',interface=Function(
            function     = CompareAnatAlign2HyperAlign,
            input_names  = ['ProjDSsave','Backproj2toSubj','usePhInsteadTseries'],
            output_names = ['DSout']),
            )
    anatVShyperCorr.inputs.usePhInsteadTseries = cfg.ProjDS.usePhInsteadTseries
    postHypalAnalysis.connect(inputspec,'ProjROIvoxel',anatVShyperCorr,'ProjDSsave')
    postHypalAnalysis.connect(getbackproj,'Backproj2toSubjOut',anatVShyperCorr,'Backproj2toSubj')

    # ================================================================ #
    ## give all the output ##
    outputspec = pe.Node(interface=util.IdentityInterface(
            fields=['Backproj2toSubjOut','CScorrToSubjOut',
                    'CompToSubjOut','CompToSubjOut_data','AnatAlign2HyperOut','Mapper','CSsubj',
                    'outnifti','outnifti_raw','surfacefiles_R','surfacefiles_L','coreg2anat',
                    'meanCS_Backproj2toSubjOut','meanCS_CScorrToSubjOut','meanCS_CompToSubjOut',
                    'meanCS_outnifti','meanCS_outnifti_raw','meanCS_CompToSubjOut_data']
            ), name='outputspec')
    postHypalAnalysis.connect(getbackproj,'Backproj2toSubjOut',outputspec,'Backproj2toSubjOut')
    postHypalAnalysis.connect(getbackproj,'CScorrToSubjOut',outputspec,'CScorrToSubjOut')
    postHypalAnalysis.connect(getbackproj,'CompToSubjOut',outputspec,'CompToSubjOut')
    postHypalAnalysis.connect(getbackproj,'CompToSubjOut_data',outputspec,'CompToSubjOut_data')
    postHypalAnalysis.connect(getbackproj,'outnifti',outputspec,'outnifti')
    postHypalAnalysis.connect(getbackproj,'outnifti_raw',outputspec,'outnifti_raw')
    postHypalAnalysis.connect(getbackproj_meanCS,'Backproj2toSubjOut',outputspec,'meanCS_Backproj2toSubjOut')
    postHypalAnalysis.connect(getbackproj_meanCS,'CScorrToSubjOut',outputspec,'meanCS_CScorrToSubjOut')
    postHypalAnalysis.connect(getbackproj_meanCS,'CompToSubjOut',outputspec,'meanCS_CompToSubjOut')
    postHypalAnalysis.connect(getbackproj_meanCS,'CompToSubjOut_data',outputspec,'meanCS_CompToSubjOut_data')
    postHypalAnalysis.connect(getbackproj_meanCS,'outnifti',outputspec,'meanCS_outnifti')
    postHypalAnalysis.connect(getbackproj_meanCS,'outnifti_raw',outputspec,'meanCS_outnifti_raw')
    postHypalAnalysis.connect(getMapperCS,'MapperOut',outputspec,'Mapper')
    postHypalAnalysis.connect(getMapperCS,'CSsubjOut',outputspec,'CSsubj')
    postHypalAnalysis.connect(anatVShyperCorr,'DSout',outputspec,'AnatAlign2HyperOut')

    return postHypalAnalysis

########################################################################
# get freesurfer overlays of backprojected nifti
########################################################################
def create_makeOutniftilistOverlays(cfgFile,name='getOutniftiOverlays',useGrpTmpl=False):

    """ 
    Just coregistration and surface projection of input niftis.
    
    
    Parameters
    ----------

    name    : name of workflow (default: preConnFlow)
    cfgFile : path to used config file (i.e. "/AbsPathToFile/hypal.cfg"

    Inputs::
    subj     : subject to process
    outnifti : nifti to overlay

    Outputs::
    surfacefiles_R : surface files for right hemisphere
    surfacefiles_L : surface files for right hemisphere
    coreg2anat     : coregistered niftis
   
    """
    getOutniftiOverlays = pe.Workflow(name=name)
    
    from HA_extFunctions import (
        niftiselector,
        )
    import nipype.interfaces.freesurfer as fs
    from nipype.interfaces import fsl
    def listsel(inlist): return inlist[0]
    def sub2subj(sub): return "sub-"+sub
    
    # ================================================================ #
    ## define workflow specific input
    
    # Infosource
    inputspec = pe.Node(interface=util.IdentityInterface(
        fields=['cfgFile','subj','outnifti',
                ]), name='inputspec')   
    # ---------------------------------------------------------------- #
    ## load cfgFile to get all the settings ##
    # (better then passing them via inputspec(!?))
    #
    print "..loading config from: ", cfgFile
    cfg = configLoader(cfgFile)
    print "..done"
    
    # ================================================================ #
    # generate list depedning on nscans to use datagrabber depending on nscans
    # without connection to nscans it gets uncoupled from the nscans dependency,
    # thus it generates everything for all subjects
    def getscanlist(nscans):
        scanlist = []
        [(scanlist.append("%03d" % scan)) for scan in range(1,nscans+1)]
        return scanlist
    if useGrpTmpl == False: datagrabber = Node(SelectFiles(cfg.RetMap.templates), name="datagrabber"+name)
    else: datagrabber = Node(SelectFiles(cfg.RetMap_grpTmpl.templates), name="datagrabber"+name)
    datagrabber.inputs.force_lists    = True
    datagrabber.inputs.sort_filelist  = True
    datagrabber.inputs.base_directory = cfg.datadir
    datagrabber.inputs.scanlist       = getscanlist(cfg.RetMap.nscans)
    getOutniftiOverlays.connect(inputspec,'subj',datagrabber,'subj')

    # ================================================================ #
    ### project the backprojected retmaps onto surface ###
    #

    ## prepare outnifti list for further steps #
    # (extract filenames for according fromSubj and toSubj from summary list "outnifti")
    niftiselector = pe.Node(name='niftiselector'+name, interface=Function(
            function     = niftiselector,
            input_names  = ['outnifti','subj','subjlist'],
            output_names = ['toSubjNiftilist']))
    niftiselector.inputs.subjlist = cfg.subjlist
    getOutniftiOverlays.connect(inputspec,'outnifti',niftiselector,'outnifti')
    getOutniftiOverlays.connect(inputspec,'subj',niftiselector,'subj')

    ### coregister phasemap to FS anatomy to project on surface ###
    #
    coreg2FS = pe.MapNode(name = 'coreg2FS'+name,
                            interface=fsl.FLIRT(apply_xfm=True),
                            iterfield=['in_file'])    
    getOutniftiOverlays.connect(datagrabber,('anatref',listsel),coreg2FS,'reference')
    getOutniftiOverlays.connect(datagrabber,('regfile',listsel),coreg2FS,'in_matrix_file')
    
    # for grptmpl, first warp image back from grpspace to subjspace,
    # then applyxfm (in postmat)
    if useGrpTmpl:
        warpGrp2subj = pe.MapNode(name = 'warpGrp2subj'+name,
                                interface=fsl.ApplyWarp(),
                                iterfield=['in_file'])    
        getOutniftiOverlays.connect(niftiselector,'toSubjNiftilist',warpGrp2subj,'in_file')
        getOutniftiOverlays.connect(datagrabber,('warpGrp2subj',listsel),warpGrp2subj,'field_file')
        getOutniftiOverlays.connect(datagrabber,('refbold3Tp2',listsel),warpGrp2subj,'ref_file')
        getOutniftiOverlays.connect(warpGrp2subj,'out_file',coreg2FS,'in_file')

    else:
        # just use applyxfm
        getOutniftiOverlays.connect(niftiselector,'toSubjNiftilist',coreg2FS,'in_file')
    #
    ####

    ## project on surface ##
    # only project, no further coreg
    proj_on_surf_lh = pe.MapNode(name = 'proj_on_surf_lh_'+name,
                        interface=fs.SampleToSurface(
                            hemi='lh',
                            subjects_dir=cfg.FSdir,
                            sampling_method="average",
                            sampling_range=(0.0,0.75,0.001),
                            sampling_units="frac",
                            surface = "white",
                            #projection_stem="",
                            reg_header=True
                            ),
                            iterfield=['source_file'])
    getOutniftiOverlays.connect(coreg2FS,'out_file',proj_on_surf_lh,'source_file')
    #getOutniftiOverlays.connect(niftiselector,'toSubjNiftilist',proj_on_surf_lh,'source_file')
    getOutniftiOverlays.connect(inputspec,('subj',sub2subj),proj_on_surf_lh,'subject_id')
 
    proj_on_surf_rh = pe.MapNode(name = 'proj_on_surf_rh_'+name,
                        interface=fs.SampleToSurface(
                            hemi='rh',
                            subjects_dir=cfg.FSdir,
                            sampling_method="average",
                            sampling_range=(0.0,0.75,0.001),
                            sampling_units="frac",
                            surface = "white",
                            #projection_stem="",
                            reg_header=True
                            ),
                            iterfield=['source_file'])
    getOutniftiOverlays.connect(coreg2FS,'out_file',proj_on_surf_rh,'source_file')
    #getOutniftiOverlays.connect(niftiselector,'toSubjNiftilist',proj_on_surf_rh,'source_file')
    getOutniftiOverlays.connect(inputspec,('subj',sub2subj),proj_on_surf_rh,'subject_id')

    # ================================================================ #
    
    ## give all the output ##
    outputspec = pe.Node(interface=util.IdentityInterface(
            fields=['surfacefiles_R','surfacefiles_L','coreg2anat']
            ), name='outputspec')    
    # the projection stuff
    getOutniftiOverlays.connect(proj_on_surf_rh,'out_file',outputspec,'surfacefiles_R')
    getOutniftiOverlays.connect(proj_on_surf_lh,'out_file',outputspec,'surfacefiles_L')
    #getOutniftiOverlays.connect(niftiselector,'toSubjNiftilist',outputspec,'coreg2anat')
    getOutniftiOverlays.connect(coreg2FS,'out_file',outputspec,'coreg2anat')

    return getOutniftiOverlays
