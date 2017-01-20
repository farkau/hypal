## FK1aug2014

#############################################
## See LICENSE for licensing and copyright ##
#############################################


def configLoader(cfgfile):
    """ loads the cfgfile and extraxts all the needed informations. 
    Then sends it back. """
    import ConfigParser
    from os.path import join
    from mvpa2.suite import h5save
    import time

    ####################################################################
    # init class to save everything inside #
    class cfg():
        def __init__(self, cfgfile):
            self.cfgfile = cfgfile
        class ProjDS():
            descr = "saves settings for processing the dataset that will be projected via hyperalignment"
        class MapperDS():
            descr = "saves settings for processing the dataset that will be used to generate the mapper for hyperalignment"
        class ProjMapper():
            descr = "here set which DS (ProjDS or MapperDS) to use to generate mapper and which to project onto other subjects"

    # ==================================================================== #
    ## load cfg file ##
    #
    #cfgfile = join('/home/fkaule/MRI/projects/Hyperalignment','hypal.cfg')
    config = ConfigParser.RawConfigParser()
    config.read(cfgfile)
    cfg.cfgfile = cfgfile


    ####################################################################
    ### initial dataset ###
    # these will saved according usage as MapperDS or ProjDS inside the 
    # cfg-class
    class RetMap():
        descr = "saves settings for RetMapping datasets"
    
    # RetMap data #
    RetMap.tasknr          = config.get('RetMap', 'tasknr')
    RetMap.TR              = config.getint('RetMap', 'TR')
    RetMap.nVol            = config.getint('RetMap', 'nVol')
    RetMap.delLead         = config.getint('RetMap', 'delLead')
    RetMap.retmapdir       = config.get('RetMap', 'retmapdir')
    RetMap.prefix          = config.get('RetMap', 'prefix')
    RetMap.brainmaskName   = config.get('RetMap', 'brainmaskName')
    RetMap.tsfilename      = config.get('RetMap', 'tsfilename')
    RetMap.cond_list       = config.get('RetMap', 'cond_list').split(',')
    RetMap.usePhInsteadTseries = config.getboolean('RetMap', 'usePh')
    RetMap.nscans          = config.getint('RetMap', 'nscans')    
    RetMap.useCohMask      = config.getboolean('RetMap', 'useCohMask')    
    RetMap.cohThres        = config.getfloat('RetMap', 'cohThres')
    RetMap.nChannels       = config.getint('RetMap', 'nChannels')
    
    class RetMap_noROI():
        descr = "settings for templates without usage of ROI"
    
    class RetMap_grpTmpl():
        descr = "settings for templeates from func in grpSpace aka grpTmpl"

    # ================================================================ #

    class Movie():
        descr = "saves settings for Movie datasets"
    
    # movie data #
    Movie.brainmaskName   = config.get('Movie', 'brainmaskName')
    Movie.tsfilename      = config.get('Movie', 'tsfilename')
    Movie.tasknr        = config.get('Movie', 'tasknr')
    Movie.nscans        = config.getint('Movie', 'nscans')
    Movie.TR            = config.getint('Movie', 'TR')
    Movie.delLead       = config.getint('Movie', 'delLead')
    Movie.cutAtEdge     = config.getint('Movie', 'cutAtEdge')
    Movie.highpass      = config.getint('Movie', 'highpass')
    Movie.lowpass       = config.getint('Movie', 'lowpass')
    Movie.brightness_threshold = config.getfloat('Movie', 'brightness_threshold')
    Movie.fwhm          = config.getint('Movie', 'fwhm')
    Movie.tempSmooth    = config.getboolean('Movie', 'temporalSmooth')
    Movie.spatialSm     = config.getboolean('Movie', 'spatialSmooth')
    Movie.despike       = config.getboolean('Movie', 'despike')
    Movie.fixTR         = config.getboolean('Movie', 'fixTR')
    Movie.removeGrandMean = config.getboolean('Movie', 'removeGrandMean')
    Movie.nVol          = -1 # set to Zero cause it varies across scans 
    Movie.rl            = config.getboolean('Movie', 'removeLast') 
    
    # ================================================================ #
    class Other():
        descr = "saves settings for Other datasets"
    
    # Other data #
    Other.brainmaskName   = config.get('Other', 'brainmaskName')
    Other.tsfilename      = config.get('Other', 'tsfilename')
    Other.tasknr        = config.get('Other', 'tasknr')
    Other.nscans        = config.getint('Other', 'nscans')
    Other.TR            = config.getint('Other', 'TR')
    Other.delLead       = config.getint('Other', 'delLead')
    Other.cutAtEdge     = config.getint('Other', 'cutAtEdge')
    Other.highpass      = config.getfloat('Other', 'highpass')
    Other.lowpass       = config.getfloat('Other', 'lowpass')
    Other.brightness_threshold = config.getfloat('Other', 'brightness_threshold')
    Other.fwhm          = config.getint('Other', 'fwhm')
    Other.tempSmooth    = config.getboolean('Other', 'temporalSmooth')
    Other.spatialSm     = config.getboolean('Other', 'spatialSmooth')
    Other.despike       = config.getboolean('Other', 'despike')
    Other.fixTR         = config.getboolean('Other', 'fixTR')
    Other.removeGrandMean = config.getboolean('Movie', 'removeGrandMean')
    Other.nVol          = -1 # set to Zero cause it varies across scans 
    Other.rl            = config.getboolean('Other', 'removeLast') 
    ####################################################################
    
    ## little helpers ##
    #
    def str2intList(strList):
        """ converts the list with strings inside into a list with in inside """
        intList = []
        [intList.append(int(ii)) for ii in strList.split(',')] 
        return intList

    # ==================================================================== #
    ## process cfg to get basic settings ##
    # and generate further settings out of it (i.e. paths) #
    #
    # major settings #
    cfg.nproc           = config.getint('major_settings', 'nproc')
    cfg.ProjMapper.ProjDS = config.get('major_settings', 'proj')
    cfg.ProjMapper.MapperDS = config.get('major_settings', 'mapper')
    cfg.makenifti           = config.getboolean('major_settings','makenifti')

    # ================================================================= #
    # load input which is needed to get further settings/paths.. #
    cfg.connectomeMethod  = config.get('analysis_variables', 'connectomeMethod')
    cfg.subjlist          = config.get('inputs', 'subjlist').split(',')
    cfg.skipSubj          = config.getboolean('analysis_variables', 'skipSubj')
    cfg.useConnectome     = config.getboolean('analysis_variables', 'useConnectome')
    cfg.usePreConnectome  = config.getboolean('analysis_variables', 'usePreConnectome')
    cfg.excludeROIvoxel   = config.getboolean('analysis_variables', 'excludeROIvoxel')
    cfg.roitmpl           = config.get('templates', 'roitmpl')
    cfg.gmtmpl            = config.get('templates', 'GMtmpl')
    cfg.onlyOneROI        = config.getboolean('templates','onlyOneROI')
    cfg.noROI             = config.getboolean('templates','noROI')
    cfg.datadir           = config.get('paths', 'datadir')
    cfg.workfolder        = config.get('paths', 'workfolder')
    cfg.hypalDir          = config.get('paths', 'hypalDir')
    cfg.sparse_radiusList = str2intList(config.get('analysis_variables', 'sparse_radiusList'))

    if config.has_option('templates', 'GM_grptmpl'):
        cfg.gmgrptmpl = config.get('templates', 'GM_grptmpl')
    else: 
        print('("templates", "GM_grptmpl") does not exist. Setting it to ""')
        cfg.gmgrptmpl = ""

    # if subjlist is 3 digits, shorten it to 2 digits (001 -> 01)
    if len(cfg.subjlist[0]) >2: cfg.subjlist = [su[1:3] for su in cfg.subjlist]

    # get further settings out of basics from hypal.cfg #
    cfg.nsubj             = len(cfg.subjlist)
    #cfg.nscans_list       = range(4,cfg.RetMap.nscans+1,4) # needs (nscans >= 4)
    cfg.toSubj_list       = range(cfg.nsubj) # keep minimal 2 subj

    if cfg.skipSubj is True:
        # the aim subj (toSubj) is left out from CS generation
        cfg.toSubj4Mapper_list  = cfg.toSubj_list
        cfg.skipSubjStr         = "_skipSubj2CS"
    else:
        # use all subj for CS
        cfg.toSubj4Mapper_list  = [-1]
        cfg.skipSubjStr         = "_allSubj2CS"

    # generate name of main dir depending on method to generate Connectome #
    if cfg.useConnectome is True:
        sRadID        = len(cfg.sparse_radiusList)*sum(cfg.sparse_radiusList)
        exclRoiKey        = "_exclRoi" + str(cfg.excludeROIvoxel)
        cfg.connectomeKey = "sRadList"+str(sRadID) + exclRoiKey
        cfg.connectomeKey = cfg.connectomeKey+"_preConn"+str(cfg.usePreConnectome)+\
                            "_" + cfg.connectomeMethod
    else:
        cfg.connectomeKey = "noConn"

    cfg.datadir         = config.get('paths', 'datadir')
    cfg.runTimeDate     = time.strftime("%c")

    # ================================================================= #
    ## define ProjDS and MapperDS dempending on which dataset shall be used: ##
    # "RetMap" or "Movie"
    if cfg.ProjMapper.ProjDS == "RetMap":
        cfg.ProjDS = RetMap
    elif cfg.ProjMapper.ProjDS == "Other":
        cfg.ProjDS = Other
    elif cfg.ProjMapper.ProjDS == "Movie":
        cfg.ProjDS = Movie
    
    if cfg.ProjMapper.MapperDS == "RetMap":
        cfg.MapperDS = RetMap
    elif cfg.ProjMapper.MapperDS == "Other":
        cfg.MapperDS = Other
    elif cfg.ProjMapper.MapperDS == "Movie":
        cfg.MapperDS = Movie
        
    ## anyway save the settings for Movie and RetMap data to use it
    #  in according processing flow
    cfg.Movie = Movie
    cfg.RetMap = RetMap
    cfg.Other = Other
    # only for templates
    cfg.RetMap_noROI = RetMap_noROI
    cfg.RetMap_grpTmpl = RetMap_grpTmpl 
    
    # ================================================================= #
    ## add templates etc. ##
    projMapperkey = "mapper"+cfg.MapperDS.tasknr+"proj"+cfg.ProjDS.tasknr    
    cfg.outdirkey   = cfg.connectomeKey + "_" + \
                        "nsubj" + str(cfg.nsubj) + \
                        cfg.skipSubjStr+"_"+ \
                        projMapperkey

    # dirs and paths #
    GM = "_"+cfg.gmtmpl    
    
    if cfg.ProjMapper.MapperDS == "RetMap": spatSm = ''
    elif cfg.MapperDS.spatialSm: spatSm = '_spatSm'+str(cfg.MapperDS.brightness_threshold)+'Bt'+str(cfg.MapperDS.fwhm)+'fwhm'
    else: spatSm = '_NOspatSm'
    
    if cfg.ProjMapper.MapperDS == "RetMap": tempSm = ''
    elif cfg.MapperDS.tempSmooth: tempSm = '_tempSm'+str(cfg.MapperDS.highpass)+'H'+str(cfg.MapperDS.lowpass)+'L'
    else: tempSm = '_NOtempSm'
    
    if cfg.ProjMapper.MapperDS == "RetMap":
        despi = ''
        rmGrM = ''
    else:
        despi = '_despi'+str(cfg.MapperDS.despike)
        rmGrM = '_rmGrM'+str(cfg.MapperDS.removeGrandMean)

    workdirkey = cfg.connectomeKey + GM + spatSm + tempSm + despi + rmGrM
    dirkey = "ROI"+cfg.roitmpl+ \
            '_mapper'+cfg.MapperDS.tasknr+ \
            '_Proj'+cfg.ProjDS.tasknr+ \
            workdirkey+ \
            '_'+config.get('analysis_settings', 'outdirAdd')
    cfg.retmapPath  = config.get('paths', 'retmapPath')
    cfg.FSdir       = config.get('paths', 'FSdir')
    cfg.workdir     = cfg.workfolder+'/wf_workdir_'+dirkey
    cfg.outdir      = cfg.hypalDir+'/wf_out_'+dirkey
    cfg.datasinkDir = join(cfg.outdir,cfg.outdirkey)
    cfg.niftidir    = join(cfg.outdir,'niftiout')
    cfg.tmpldir     = join(cfg.datadir,'templates','grpbold7Tad')
    cfg.brainmaskAll= config.get('templates', 'brainmaskAll') 

    # if only one ROI, no differentiation between left and right ROI
    if cfg.onlyOneROI == False:
        cfg.roiRname    = 'Right'+cfg.roitmpl
        cfg.roiLname    = 'Left'+cfg.roitmpl
    else:
        cfg.roiRname    = cfg.roitmpl
        cfg.roiLname    = cfg.roitmpl
    
    cfg.gmRname    = 'right'+cfg.gmtmpl+".nii.gz"
    cfg.gmLname    = 'left'+cfg.gmtmpl+".nii.gz" 
    cfg.gmgrpRname    = 'right'+cfg.gmgrptmpl+".nii.gz"
    cfg.gmgrpLname    = 'left'+cfg.gmgrptmpl+".nii.gz" 

    # mostly for the anatCdist
    RetMap_grpTmpl.templates = \
        {"func": "aligned/sub-{subj}/in_grpbold3Tp2/sub-{subj}_task-retmap_run-{scanlist!s}_bold.nii.gz",
        #"regfile": "tnt/sub-{subj}/bold3Tp2/in_t1w/xfm_6dof_occ.mat",
        "regfile": "tnt/sub-{subj}/bold3Tp2/in_t1w/xfm_6dof.mat",
        "rightROI": "aligned/sub-{subj}/in_grpbold3Tp2/ROI/"+cfg.roiRname+".nii.gz",
        "leftROI": "aligned/sub-{subj}/in_grpbold3Tp2/ROI/"+cfg.roiLname+".nii.gz",
        "rightGM": "aligned/sub-{subj}/in_grpbold3Tp2/ROI/"+cfg.gmRname,
        "leftGM": "aligned/sub-{subj}/in_grpbold3Tp2/ROI/"+cfg.gmLname,
        "anatref": "tnt/sub-{subj}/t1w/brain.nii.gz",
        "brainmask": "tnt/sub-{subj}/bold3Tp2/in_grpbold3Tp2/brain_mask.nii.gz",
        "warpGrp2subj": "tnt/sub-{subj}/bold3Tp2/in_grpbold3Tp2/tmpl2subj_warp.nii.gz",
        "refbold3Tp2": "tnt/sub-{subj}/bold3Tp2/brain.nii.gz",
        } 

    RetMap.templates = \
        {"func": "aligned/sub-{subj}/in_bold3Tp2/sub-{subj}_task-retmap_run-{scanlist!s}_bold.nii.gz",
        "regfile": "tnt/sub-{subj}/bold3Tp2/in_t1w/xfm_6dof.mat",
        "rightROI": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.roiRname+".nii.gz",
        "leftROI": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.roiLname+".nii.gz",
        "rightGM": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.gmRname,
        "leftGM": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.gmLname,
        "anatref": "tnt/sub-{subj}/t1w/brain.nii.gz",
        "brainmask": "tnt/sub-{subj}/bold3Tp2/brain_mask.nii.gz",
        }
        
    Movie.templates  = \
        {"func": "aligned/sub-{subj}/in_bold3Tp2/sub-{subj}_task-avmovie_run-{scanlist!s}_bold.nii.gz",
        "regfile": "tnt/sub-{subj}/bold3Tp2/in_t1w/xfm_6dof.mat",
        "rightROI": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.roiRname+".nii.gz",
        "leftROI": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.roiLname+".nii.gz",
        "rightGM": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.gmRname,
        "leftGM": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.gmLname,
        "anatref": "tnt/sub-{subj}/t1w/brain.nii.gz",
        "brainmask": "tnt/sub-{subj}/bold3Tp2/brain_mask.nii.gz",
        }

    Other.templates  = \
        {"func": "aligned/sub-{subj}/in_bold3Tp2/sub-{subj}_task-aomovie_run-{scanlist!s}_bold.nii.gz",
        "regfile": "tnt/sub-{subj}/bold3Tp2/in_t1w/xfm_6dof.mat",
        "rightROI": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.roiRname+".nii.gz",
        "leftROI": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.roiLname+".nii.gz",
        "rightGM": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.gmRname,
        "leftGM": "aligned/sub-{subj}/in_bold3Tp2/ROI/"+cfg.gmLname,
        "anatref": "tnt/sub-{subj}/t1w/brain.nii.gz",
        "brainmask": "tnt/sub-{subj}/bold3Tp2/brain_mask.nii.gz",
        }
 
    # overwrite brainmask depending on mapperDS: if 7T (task001) 
    # then different mask then 3T. This way the chance is much higher 
    # selecting the same voxel for the mapperROIvoxel and projROIvoxel
    if cfg.MapperDS.tasknr == '001': RetMap.templates['brainmask'] = Other.templates['brainmask']
 
    ## stop if nscans < 4
    #if (RetMap.usePhInsteadTseries == False) and (RetMap.nscans != 4):
        #print "ERROR.. choose RetMap nscans = 4"
        #from sys import exit
        #exit()
    return cfg
