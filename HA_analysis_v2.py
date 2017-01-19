# FK2014: Meta pipeline for Hyperalignment 

#############################################
## See LICENSE for licensing and copyright ##
#############################################
# ==================================================================== #
from HA_configLoader_v2 import configLoader
from HA_extFunctions import (
    make_Connectome,
    parser,
    get_refMapperDS,
        )
from HA_extWFlows import (
    create_rmFlow,
    create_moviePPflow,
    create_doHyperHyper,
    create_postHypalAnalysis,
    create_makeOutniftilistOverlays,
    )

from mvpa2.base import debug
import os 
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio           # Data i/o
from nipype.interfaces.utility import Function
from nipype import Node, Workflow, SelectFiles
import sys # to get the cfgFileName as input

# ==================================================================== #
### set Variables and paths ###
## set all the settings in hypal.cfg (the cfgFile) #
#
cfgFileName = sys.argv[1]
cfgFile = os.path.join('/home/fkaule/MRI/projects/Hyperalignment',cfgFileName)      
print "..loading config from: ", cfgFile
cfg = configLoader(cfgFile)
print "..done"

# set up workflow #
wf          = pe.Workflow(name= "wf_HYPALanalysis")
wf.base_dir = cfg.workdir
# generate output folder if not existing #
if not os.path.exists(cfg.datasinkDir): os.makedirs(cfg.datasinkDir)

# make backup of cfg-file #
from shutil import copyfile
cfgDestPath = os.path.join(cfg.datasinkDir,"cfgBackup.cfg")
copyfile(cfgFile,cfgDestPath) 

########################################################################
#### infosources, datagrabber etc. ####
#

# nscans to use for MapperDS. Depends on which DS used (RetMap or Movie)
info_nscans = pe.Node(name='info_nscans',interface=util.IdentityInterface(
        fields    = ['nscans']),
        iterables = [('nscans', ([cfg.MapperDS.nscans]))]
        )
info_subjects = pe.Node(name='info_subjectsMain',interface=util.IdentityInterface(
        fields    = ['subj']),
        iterables = [('subj', cfg.subjlist)])
########################################################################
#### Decissions on workflow settings ####
#

### Get Data ####
## preProc only if needed, thus first check if Movie and/or RetMap used ##
# and skip one if possible. Here just data generation and collection for 
# further steps:
if ((cfg.ProjMapper.ProjDS == "RetMap") or (cfg.ProjMapper.MapperDS == "RetMap")):
    print ".. using RetMapFlow"
    retMapFlow = create_rmFlow(cfgFile,name="retMapFlow")
    wf.connect(info_subjects,'subj',retMapFlow,'inputspec.subj')

if ((cfg.ProjMapper.ProjDS == "Movie") or (cfg.ProjMapper.MapperDS == "Movie")):
    print ".. using MovieFlow"
    movieFlow = create_moviePPflow(cfgFile,condition="Movie",name="movieFlow")
    wf.connect(info_nscans,'nscans',movieFlow,'inputspec.nscans')
    wf.connect(info_subjects,'subj',movieFlow,'inputspec.subj')
    # distribute the roimasks if generated AND needed somewhere else
    if ((cfg.ProjMapper.ProjDS == "RetMap") or (cfg.ProjMapper.MapperDS == "RetMap")):
        wf.connect(retMapFlow,'outputspec.roimasks',movieFlow,'inputspec.roimasks')
    else: movieFlow.inputspec.roimasks = ""
    
if ((cfg.ProjMapper.ProjDS == "Other") or (cfg.ProjMapper.MapperDS == "Other")):
    print ".. using MovieFlow for OtherDataset"
    otherFlow = create_moviePPflow(cfgFile,condition="Other",name="otherFlow")
    wf.connect(info_nscans,'nscans',otherFlow,'inputspec.nscans')
    wf.connect(info_subjects,'subj',otherFlow,'inputspec.subj')    
    # distribute the roimasks if generated AND needed somewhere else
    if ((cfg.ProjMapper.ProjDS == "RetMap") or (cfg.ProjMapper.MapperDS == "RetMap")):
        wf.connect(retMapFlow,'outputspec.roimasks',otherFlow,'inputspec.roimasks')
    else: movieFlow.inputspec.roimasks = ""

### which dataset should be projected ###
#
flowparser = pe.Node(name='flowparser',interface=Function(
    function     = parser,
    input_names  = ['in_DSout','in_ROIvoxel'],
    output_names = ['DSout_ROIvoxel_list']))
def listsel(inlist,pos): return inlist[pos]

# Get ProjDS from RetMapping: 
ProjSource = flowparser.clone("ProjSource")
if cfg.ProjMapper.ProjDS == "RetMap":
    wf.connect(retMapFlow,'outputspec.ROIvoxel',ProjSource,'in_ROIvoxel')
    wf.connect(retMapFlow,'outputspec.DSout',ProjSource,'in_DSout')
# Get ProjDS from Movie:
elif cfg.ProjMapper.ProjDS == "Movie":
    wf.connect(movieFlow,'outputspec.DSout',ProjSource,'in_DSout')
    wf.connect(movieFlow,'outputspec.ROIvoxel',ProjSource,'in_ROIvoxel')
# Get ProjDS from Other dataset:
elif cfg.ProjMapper.ProjDS == "Other":
    wf.connect(otherFlow,'outputspec.ROIvoxel',ProjSource,'in_ROIvoxel')
    wf.connect(otherFlow,'outputspec.DSout',ProjSource,'in_DSout')
     
print ">> MapperDS (dataset to train Hypal-mappers): ", cfg.ProjMapper.MapperDS 
# Get ProjDS from RetMapping: 
MapperSource = flowparser.clone("MapperSource")
if cfg.ProjMapper.MapperDS == "RetMap":
    wf.connect(retMapFlow,'outputspec.ROIvoxel',MapperSource,'in_ROIvoxel')
    wf.connect(retMapFlow,'outputspec.DSout',MapperSource,'in_DSout')
# Get ProjDS from Movie:     
elif cfg.ProjMapper.MapperDS == "Movie":
    wf.connect(movieFlow,'outputspec.ROIvoxel',MapperSource,'in_ROIvoxel')
    wf.connect(movieFlow,'outputspec.DSout',MapperSource,'in_DSout')
elif cfg.ProjMapper.MapperDS == "Other":
    MapperSource = flowparser.clone("MapperSource")
    wf.connect(otherFlow,'outputspec.ROIvoxel',MapperSource,'in_ROIvoxel')
    wf.connect(otherFlow,'outputspec.DSout',MapperSource,'in_DSout')

########################################################################
#### Whole Hyperalignment ####
#

### do core analysis ###
#
hypalAnalyFlow = create_postHypalAnalysis(cfgFile,name="hypalAnalyFlow") # load final analysis
wf.connect(info_nscans,'nscans' ,hypalAnalyFlow,'inputspec.nscans') 
wf.connect(ProjSource,('DSout_ROIvoxel_list',listsel,1),hypalAnalyFlow,'inputspec.ProjROIvoxel') 

### usage of Connectome ###
# ==================================================================== #
## with Connectome ##
if cfg.useConnectome is True:
    # ---------------------------------------------------------------- #
    ## define nodes for general Conenctome use ##
    #
    ## iterate across different sparse_radii ##
    info_sparseRad = pe.Node(name='info_sparseRad',interface=util.IdentityInterface(
        fields    = ['sparse_radius']),
        iterables = [('sparse_radius', cfg.sparse_radiusList)])
    
    # ---------------------------------------------------------------- #
    ## Without usage of PreConnectome, use initial Connectome to get
    # Mapper for Hypal    
    if cfg.usePreConnectome is False: 
        ## generate Connectome ##
        getConnectome = pe.Node(name='getConnectome',interface=Function(
            function     = make_Connectome,
            input_names  = ['sparse_radius','ROIvoxel',
                            'refDS','brainmask','refDSisMean',
                            'nproc','connectomeMethod'],
            output_names = ['ConnectomeOut','intersectionMaskPath']))
        getConnectome.inputs.brainmask        = cfg.brainmaskAll
        getConnectome.inputs.usePreConnectome = cfg.usePreConnectome
        getConnectome.inputs.connectomeMethod = cfg.connectomeMethod
        wf.connect(info_sparseRad,'sparse_radius',getConnectome,'sparse_radius')
        wf.connect(MapperSource,('DSout_ROIvoxel_list',listsel,1),getConnectome,'ROIvoxel')
        
        # to exclude ROIvoxel from DSRetMap #
        # to avoid overfitting of ROIvoxel to itself(?): #
        if cfg.excludeROIvoxel:
            getrefMapperDS = pe.Node(name='getrefMapperDS',interface=Function(
            function     = get_refMapperDS,
            input_names  = ['dsAll', 'nsubj'],
            output_names = ['DSout']))
            getrefMapperDS.inputs.nsubj   = cfg.nsubj
            getrefMapperDS.inputs.roitmpl = cfg.roitmpl
            wf.connect(MapperSource,('DSout_ROIvoxel_list',listsel,0),getrefMapperDS,'dsAll')
            wf.connect(getrefMapperDS,'DSout',getConnectome,'refDS')
        else:
            wf.connect(MapperSource,('DSout_ROIvoxel_list',listsel,0),getConnectome,'refDS')     
        wf.connect(getConnectome,'ConnectomeOut',hypalAnalyFlow,'inputspec.getMapperDSsave')    
    # ---------------------------------------------------------------- #
    ## Use Connectome from Hyperaligned Spheres ##
    if cfg.usePreConnectome is True:
        hyperhyperflow = create_doHyperHyper(cfgFile,name="hyperhyperflow")
        wf.connect(info_sparseRad,'sparse_radius',hyperhyperflow,'inputspec.sparse_radius')
        wf.connect(MapperSource,('DSout_ROIvoxel_list',listsel,0),hyperhyperflow,'inputspec.Mapper_refDS')    
        wf.connect(MapperSource,('DSout_ROIvoxel_list',listsel,1),hyperhyperflow,'inputspec.Mapper_ROIvoxel')    
        wf.connect(hyperhyperflow,'outputspec.Connectome_ProjDSvsHypalSphere',hypalAnalyFlow,'inputspec.getMapperDSsave')

# ==================================================================== #    
## without Connectome ##
elif cfg.useConnectome is False:   
    print "..not using Connectome"
    ## if no Connectome used just parse MapperDS to get Mapper for Hypal
    wf.connect(MapperSource,('DSout_ROIvoxel_list',listsel,1),hypalAnalyFlow,'inputspec.getMapperDSsave')
   
########################################################################
### save everything important for further analysis ###
#
                    
## datasink ##
datasink = pe.Node(interface=nio.DataSink(
    parameterization=False), name="datasink",overwrite=True)
datasink.inputs.base_directory = cfg.datasinkDir
if cfg.useConnectome:
    def getContainerName(innumber): return "srad"+str(innumber)
    wf.connect(info_sparseRad,('sparse_radius',getContainerName),datasink,'container')
    # PreConn used, save Spheres projected into thier subject specific Commonspace
    if cfg.usePreConnectome:
        wf.connect(hyperhyperflow,'outputspec.sphereCommonspace',datasink,'PreConn_sphereCommonspace')
        wf.connect(hyperhyperflow,'outputspec.Connectome_ProjDSvsHypalSphere',datasink,'HyperConnectome')
    else:
        wf.connect(getConnectome,'ConnectomeOut',datasink,'Connectome')

else: 
    datasink.inputs.container = "noConn"

wf.connect(hypalAnalyFlow,'outputspec.meanCS_Backproj2toSubjOut',datasink,'meanCS_BackprojectedToSubj')
wf.connect(hypalAnalyFlow,'outputspec.meanCS_CScorrToSubjOut',datasink,'meanCS_CorrelationInCommonspace')
wf.connect(hypalAnalyFlow,'outputspec.meanCS_CompToSubjOut',datasink,'meanCS_CompToSubjOut')
wf.connect(hypalAnalyFlow,'outputspec.meanCS_CompToSubjOut_data',datasink,'meanCS_CompToSubjOut_data')
wf.connect(hypalAnalyFlow,'outputspec.CompToSubjOut',datasink,'CompToSubjOut')
wf.connect(hypalAnalyFlow,'outputspec.CompToSubjOut_data',datasink,'CompToSubjOut_data')
wf.connect(hypalAnalyFlow,'outputspec.Backproj2toSubjOut',datasink,'BackprojectedToSubj')
wf.connect(hypalAnalyFlow,'outputspec.CScorrToSubjOut',datasink,'CorrelationInCommonspace')
wf.connect(hypalAnalyFlow,'outputspec.AnatAlign2HyperOut',datasink,'AnatAlignVSHypAl')
wf.connect(hypalAnalyFlow,'outputspec.Mapper',datasink,'Mapper')
wf.connect(hypalAnalyFlow,'outputspec.CSsubj',datasink,'CSsubj')

# save overlays from retmapFlow if they are generated .. #
if (cfg.ProjMapper.MapperDS == "RetMap") or (cfg.ProjMapper.ProjDS == "RetMap"):
    ## datasink just for retmap output ##
    datasinkRM = pe.Node(interface=nio.DataSink(
        parameterization=False), name="datasinkRM",overwrite=True)
    datasinkRM.inputs.base_directory = cfg.datasinkDir
    datasinkRMall = pe.Node(interface=nio.DataSink(
        parameterization=False), name="datasinkRMall",overwrite=True)
    datasinkRMall.inputs.base_directory = cfg.datasinkDir
    
    def getContainerName(innumber): return 'sub'+innumber

    wf.connect(info_subjects,('subj',getContainerName),datasinkRM,'container') 
    wf.connect(retMapFlow,'outputspec.phaseFiles',datasinkRM,'phaseFiles')
    wf.connect(retMapFlow,'outputspec.surfacefiles_R',datasinkRM,'surfR')
    wf.connect(retMapFlow,'outputspec.surfacefiles_L',datasinkRM,'surfL')
    wf.connect(retMapFlow,'outputspec.coreg2anat',datasinkRM,'coreg2anat')
    wf.connect(retMapFlow,'outputspec.DSout',datasinkRMall,'retmapFullDSList')
    wf.connect(retMapFlow,'outputspec.ROIvoxel',datasinkRMall,'retmapROIvoxelList')

# save the retmap overlays from the projected data #
if cfg.ProjMapper.ProjDS == "RetMap":
    # get overlays
    getOutniftiOverlays = create_makeOutniftilistOverlays(cfgFile,name="getOutniftiOverlays") 
    getOutniftiOverlays_meanCS = create_makeOutniftilistOverlays(cfgFile,name="getOutniftiOverlays_meanCS")
    wf.connect(info_subjects,'subj',getOutniftiOverlays,'inputspec.subj')
    wf.connect(info_subjects,'subj',getOutniftiOverlays_meanCS,'inputspec.subj')
    wf.connect(hypalAnalyFlow,'outputspec.outnifti_raw',getOutniftiOverlays,'inputspec.outnifti')
    wf.connect(hypalAnalyFlow,'outputspec.meanCS_outnifti_raw',getOutniftiOverlays_meanCS,'inputspec.outnifti')
    
    wf.connect(getOutniftiOverlays,'outputspec.surfacefiles_R',datasink,'backproj_surfR')
    wf.connect(getOutniftiOverlays,'outputspec.surfacefiles_L',datasink,'backproj_surfL')
    wf.connect(getOutniftiOverlays,'outputspec.coreg2anat',datasink,'backproj_coreg2anat')
    wf.connect(getOutniftiOverlays_meanCS,'outputspec.surfacefiles_R',datasink,'meanCS_backproj_surfR')
    wf.connect(getOutniftiOverlays_meanCS,'outputspec.surfacefiles_L',datasink,'meanCS_backproj_surfL')
    wf.connect(getOutniftiOverlays_meanCS,'outputspec.coreg2anat',datasink,'meanCS_backproj_coreg2anat')

if (cfg.ProjMapper.MapperDS == "Movie") or (cfg.ProjMapper.ProjDS == "Movie"):
    datasinkMV = pe.Node(interface=nio.DataSink(
        parameterization=False), name="datasinkMV",overwrite=True)
    datasinkMV.inputs.base_directory = cfg.datasinkDir
    wf.connect(movieFlow,'outputspec.ROIvoxel',datasinkMV,'movieROIvoxelList')
    wf.connect(movieFlow,'outputspec.DSout',datasinkMV,'movieFullDSList')
if (cfg.ProjMapper.MapperDS == "Other") or (cfg.ProjMapper.ProjDS == "Other"):
    datasinkMV = pe.Node(interface=nio.DataSink(
        parameterization=False), name="datasinkMV",overwrite=True)
    datasinkMV.inputs.base_directory = cfg.datasinkDir
    wf.connect(otherFlow,'outputspec.ROIvoxel',datasinkMV,'movieROIvoxelList')
    wf.connect(otherFlow,'outputspec.DSout',datasinkMV,'movieFullDSList')
#######################################################################
### run workflow ###
#
if __name__ == '__main__':
    wf.run(plugin = 'CondorDAGMan')
    #wf.run() #run local (on medusa)
    wf.write_graph(graph2use='flat') 
