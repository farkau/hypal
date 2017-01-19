# FK2015: Meta pipeline to get anatomy based reference for the Hyperalignment

#############################################
## See LICENSE for licensing and copyright ##
#############################################

# ==================================================================== #
from HA_configLoader_v2 import configLoader
from HA_extFunctions import (
        getAnatCdist
        )
from HA_extWFlows import (
    create_rmFlow,
    create_makeOutniftilistOverlays,
    )

from mvpa2.base import debug
import os 
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio           
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
wf          = pe.Workflow(name= "wf_AnatCdist")
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
        iterables = [('nscans', ([cfg.MapperDS.nscans]))])
info_subjects = pe.Node(name='info_subjectsMain',interface=util.IdentityInterface(
        fields    = ['subj']),
        iterables = [('subj', cfg.subjlist)])
########################################################################
#### Get Data #####
#
print ".. using RetMapFlow"
retMapFlow = create_rmFlow(cfgFile,name="retMapFlow",useGrpTmpl=True)
wf.connect(info_subjects,'subj',retMapFlow,'inputspec.subj')

########################################################################
### Get meanAnat(other) ###
#
getCdistAnat = pe.Node(name='getCdistAnat',interface=Function(
    function     = getAnatCdist,
    input_names  = ['DSlist','makenifti'],
    output_names = ['AnatCdistSummary','DSmeanOther','outnifti',
                    'outnifti_raw','CompToSubjOut_data']),)
getCdistAnat.inputs.makenifti = cfg.makenifti
wf.connect(retMapFlow,'outputspec.ROIvoxel',getCdistAnat,'DSlist')

########################################################################
### save everything important for further analysis ###
#
                    
## datasink ##
datasink = pe.Node(interface=nio.DataSink(
    parameterization=False), name="datasink",overwrite=True)
datasink.inputs.base_directory = cfg.datasinkDir
datasink.inputs.container = "noConn"

# save overlays from retmapFlow if they are generated .. #
if (cfg.ProjMapper.MapperDS == "RetMap") or (cfg.ProjMapper.ProjDS == "RetMap"):
    ## datasink just for retmap output ##
    datasinkRM = pe.Node(interface=nio.DataSink(
        parameterization=False), name="datasinkRM",overwrite=True)
    datasinkRM.inputs.base_directory = cfg.datasinkDir
    def getContainerName(innumber): return 'sub'+innumber

    wf.connect(info_subjects,('subj',getContainerName),datasinkRM,'container') 
    
    wf.connect(retMapFlow,'outputspec.phaseFiles',datasinkRM,'phaseFiles')
    wf.connect(retMapFlow,'outputspec.surfacefiles_R',datasinkRM,'surfR')
    wf.connect(retMapFlow,'outputspec.surfacefiles_L',datasinkRM,'surfL')
    wf.connect(retMapFlow,'outputspec.coreg2anat',datasinkRM,'coreg2anat')
    wf.connect(retMapFlow,'outputspec.DSout',datasink,'retmapFullDSList')
    wf.connect(retMapFlow,'outputspec.ROIvoxel',datasink,'retmapROIvoxelList')

    wf.connect(getCdistAnat,'AnatCdistSummary',datasink,'AnatCdistSummary')
    wf.connect(getCdistAnat,'DSmeanOther',datasink,'DSmeanOther')
    wf.connect(getCdistAnat,'CompToSubjOut_data',datasink,'CompToSubjOut_data')
    wf.connect(getCdistAnat,'outnifti_raw',datasink,'outnifti')

# save the retmap overlays from the projected data #
if cfg.ProjMapper.ProjDS == "RetMap":
    # get overlays
    getOutniftiOverlays = create_makeOutniftilistOverlays(cfgFile,name="getOutniftiOverlays",useGrpTmpl=True) 
    wf.connect(info_subjects,'subj',getOutniftiOverlays,'inputspec.subj')
    wf.connect(getCdistAnat,'outnifti_raw',getOutniftiOverlays,'inputspec.outnifti')
    wf.connect(getOutniftiOverlays,'outputspec.surfacefiles_R',datasink,'backproj_surfR')
    wf.connect(getOutniftiOverlays,'outputspec.surfacefiles_L',datasink,'backproj_surfL')
    wf.connect(getOutniftiOverlays,'outputspec.coreg2anat',datasink,'backproj_coreg2anat')


#######################################################################
##### run workflow ###
#
if __name__ == '__main__':
    wf.run(plugin = 'CondorDAGMan')
    #wf.run() #run local
    wf.write_graph(graph2use='flat') 
