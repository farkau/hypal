# FK2015: all the stuff to check how/if the hypal worked out:
# - do the plots!
# - write the summary for the import into Latex

#############################################
## See LICENSE for licensing and copyright ##
#############################################

### import stuff to do magic ###
#
import os
import commands
import matplotlib.pyplot as plt
from mvpa2.suite import h5load
import os
import numpy as np
from numpy import unique, size, sort
from HA_configLoader_v2 import configLoader
from HA_extFunctions import _findfile as findfile
from HA_extFunctions import _getCosCdist, _ds2deg, _circdescribe
import scipy.stats as sps
import copy
from matplotlib.mlab import find
import pandas as pd
import matplotlib
import scipy.spatial.distance as sd

# nicer plotting
import seaborn as sns
sns.set(style="whitegrid", color_codes=True)


# ==================================================================== #
## little helpers ##
#
def circCdist(a,b):
    """
    do circstat as used for the other cdists for hypal comparison
    """
    
    a_cos = np.cos(np.deg2rad(a))
    b_cos = np.cos(np.deg2rad(b))

    # correlate transformed cos values (nChannels for each ph value)
    # (float to avoid conversion error)
    tmp = sd.cdist(a_cos.reshape(1,len(a_cos)),b_cos.reshape(1,len(b_cos)),'co').squeeze()
    #tmp = sd.cdist(a_cos,b_cos,'correlation').squeeze()
    return float(tmp)

def cdistVal(x,y):
    import scipy.spatial.distance as sd
    return float(sd.cdist(x.reshape(1,len(x)),y.reshape(1,len(y)),'co').squeeze())

def fileloader(searchdir, tmpl, startline=1, cond=""):
    """
    gets all files in "searchdir" fitting template "tmpl"
    and loads them into a list via h5load. Ignore "ignoreLine" while
    loading 
    """
    
    from HA_extFunctions import _findfile as findfile
    from mvpa2.suite import h5load
    Files = findfile(searchdir, tmpl)
    outdata = []
    
    ## sort filenames according "human sorting" 
    # (1,2,10 instead of 1,10,2)
    #TODO: merge both versions into one
    def getFromSu_meanVals(name):
        return int(name.split('leavesubjnr')[-1].split('of')[0])
    
    # separate search mask (don't know how to pass inside sort .. )
    def getFromSu_csMeanData(name):
        return int(name.split('of')[-1].split('subj.h5')[0])
    
    if (size(Files) == 0): # if no files found 
        Files.sort()
    elif  cond=="csMeanData":
        # if data should be loaded, use different search pattern
        Files.sort(key=getFromSu_csMeanData) 
    elif ("leavesubjnr" in Files[0]) :
        Files.sort(key=getFromSu_meanVals) 
    else: 
        # .. but only if needed
        Files.sort()

    out0, out1 = [],[]
    for filename in Files:
        #print 'DEBUG: .. loading: ', filename
        tmp = h5load(filename)
        if cond=="anatVShypal":
            for tt in tmp[startline:]:
                out0.append([tt[0][0],tt[3],tt[4]])
                out1.append([tt[0][1],tt[3],tt[4]])
        elif cond=="meanCS":
            out0.append([tmp[1][0][0],tmp[1][1]])
            out1.append([tmp[1][0][1],tmp[1][1]])
        elif cond=="csMeanData":
            # load the data according toSubj for later scatterplot 
            # of dsback VS dsorg
            out0.append(tmp[1]) 
        elif cond=="anat":
            # for anat comparison: [dsMean_other_masked, ds_masked,toSubj]
            out0 = tmp[startline:] # is still in (cosChannels, nVoxel)
            out1 = [] #placeholder .. wahtfor??
        else:                
            for tt in tmp[startline:]:
                out0.append([tt[0][0],tt[1],tt[2]])
                out1.append([tt[0][1],tt[1],tt[2]])
          
    return out0, out1
# -------------------------------------------------------------------- #
def plot_cdist(inmat,nsubj,titleadd, sparse_radius):
    """
    plot backprojected data for all subjects with all  combinations
    """
    
    import numpy as np
    import matplotlib.pyplot as plt

    # prpeare plotting
    subjrange = range(nsubj)
    min_val, max_val, diff = 0., nsubj , 1.
    N_points = (max_val - min_val) / diff

    fig, ax = plt.subplots()
    imshow_data = np.around(inmat, decimals = 2)
    cax = ax.imshow(imshow_data, interpolation='nearest',vmin = 0,vmax = 1)
    cbar = fig.colorbar(cax, ticks=[0,0.5,1])

    #text portion
    ind_array = np.arange(min_val, max_val, diff)
    x, y = np.meshgrid(ind_array, ind_array)
    for x_val, y_val, inval in zip(x.flatten(), y.flatten(),imshow_data.flatten()):
        c = inval
        ax.text(x_val, y_val, c, va='center', ha='center')

    #set tick marks for grid
    ax.set_title("cdist "+titleadd+" for sparse_radius %i for projections" % sparse_radius)
    ax.set_xlabel("from Subject")
    ax.set_ylabel("to Subject")
    ax.set_xticks(np.arange(min_val-diff/2, max_val-diff/2))
    ax.set_yticks(np.arange(min_val-diff/2, max_val-diff/2))
    ax.set_xticklabels(subjrange,ha = "left")
    ax.set_yticklabels(subjrange,va = "bottom")
    ax.set_xlim(min_val-diff/2, max_val-diff/2)
    ax.set_ylim(min_val-diff/2, max_val-diff/2)
    ax.grid()
    plt.show()
    return
# -------------------------------------------------------------------- #
def plot_cdistNegative(inmat,nsubj,titleadd, sparse_radius):
    """
    plot backprojected data for all subjects with all  combinations
    """
    
    import numpy as np
    import matplotlib.pyplot as plt

    # prpeare plotting
    subjrange = range(nsubj)
    min_val, max_val, diff = 0., nsubj , 1.
    N_points = (max_val - min_val) / diff

    fig, ax = plt.subplots()
    imshow_data = np.around(inmat, decimals = 2)
    cax = ax.imshow(imshow_data, interpolation='nearest',vmin = -1,vmax = 1)
    cbar = fig.colorbar(cax, ticks=[-1,0,1])

    #text portion
    ind_array = np.arange(min_val, max_val, diff)
    x, y = np.meshgrid(ind_array, ind_array)
    for x_val, y_val, inval in zip(x.flatten(), y.flatten(),imshow_data.flatten()):
        c = inval
        ax.text(x_val, y_val, c, va='center', ha='center')

    #set tick marks for grid
    ax.set_title("cdist "+titleadd+" for sparse_radius %i for projections" % sparse_radius)
    ax.set_xlabel("from Subject")
    ax.set_ylabel("to Subject")
    ax.set_xticks(np.arange(min_val-diff/2, max_val-diff/2))
    ax.set_yticks(np.arange(min_val-diff/2, max_val-diff/2))
    ax.set_xticklabels(subjrange,ha = "left")
    ax.set_yticklabels(subjrange,va = "bottom")
    ax.set_xlim(min_val-diff/2, max_val-diff/2)
    ax.set_ylim(min_val-diff/2, max_val-diff/2)
    ax.grid()
    plt.show()
    return
# -------------------------------------------------------------------- #
def meanExclDiag(A):
    """
    exclude the values from the diagonal of matrix A for calculating the mean.
    
    This way ignore the backprojection to self. Just take care about
    correlation of projections into other subjects
    """
    
    import numpy as np
    sumExclDiag = np.sum((np.triu(A),np.tril(A))) - 2*(np.diag(A).sum())
    nPointsExclDiag = np.array(A.shape).prod() - A.shape[0]
    return sumExclDiag / nPointsExclDiag

def getCScorrTRsr(csCorrTR,cfg,sparse_radius, scan): 
    """ get matrix for plotting """
    
    csCorrTRsr = np.round(csCorrTR[:,0].astype(np.float64).reshape(cfg.nsubj,cfg.nsubj),2)
    csCorrTRsrMean = np.round(meanExclDiag(csCorrTRsr),2)
    csCorrTRAllsr = [csCorrTRsr,sparse_radius,csCorrTRsrMean,scan]
    return csCorrTRAllsr

def plotCScorr(CSCorr,condition,cfg, topResults, doPlots):
    """ plot only best result """

    # sum up all correlations of refDSvoxel with ROIvoxel and order them by
    # thier sum of cdist inside dsm. The highest cdist across all subj
    #  is rank 1, the best contributing voxel is last rank.
    ranks = []
    ranks = sps.rankdata(CSCorr[0][0])
    ranksSorted = copy.deepcopy(ranks)
    ranksSorted.sort()

    for rank, csCorr in zip(ranks, CSCorr):
        if rank <= (ranksSorted[:topResults].max()): #take only best "topResults"
            toPlot = csCorr[0] 
            sparse_radius = csCorr[1]
            csmean = csCorr[2]
            plotkey = condition+" tSeries in CS of each singleSubject (NOT with meanDS !)" + \
                        "\n useConn"+ str(cfg.useConnectome) + \
                        " preConn"+ str(cfg.usePreConnectome) + \
                        " excludeROI"+ str(cfg.excludeROIvoxel) + \
                        " MapperDS: "+str(cfg.ProjMapper.MapperDS) +\
                         " nscanMapper" + str(csCorr[3])
            
            if doPlots: plot_cdist(toPlot,
                        cfg.nsubj,
                        plotkey,
                        sparse_radius
                        )
            print plotkey, "sparse_radius: ", sparse_radius
            print ">> mean, rank",csmean, rank 
    return sparse_radius

def loadCSmean(dirname,filename,outPath,cfg,cond=""):
    import os
    import numpy as np
    # check which dataset is projected and adopt nscans accordingly
    nscans = cfg.MapperDS.nscans
    if cfg.useConnectome == True: sr_list = cfg.sparse_radiusList
    else: sr_list = [-1]
    
    ### load all files  or each sparse_radius and toSubj ###
    ## for Correlation in Commonspace ##
    #
    CSCorrAllsr_ECC, CSCorrAllsr_POL = [],[]
    for scan in range(1,nscans+1):
        for sparse_radius in sr_list:
            if cfg.useConnectome == True: srad = "srad"+str(sparse_radius)
            else: srad = "noConn"
            searchPath = os.path.join(outPath,srad,dirname) 
            CSCorr1Tmpl = filename+"_"+str(scan)+'scans_leavesubjnr*of*.h5'
            #print "DEBUG: .. loading: ", os.path.join(outPath,CSCorr1Tmpl)
            # only load last line for each, because fromSubjCS is same for all
            # (because meanCS), only toSubj varies
            csCorr_ECC, csCorr_POL = fileloader(searchPath, CSCorr1Tmpl,cfg.nsubj,cond)            
            csCorr_ECC = np.array(csCorr_ECC)
            csCorr_POL = np.array(csCorr_POL)

            if len(csCorr_ECC) > 0:
                CSCorrAllsr_ECC.append([csCorr_ECC,sparse_radius,np.round(np.mean(csCorr_ECC[:,0]),2)])
                CSCorrAllsr_POL.append([csCorr_POL,sparse_radius,np.round(np.mean(csCorr_POL[:,0]),2)])
    return CSCorrAllsr_ECC, CSCorrAllsr_POL  

def getMeanANDMinCScorr(CScorr):
    """
    get mean cdist, then get mean for each sparse_radius and then
    collect cdist for all subjects for this best sparse_radius
    """
    
    import numpy as np
        
    CScorr_bestSR, CScorr_mean = [],[]

    # get mean cdist
    [CScorr_mean.append(val[2]) for val in CScorr]
    CScorr_mean = np.array(CScorr_mean)

    # get srad of best (min) cdist ..
    minIdx = find(CScorr_mean == CScorr_mean.min())[0]
    srad_best = CScorr[minIdx][1]

    return CScorr[minIdx],srad_best

def makePhaseCorrPlots(cfg, outPath, basedir, anatData=False, saveImg=False, doSingleSuPlots=False, picsavedir=""):
    from scipy.stats import linregress # for linear regression
    from pylab import polyfit, polyval
    
    # when I save the plots anyway and want to keep the ratios, so don't show the plots ..
    if saveImg:
        plt.ioff()
    
    SHAPE = []
    ## do scatter plots ##
    #
    if cfg.useConnectome == True: sr_list = cfg.sparse_radiusList
    else: sr_list = [-1]
    
    if anatData: 
        dirname = "CompToSubjOut_data"
        dataname = 'CorrToSubj2fromSubj_data_nSubj*.h5'
        cond = "anat"
        startLine = 1  
    else: 
        dirname = "meanCS_CompToSubjOut_data"    
        dataname = 'CorrToSubj2fromSubj_data_*scans_leavesubjnr*.h5'
        cond = "csMeanData"
        startLine = 0
     
     
    DF = pd.DataFrame()
    for sparse_radius in sr_list:
        if cfg.useConnectome == True: srad = "srad"+str(sparse_radius)
        else: srad = "noConn"
    
        # load CsmeanData
        searchPath = os.path.join(outPath,srad,dirname) 

        [Data, __] = fileloader(searchPath, dataname, startLine, cond)

        ## plot scatter for all meanCS (each one for each subj excluded)
        CompToSubj,Circdesc = [],[]
        for toSubj in xrange(len(Data)):    
            
            dsback = Data[toSubj][0][0]
            dsorg  = Data[toSubj][0][1]
            
            postHypal_ph = _ds2deg(dsback)
            refds_ph = _ds2deg(dsorg)

            ## voxelwise cdist for all subj: 
            # compare "pre vs pre" and "pre vs post" to get the qualtiy of
            # alignment of anatomical and hyperaligment
            channels = refds_ph.a.deg2cos_channels
            hypalcdist = _getCosCdist(dsorg, dsback, channels)
            ## only use circdesc of 2nd ds (fromSuDsHypal)
            #circdesc = [_circdescribe(postHypal_ph.S[0]), _circdescribe(postHypal_ph.S[1])]
            
            xP,yP = refds_ph[0,:].S.T, postHypal_ph[0,:].S.T    
            xE,yE = refds_ph[1,:].S.T, postHypal_ph[1,:].S.T    
            ### get linear regeression 
            ## from [ http://www.wired.com/2011/01/linear-regression-with-pylab/ ]
            #(slP,intP) = polyfit(xP.squeeze(),yP.squeeze(),1)
            #yLinRegP = polyval([slP,intP],xP) 
            #(slE,intE) = polyfit(xE.squeeze(),yE.squeeze(),1)
            #yLinRegE = polyval([slE,intE],xE) 

            # transform to pandas dataset

            d = zip( xE.squeeze().T, yE.squeeze().T, xP.squeeze().T, yP.squeeze().T, [toSubj]*xE.shape[0], [sparse_radius]*xE.shape[0] )
            df = pd.DataFrame(d,columns=['x_ecc','y_ecc','x_pol','y_pol', 'toSubj','sparse_radius'])
            # now attach the df to generate the allInOneDataFrame .. yeah!!
            DF = pd.concat([DF, df])

    # reset index. Otherwise it will just repeat the index for every sub-df
    DF.reset_index(drop=True)

    # make nice plots
    ticks = range(0,361,60)
    sns.set_palette('dark')
    sns.set(style="darkgrid", color_codes=True, font_scale=2)
    
    # bigger legend
    matplotlib.rc("legend", 
                    fontsize=20, shadow=True, frameon=True,
                    numpoints=3) 
    
    # all subj plots in one plot
    statName = "correlation distance"
    with sns.color_palette("Blues_d"):
        title = ("Polar angle acorss all Subjects")

        gP= sns.jointplot(x="x_pol", y="y_pol", data=DF, kind="reg",
                n_boot=10000, ci=99,
                xlim=(0, 360), ylim=(0, 360),
                color="b", size=7,
                line_kws={'color': 'black'},
                stat_func=None
                )
        # add my circStat instead
        gP.annotate(circCdist, template="{stat}: {val:.2f}",
                    stat=statName, loc="upper left", frameon=True, fontsize=20)
        gP.set_axis_labels('target phase values [deg]','transfered phase values [deg]')
        gP.ax_marg_y.set_yticks(ticks)
        gP.ax_marg_x.set_xticks(ticks)
        if saveImg: gP.savefig(os.path.join(picsavedir + '/polPhaseCorr.svg'),bbox_inches='tight')

        if doSingleSuPlots:
            print ".. doing linear regression: 1/2"
            pP = sns.lmplot(x="x_pol", y="y_pol", data=DF, \
                    ci=95, robust=False, 
                    col="toSubj", 
                    col_wrap=2, size=3,
                    palette="muted").set(
                    xlim=(0, 360), 
                    ylim=(0, 360),
                    xticks=ticks, 
                    yticks=ticks,
                    xlabel="target pol",
                    ylabel="transfered pol"
                    )
            pP.fig.suptitle(title)
            if saveImg: pP.savefig(os.path.join(picsavedir + '/polAllSubj.png'),bbox_inches='tight')
            
    with sns.color_palette("Reds_d",desat=1):
        title = ("Eccentricity acorss all Subjects")
        # non-circular statistics for ECC
        cdistVal = lambda x,y : float(sd.cdist(x.reshape(1,len(x)),y.reshape(1,len(y)),'correlation').squeeze())
        
        gE = sns.jointplot(x="x_ecc", y="y_ecc", data=DF, kind="reg", 
                n_boot=10000, ci=99, 
                xlim=(0, 360), ylim=(0, 360), 
                color="r", size=7,
                #line_kws={'color': sns.xkcd_rgb['hot magenta']}
                line_kws={'color': 'black'},
                stat_func=None
                )
        # add my circStat instead
        gE.annotate(cdistVal, template="{stat}: {val:.2f}",
                    stat=statName, loc="upper left", frameon=True, fontsize=20)
        gE.set_axis_labels('target phase values [deg]','transfered phase values [deg]')
        gE.ax_marg_y.set_yticks(ticks)
        gE.ax_marg_x.set_xticks(ticks)
        if saveImg: gE.savefig(os.path.join(picsavedir + '/eccPhaseCorr.svg'),bbox_inches='tight')
        
        if doSingleSuPlots:
            print ".. doing linear regression: 2/2"
            pE = sns.lmplot(x="x_ecc", y="y_ecc", data=DF, \
                    ci=95, robust=False, 
                    col="toSubj", 
                    col_wrap=2, size=3,
                    palette="muted").set(
                    xlim=(0, 360), 
                    ylim=(0, 360),
                    xticks=ticks, 
                    yticks=ticks,
                    xlabel="target ecc",
                    ylabel="transfered ecc"
                    )
            pE.fig.suptitle(title)
            if saveImg: pE.savefig(os.path.join(picsavedir + '/eccAllSubj.png'),bbox_inches='tight')

    # close all figures after saving
    if saveImg: 
        print ".. saved png to: ", picsavedir
        plt.close("all")

    return

# ==================================================================== #

def plotdata(basedir, doPlots=False,doPhasePlots=False,doSingleSuPlots=False, saveImg=False):
    """ 
    plot data
    
    basedir should be the outdir from the hypal processing sth like 
    "wf_out_ROIV1_JuelThres99_2tmpl_sRadList87_exclRoiTrue_NOpreConn_nsubj5_allSubj2CS_conditionorTestingConn"
    """
    
    #TODO:
    # - take mean cdist to plot only top3 at the end
    print "####### "+basedir+" #########"
    # ==================================================================== #
    ## load cfg file ##
    hypalDir = "/home/fkaule/MRI/projects/Hyperalignment/"
    cfgFilePath = os.path.join(hypalDir + basedir)
    print ".. searchPath:  %s" % cfgFilePath
    cfgFile = findfile(cfgFilePath,'cfgBackup.cfg')[0]
    outPath = os.path.dirname(cfgFile)
    print "..loading config from ", cfgFile
    cfg = configLoader(cfgFile)
    #print "DEBUG: ..done"

    topResults = 3 # plot only the top X results
    print " useConn"+str(cfg.useConnectome)+" preConn"+str(cfg.usePreConnectome)+" exclROIvoxel"+str(cfg.excludeROIvoxel)
    print 'smaller is better cause cdist-co is used'
    subjrange = range(cfg.nsubj)    
    # check which dataset is projected and adopt nscans accordingly
    nscans = cfg.MapperDS.nscans
    if cfg.useConnectome == True: sr_list = cfg.sparse_radiusList
    else: sr_list = [-1]

    ### load all files  or each sparse_radius and toSubj ###
    ### for Correlation in Commonspace ##
    #
    CSCorrTRAllsr_ECC, CSCorrTRAllsr_POL = [],[]
    for scan in range(1,nscans+1):
        for sparse_radius in sr_list:
            if cfg.useConnectome == True: srad = "srad"+str(sparse_radius)
            else: srad = "noConn"
            searchPath = os.path.join(outPath,srad,"CorrelationInCommonspace") 
            CSCorrTR1Tmpl = "CScorrToSubj_"+str(scan)+'scans_leavesubjnr*of*.h5'
            #print "DEBUG: .. loading: ", os.path.join(outPath,CSCorrTR1Tmpl)
            csCorrTR_ECC, csCorrTR_POL = fileloader(searchPath, CSCorrTR1Tmpl,1)            
            csCorrTR_ECC = np.array(csCorrTR_ECC)
            csCorrTR_POL = np.array(csCorrTR_POL)
            
            if len(csCorrTR_ECC) > 0:
                CSCorrTRAllsr_ECC.append(getCScorrTRsr(csCorrTR_ECC,cfg,sparse_radius,scan))
                CSCorrTRAllsr_POL.append(getCScorrTRsr(csCorrTR_POL,cfg,sparse_radius,scan))
                
    # do the plots of correlation in Commonspace #
    # .. just do all the stuff before to get sradBest .. nothing else if you don't need the plots (thus, if doPlots == False)
    sradBest = plotCScorr(CSCorrTRAllsr_ECC,"ECC",cfg,topResults,doPlots) # sradBest same for POL and ECC (.. never got different)
    plotCScorr(CSCorrTRAllsr_POL,"POL",cfg,topResults,doPlots)
        
    """
    ================================================================ 
    get difference between the anatomical and hyperaligned datasets, to see the improvement:
    POSITIVE values Improvement by HYPAL: cdist(org vs Anataligned) > cdist(org vs Hypaligned) 
    NEGATIVE values bad HYPAL: cdist(org vs Anataligned) < cdist(org vs Hypaligned)
    """
    if cfg.useConnectome == True: srad = "srad"+str(sradBest)
    else: srad = "noConn"
    
    ## do that for all different nScans ##
    SScorrAll_ECC, SScorrAll_POL = loadCSmean('meanCS_CompToSubjOut','CorrToSubj2fromSubj',outPath,cfg,"meanCS")

    scansPerNScans = len(sr_list)
    nscans = len(SScorrAll_ECC)/scansPerNScans
    for scan in range(0,nscans*2,2):

        ## Cdist(backprojected mean, orgSubj)@subjectSpace ##    
        ECCcorr_SS, ECCbestSRad_SS = getMeanANDMinCScorr(SScorrAll_ECC[scan:scan+2])
        POLcorr_SS, POLbestSRad_SS = getMeanANDMinCScorr(SScorrAll_POL[scan:scan+2])

        ## give overview ##
        print "\n ### nScans: ", str(scan/2 +1),"(max nScans if just one nScans tested) ##"
        print "\n Summary of correlation with meanDS: \nCdist@SS >> mean ECC; POL: ", ECCcorr_SS[2], POLcorr_SS[2],"\n"
        print ">> Cdist(backprojected mean, orgSubj)@SubjSpace"
        print "ECCcorr: \n", ECCcorr_SS
        print "\n POLcorr: \n", POLcorr_SS

    # ================================================================ #    
    
    if doPhasePlots: makePhaseCorrPlots(cfg, outPath, basedir, False, saveImg, doSingleSuPlots, cfgFilePath)  
    
    # return meanECC, meanPOL, MapperDScond, useConn, preConn, roitmpl, [ECCcorr_all4sd, POLcorr_all4sd]
    return [np.round(ECCcorr_SS[2],2), np.round(POLcorr_SS[2],2), cfg.ProjMapper.MapperDS, "useConn"+str(cfg.useConnectome), "preConn"+str(cfg.usePreConnectome), cfg.roitmpl, [ECCcorr_SS[0][:,0], POLcorr_SS[0][:,0]]]

########################################################################
def anatSummary(wfDir, doPhasePlots=False, doSingleSuPlots=False, saveImg=False, picsavedir=""):  
    """
    get anatSummary
    """
    
    hypaldir = '/home/fkaule/MRI/projects/Hyperalignment/'
    basedir = hypaldir + wfDir + '/*/noConn/'
    anatSummary = commands.getoutput('ls '+basedir+'AnatCdistSummary/AnatCdistSummary_nSubj8.h5')

    # do the save plots as for the hypal Phase Correlations
    ## load cfg file ##
    cfgFilePath = os.path.join(hypaldir + wfDir)
    cfgFile = findfile(cfgFilePath,'cfgBackup.cfg')[0]
    outPath = os.path.dirname(cfgFile)
    print "..loading config from ", cfgFile
    cfg = configLoader(cfgFile)
    if doPhasePlots: makePhaseCorrPlots(cfg, outPath, wfDir, True, saveImg, doSingleSuPlots, cfgFilePath)
    tmpACd = h5load(anatSummary)
    tmpACd = tmpACd[1:]

    print tmpACd

    # extract POL and ECC anatCdist(ds,MeanOther)
    eccACD = [s[0][0] for s in tmpACd]
    polACD = [s[0][1] for s in tmpACd]

    eccACD_avg = np.mean(eccACD)
    polACD_avg = np.mean(polACD)

    print "mean anatCdist(ds,DSmeanOther) \n ECCcdist: ",eccACD_avg,"\n POLcdist: ",polACD_avg
    return [eccACD_avg, polACD_avg, cfg.roitmpl, [eccACD, polACD]]

########################################################################
def plotSummary(DSSS, DSSS_anat, plotType="box"):
    """
    summarize the plots
    
    plotType can be: box, violin, bar
    """
    
    from numpy import mean as mean
    import numpy as np
    plt.ion()
    ### use pandas and seaborn for nice plots ###
    #
    # convert to pandas
    DF = convertDSSS2pandasDFrame(DSSS, DSSS_anat)
    # Draw a nested barplot to show the cdist for the different 
    # methods (conn, noConn..) across for the different mapper datasets (Retmap, movie ..)
    # also nice kind can be: box, violin, bar
    DFplot = DF
    
    sns.set_context("talk")
    g = sns.factorplot(x="UsedDataset", y="cdist", hue="method", data=DFplot, 
                    kind=plotType,
                    palette="muted",
                    col="ecc_pol", 
                    n_boot=10000, ci=99,
                    )
    # Use semantically meaningful titles for the columns
    titles = ["eccentricity", "polar angle"]
    sns.despine(left=True, bottom=True)
    g.despine(left=True)
    g.set_ylabels("correlation distance of phase values")
    g.set_xlabels("training dataset")
    g.set_xticklabels(["Retinotopy","Audio-visual\nmovie","Audio-movie","Anatomical"], rotation=45)    
    for ax, title in zip(g.axes.flat, titles):
        # Set a different title for each axes
        ax.set(title=title)
        # Make the grid horizontal instead of vertical
        ax.xaxis.grid(False)
        ax.yaxis.grid(True)

    return

########################################################################
def getLatexVars_DSSS(DSSS, outsavepath):
    """
    get latex readable variables

    take the DSSS list and convert it into a variable list for tex
    ie. \newcommand{\OtheruseConnTruepreConnFalseleftGMdil3}{0.92}
    FRK2015
    """
    
    # TODO: 
    # - fix problem: multiple datapints from exactly same condition will get same varname 
    #     -> bad for tex!!
    
    # save to a texfile
    outfile = open(outsavepath,"w")

    DF = convertDSSS2pandasDFrame(DSSS, [])
    
    ## write the mean acorss both hemispheres ##
    # iterater across all the combinations and calculate mean #
    # (this way have the option to calc SEM or so if needed
    for Mapper in pd.unique(DF[DF.keys()[0]]):
        for method in pd.unique(DF[DF.keys()[1]]):
            for ecc_pol in pd.unique(DF[DF.keys()[4]]):
                # could give sth like: method, ecc_pol, Mapper = "Conn", "ecc", "RetMap"
                selector = 'method == "%s" & ecc_pol == "%s" & UsedDataset == "%s"' % (method, ecc_pol, Mapper)
                meanval = DF.query(selector).mean()
                
                if method == "normal Hyperalignment": m="normalHypal"
                else: m = method
                # will give sth like 'RetMapuseConnTruepreConnFalserightGMdil'
                varname = '%s%s%sMean' % (ecc_pol, Mapper, m)
                
                texcmd = '\\newcommand{\%s}{%.2f}\n' % (varname, meanval)
                outfile.write(texcmd) ## write to file
    
    ## .. and here the single values foer left and right hemisphere ##
    for line in DSSS:
        ## get varname
        varname = ''.join(line[2:6])
        varname = ''.join(varname.split('_')) # remove remaining underlines in varname
        varname = ''.join(varname.split('3')) # remove number3 of dil3 from ROI
        
        # \RetMapuseConnFalsepreConnFalserightGMdilECC
        
        ecc, pol = line[0], line[1]
        texcmd = '\\newcommand{\%sECC}{%.2f}\n' % (varname, ecc)
        outfile.write(texcmd) ## write to file
        texcmd = '\\newcommand{\%sPOL}{%.2f}\n' % (varname, pol)
        outfile.write(texcmd) ## write to file
        
    outfile.close()    
    print "Saved in %s" % outsavepath
    return outsavepath

########################################################################
def getLatexVars_DSSSanat(DSSS_anat, outsavepath):
    """
    get latex readable variables

    take the DSSS list and convert it into a variable list for tex
    ie. \newcommand{\OtheruseConnTruepreConnFalseleftGMdil3}{0.92}
    FRK2015
    """
    
    # TODO: 
    # - fix problem: multiple datapints from exactly same condition will get same varname 
    #     -> bad for tex!!
        
    # save to a texfile
    outfile = open(outsavepath,"w")

    for line in DSSS_anat:
        ## get varname
        varname = ''.join(line[2])
        varname = ''.join(varname.split('_')) # remove remaining underlines in varname
        varname = ''.join(varname.split('3')) # remove number3 of dil3 from ROI

        ecc, pol = line[0], line[1]
        texcmd = '\\newcommand{\Anat%sECC}{%.2f}\n' % (varname, ecc)
        outfile.write(texcmd) ## write to file
        texcmd = '\\newcommand{\Anat%sPOL}{%.2f}\n' % (varname, pol)
        outfile.write(texcmd) ## write to file
    # also mean
    eccMean = np.mean([DSSS_anat[1][3][0], DSSS_anat[0][3][0]]) # for the mean(eccL, eccR)
    polMean = np.mean([DSSS_anat[1][3][1], DSSS_anat[0][3][1]]) # for the mean(polL, polR)
    
    texcmd = '\\newcommand{\AnatBothMeanECC}{%.2f}\n' % (eccMean)
    outfile.write(texcmd) ## write to file
    texcmd = '\\newcommand{\AnatBothMeanPOL}{%.2f}\n' % (polMean)
    outfile.write(texcmd) ## write to file

    outfile.close()    

    return outsavepath

########################################################################
def convertDSSS2pandasDFrame(DSSS, DSSS_anat=[]):
    """
    FK2016: Take the DSSS dataset and convet it to pandas for the ie. barplots
    
    usage:
    DF = convertDSSS2pandasDFrame(DSSS)
    
    sns.barplot(x="method", y="cdist",  data=DF[DF['Mapper'] == 'RetMap']);
    """
    
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    import copy
    
    # merge DSSS with DSSS_anat 
    DSSS_hypal = copy.copy(DSSS)
    [DSSS_hypal.append(ds) for ds in DSSS_anat]
    
    DF = pd.DataFrame()
    
    for DS in DSSS_hypal:
        # exclude the mean because I will calculate them (in the plots) anyway
        ds = DS[2:] 
        mapper = ds[0]
        
        # get the hypal condition.
        # for anat replace some things because they are at different positions.
        c = ds[1:3]
        if c == ['useConnFalse', 'preConnFalse']: 
            cond = 'normal Hyperalignment'
            roi = ds[3]
            vals = ds[4]
        elif c == ['useConnTrue', 'preConnTrue']: 
            cond = 'PreConnectome'
            roi = ds[3]
            vals = ds[4]
        elif c == ['useConnTrue', 'preConnFalse']: 
            cond = 'Connectome'
            roi = ds[3]
            vals = ds[4]
        elif type(round(c[0][0][0])) is float: 
            cond = 'anatomical alignment'
            mapper = 'Anatomy'
            roi = ds[0]
            vals = ds[1]
        else: print "ERROR: check your Conditions (like preConnTrue, useConnTrue)"

        # number of values for ecc and for pol
        cdistVals = np.array(vals)
        nEle = cdistVals[1].shape[0] 

        # should give _sth_ (just rough!) like:
        # [RetMap noConn leftGM_dil3 0.5655546 ecc]
        d = np.vstack([
            zip( [mapper]*nEle, [cond]*nEle, [roi]*nEle, cdistVals[0,:], ['ecc']*nEle ),
            zip( [mapper]*nEle, [cond]*nEle, [roi]*nEle, cdistVals[1,:], ['pol']*nEle )
            ])
        df = pd.DataFrame(d,columns=['UsedDataset','method','ROI','cdist','ecc_pol'])

        # convert str2num to be able to calc in plots
        df['cdist'] = pd.to_numeric(df['cdist']) 

        # now attach the df to generate the allInOneDataFrame .. yeah!!
        DF = pd.concat([DF, df])

    # reset index. Otherwise it will just repeat the index for every sub-df
    DF.reset_index(drop=True)

    return DF
