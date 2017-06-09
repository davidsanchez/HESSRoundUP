import Loggin,numpy,string,os,glob
import RunListMaker

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler

class MergeMaker(Loggin.Message):
    def __init__(self,period= "test",cut="Stereo",Library="RoundUp/Library",Maps="RoundUp/Maps"):
        super(MergeMaker,self).__init__()
        Loggin.Message.__init__(self)

        self.info("Period "+period)
        self.info("Config "+cut)
        self.MergeTemplate    = os.environ["HESSUSER"]+"/RoundUpCode/gui/template/MergeTemplate"+cut
        self.clustered_runlists = None
        self.period = period
        self.cut = cut
        self.Library = os.environ["HESSUSER"]+"/"+Library
        self.Maps = os.environ["HESSUSER"]+"/"+Maps
        self.firstrun = 0
        os.system("mkdir -p "+self.Maps)
        os.system("mkdir -p "+self.Maps+"/"+self.period)

    def ClusterRun(self,runlist):
        RLMaker = RunListMaker.RunListMaker(runlist)
        RLMaker.LoadRunlist()
        RLMaker.SelectRunFromRunlist()
        centers  = []


        self.info("Clustering runs")

        for key in RLMaker.runlist_info:
           gg = 0
           ra = RLMaker.runlist_info[key]["RaTarget"]#+RLMaker.runlist_info[key]["RaOffset"]
           dec = RLMaker.runlist_info[key]["DecTarget"]#+RLMaker.runlist_info[key]["DecOffset"]

           if  ra> 180:
              gg= ra-360
           else:
              gg= ra
           centers.append([ra,dec])
   
# n_samples put to correspond to the number of runs in the file
# shuffle = 0, otherwise it randomize the positions...
        X, labels_true = make_blobs(n_samples=len(RLMaker.runlist_info), centers=centers, cluster_std=0.001, shuffle=0,
                            random_state=0)

        print('Length of X array = %i' % len(X))

# Compute DBSCAN
# eps = 2.0 should be the search radius
        db = DBSCAN(eps=3.0, min_samples=1, metric='euclidean').fit(X)
        core_samples_mask = numpy.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

        print('Estimated number of clusters: %d' % n_clusters_)
        print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
        print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
        print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))

##############################################################################
# Print results aself.periodnd make run lists
        from collections import defaultdict

        clustered_runlists = defaultdict(list)

        j = 0
        for key in RLMaker.runlist_info:
           clustered_runlists[labels[j]].append(int(key))
           j+=1

        self.runlist = RLMaker.runlist_info
        self.clustered_runlists = clustered_runlists

    def CheckCluster(self,firstrun = 0, doall = False):
       if self.clustered_runlists == None:
          self.warning("No cluster found")
          return
       for key in self.clustered_runlists.keys():
          print "\n"
          if max(self.clustered_runlists[key])<firstrun and not(doall):
             del self.clustered_runlists[key]	

          for i in xrange(len(self.clustered_runlists[key])):
             run = self.clustered_runlists[key][i]
             print run," ",self.runlist[run]["TargetName"]," ",self.runlist[run]["DecTarget"]," ",self.runlist[run]["RaTarget"]

    def plot(self):
        import matplotlib.pyplot as plt

        glong =  []
        glat =  []
        for key in self.runlist.keys():
          if self.runlist[key]["GLongTarget"] > 180:
             glong.append(self.runlist[key]["GLongTarget"]-360)
          else :
             glong.append(self.runlist[key]["GLongTarget"])
          glat.append(self.runlist[key]["GLatTarget"])
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection="hammer")
        ax.scatter(-numpy.radians(glong),numpy.radians(glat))
        ax.set_xticklabels(['150','120','90','60','30','0','-30','-60','-90','-120','-150'])
        ax.grid(True)
#        fig.show()
        fig.savefig("Obs_"+self.period+"_"+self.cut+".png")
        self.info("Fig save to Obs_"+self.period+"_"+self.cut+".png")
        os.system("display Obs_"+self.period+"_"+self.cut+".png")


    def MakeRootFileList(self,OnlyPeriod = True, doall = False):
#        get first run
        path = self.Library+"/"+self.period+"/*"+self.cut+"*.root"
        rl = []
        for f in glob.glob(path):
               rl.append(int(f[f.find(self.cut)+len(self.cut)+1:].replace(".root","\n")))
        self.firstrun = min(rl)
        print "first run of the period "+str(self.firstrun)

        path = self.Library

        if OnlyPeriod:
           path += "/"+self.period+"/*"+self.cut+"*.root"
           self.info("Merge only the current period "+self.period)
           out =  self.Library+"/"+self.period+"/"+self.cut+".list"

        else :
           path += "/*/*"+self.cut+"*.root"
           self.info("Merge all")
           out =  self.Library+"/Runlist_"+self.cut+".list"

        savefile = open(out,"w")
        for f in glob.glob(path):
          print f
          words = f[f.find(self.cut)+len(self.cut)+1:].replace(".root","\n")
          savefile.write(words)

        savefile.close()
        self.ClusterRun(out)
        self.CheckCluster(self.firstrun, doall=doall) 

    def ReadTemplate(self,template):
        template = open(template,"r").readlines()
#        template = open(self.AnalysisTemplate,"r").readlines()
        AT = ""
        for line in template:
             AT += line
        return AT

    def RunMerge(self,OnlyPeriod = True):
        out = self.Maps
        if OnlyPeriod:
           out+="/"+self.period+"/"

        self.info("save folder "+out)
        for key in self.clustered_runlists.keys():
            if  len(self.clustered_runlists[key]) == 0:
                continue
            AT = self.ReadTemplate(self.MergeTemplate)


            balise = '{0}'
            if OnlyPeriod:
               AT = AT.replace(balise,self.Maps+'/'+self.period+"/out")
               os.system("mkdir -p "+self.Maps+'/'+self.period+"/out")
            else:
               AT = AT.replace(balise,self.Maps+'/out')
               os.system("mkdir -p "+self.Maps+"/out")

            
            cmd = ""
            runlist = ""
            for run in  self.clustered_runlists[key]:
               target = self.runlist[run]["TargetName"].replace(" ","").replace("-","m").replace("+","p")+"_"+self.cut
               cmd+='ln -s '+self.Library+'/*/*'+str(self.cut)+'*'+str(run)+'*.root Results_'+target+"_run_"+str(run)+".root\n"
               runlist += "add_run("+str(int(run))+");\n"    

            balise = '{1}'
            AT = AT.replace(balise,cmd)

            balise = '{2}'
            AT = AT.replace(balise,str(key))
 
            balise = '{3}'
            AT = AT.replace(balise,target)

            balise = '{4}'
            AT = AT.replace(balise,out)

            balise = '{5}'
            AT = AT.replace(balise,runlist)

            self.info("Work on Target "+self.runlist[run]["TargetName"])

            if OnlyPeriod:
               path = self.Maps+'/'+self.period+"/"
               scriptname = "Merge"+str(int(key))+"_"+self.cut+".csh"
            else:
               path = self.Maps+"/Scripts/"
               scriptname = "Merge"+str(int(key))+"_"+self.cut+".csh"

            
            self.info("Saving Script to "+path+scriptname)
            script = open(path+scriptname,"w")
            script.write(AT)
            script.close()
            self.info("submiting Job")
            os.system("chmod +x "+path+scriptname)
            os.system("xterm  -e "+path+"./"+scriptname)

if __name__=='__main__':
   MM = MergeMaker("P2015_08",cut="Stereo")
   MM.MakeRootFileList(False)
   #MM.RunMerge(False)
#   MM.ClusterRun("allrun.list")
#   MM.CheckCluster(116000)

#   MM.plot()
