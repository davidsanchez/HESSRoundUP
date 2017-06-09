import Loggin,numpy,string,os


class ScriptMaker(Loggin.Message):
#    def __init__(self,Source,prod = "Prod24",out="/afs/in2p3.fr/home/d/dsanchez/sps",submit=False):
    def __init__(self,runlist_info,period,cut="Mono",Library="RoundUp/Library"):
        super(ScriptMaker,self).__init__()
        Loggin.Message.__init__(self)

        self.info("Period "+period)
        self.info("Config "+cut)

        self.AnalysisTemplate = os.environ["HESSUSER"]+"/RoundUpCode/gui/template/AnalyseTemplate"+cut

        self.runlist_info = runlist_info## contains runId and ra/dec and sourcename
        self.Library = os.environ["HESSUSER"]+"/"+Library
        self.period = period
        self.cut = cut
        os.system("mkdir -p "+self.Library+'/'+self.period)


    def ReadTemplate(self,template):
        template = open(template,"r").readlines()
#        template = open(self.AnalysisTemplate,"r").readlines()
        AT = ""
        for line in template:
             AT += line
        return AT


    def MakeAnalysisScript(self):
        self.info("Run analysis")

        for key in self.runlist_info :
            AT = self.ReadTemplate(self.AnalysisTemplate)

            balise = '{0}'
            AT = AT.replace(balise,self.runlist_info[key]["TargetName"])

            balise = '{1}'
            AT = AT.replace(balise,str(self.runlist_info[key]["RaTarget"]))

            balise = '{2}'
            AT = AT.replace(balise,str(self.runlist_info[key]["DecTarget"]))

            balise = '{3}'
            AT = AT.replace(balise,str(int(key)))

            balise = '{4}'
            AT = AT.replace(balise,self.Library+'/'+self.period+'/Results_run_'+self.cut+'_'+str(int(key))+'.root')

            balise = '{5}'
            AT = AT.replace(balise,self.Library+'/'+self.period+"/out/")
            os.system("mkdir -p "+self.Library+'/'+self.period+"/Scripts/")
            os.system("mkdir -p "+self.Library+'/'+self.period+"/out/")
            scriptname = self.Library+'/'+self.period+"/Scripts/analysis"+str(int(key))+"_"+self.cut+".csh"
            self.info("Saving Script to "+scriptname)
            script = open(scriptname,"w")
            script.write(AT)
            script.close()
            self.info("submiting Job")
            os.system("qsub "+scriptname)

        self.info("Finished")

    def CheckJobs(self,Nrun=0):
        os.system("qstat -u "+os.environ['USER']+" > joblist.dat")
        N = len(open("joblist.dat","r").readlines())
        if N>Nrun+2:
             self.info("Still waiting for "+str(N-2)+" job(s) to finish")
             return True
        return False


    def CheckAna(self,runana = False):
       folder = self.Library+'/'+self.period
       new_list = {}
       goodrun = 0
       for key in self.runlist_info :
           rootname = folder+"/Results_run_"+self.cut+"_"+str(key)+".root"
           isthere = os.path.isfile(rootname)
           if not(isthere):
                self.warning(rootname+' not found')
                new_list[key] = self.runlist_info[key]
           else:
                goodrun +=1
       self.success(str(goodrun)+" over "+str(len(self.runlist_info))+" have been successfully analysed")
       if runana:
          self.runlist_info = new_list
          self.MakeAnalysisScript()

if __name__=='__main__':
   import RunListMaker
   runlist = RunListMaker.RunListMaker(runlist = "../RoundUp/Library/P2015_07/P2015_07_Mono.list.list")
   runlist.LoadRunlist()
   runlist.SelectRunFromRunlist()
   runlist.SaveRunFromRunlist()


   SMaker = ScriptMaker(runlist.runlist_info,"P2015_07",cut="Mono")
   SMaker.MakeAnalysisScript()
