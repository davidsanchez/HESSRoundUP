import os,sys
import traceback
import gtk
import math

try:
    sys.path.append(os.environ["HESSUSER"])
except KeyError:
    pass
sys.path.append(os.environ["HESSROOT"])

from dbgui.db import MySQLconnection, SimpleTable
#from display.python import gtk_supp, ErrorDialogs, Errors, Helper, Storage

import Loggin,numpy,string,os


class RunListMaker(Loggin.Message):
    def __init__(self,runlist = "test.list",Library="",outrunlist=""):
        super(RunListMaker,self).__init__()
        Loggin.Message.__init__(self)
        self.outrunlist = outrunlist
        try :
            self.runlist = numpy.genfromtxt(runlist,unpack=True)
            #get the 1 row, assuming that this is the run number
            if self.runlist.size != len(self.runlist):
               self.runlist = self.runlist[0]
        except:
            self.error("Can fill the run list")
            self.runlist = []
 

    def LoadRunlist(self):
            self.info("Loading Run Quality list")
            # read the analysis section of the database
            ana_connection = MySQLconnection(config="analysis")
            ana_connection.connect()

            # Read run quality data to fill the RW database
            runquality_table = SimpleTable("RunQuality", ana_connection)
            runquality_set_data = runquality_table.get_sets(as_dictionary=1)
            runquality_data = runquality_table.read_data_as_dictionary(set = None,columns = runquality_table.fieldsnames + ["RunNumber",])

            ana_connection = MySQLconnection(config="analysis")
            ana_connection.connect()

            # Read run quality data to fill the RW database
            runquality_table = SimpleTable("RunQuality", ana_connection)
            runquality_set_data = runquality_table.get_sets(as_dictionary=1)
            runquality_data = runquality_table.read_data_as_dictionary(set = None,columns = runquality_table.fieldsnames + ["RunNumber",])

            # loop over the runquality_data
            self.runsQualityDB = {}
            mydic = {}
            for data in runquality_set_data:
                RunId = int(data['RunNumber'])
#                if RunId == 90106 or RunId == 90108:
#                   for k in data.keys():
#                        print k," ",data[k]
#                   print 
                self.runsQualityDB[RunId] = {}
                mydic[RunId] = {}
                mydic[RunId]["TargetName"]  = data["Target"]
                mydic[RunId]["RaTarget"]  = data["Nominal_Target_Ra"]-data["Nominal_WoobleOffset_Ra"]
                mydic[RunId]["DecTarget"] = data["Nominal_Target_Dec"]-data["Nominal_WoobleOffset_Dec"]
                mydic[RunId]["RaOffset"]  = data["Nominal_WoobleOffset_Ra"]
                mydic[RunId]["DecOffset"] = data["Nominal_WoobleOffset_Dec"]
                mydic[RunId]["TelsInRun"] = data["TelsInRun"]
                mydic[RunId]["RunDate"] = data["RunDate"]
                mydic[RunId]["RunDate"] = data["RunDate"]
                mydic[RunId]["GLongTarget"] = data["Nominal_Target_GLong"]
                mydic[RunId]["GLatTarget"] = data["Nominal_Target_GLat"]
                self.runsQualityDB[RunId] = mydic[RunId]

    def SelectRunFromRunlist(self):
        self.runlist_info = {}
        for RunId in self.runlist:
          if RunId in self.runsQualityDB.keys():
             self.runlist_info[int(RunId)] = self.runsQualityDB[RunId] 

    def SaveRunFromRunlist(self):
        self.SelectRunFromRunlist()
        if self.outrunlist == "":
             self.outrunlist = "RunList_info.list"
        fil = open(self.outrunlist,"w")
        for key  in self.runlist_info.iterkeys():
#           print str(self.runlist_info[key])
           fil.write(str(int(key))+" "+str(self.runlist_info[key])+"\n")
        fil.close()

if __name__=='__main__':
   runlist = RunListMaker(runlist = "../RoundUp/PKS2155_NewP.list")
   runlist.LoadRunlist()
   runlist.SelectRunFromRunlist()
   #runlist.SaveRunFromRunlist()

