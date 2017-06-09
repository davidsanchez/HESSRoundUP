# author David Sanchez david.sanchez@lapp.in2p3.fr

'''
'''
import os,sys,glob
import traceback
import gtk,numpy
import math

from parisrunquality.gui import find_target
from dbgui.db import MySQLconnection, SimpleTable
from display.python import gtk_supp, ErrorDialogs, Errors, Helper, Storage, SignalWidgets
import ClusterAndMerge as CM


class AnalyseWindow(gtk.VBox):
    '''
    '''
    def __init__(self, mainwindow):
        '''
        Initialize the window
        @param mainwindow: reference to the parent window
        '''
        gtk.VBox.__init__(self)
        self.mainwindow = mainwindow
        self.config = mainwindow.config
        self.section = type(self).__name__
#        h1 = gtk.VBox()

#        frame = gtk.Frame(label="Analysis")
#        self.pack_start(frame, expand=False, fill=False, padding=0)
#        v = gtk.VBox()
#        frame.add(v)
        frame = gtk.Frame("Run analysis")


        table = gtk.Table(3,3, True)
        table.set_size_request(200, 150)
        table.show()

        frame.add(table)
        
        self.pack_start(frame, expand=False, fill=False, padding=0)
#        v = gtk.VBox()
        frame.show()


        label = gtk.Label("run list")
        table.attach(label, 0, 1, 0,1, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)

        self.runlist = gtk.FileChooserButton("run list")
        self.runlist.set_title("run list")
        self.runlist.set_filename("")
        self.runlist.show()

        table.attach(self.runlist, 1,2, 0,1, gtk.FILL,  gtk.FILL, 5, 1)

        plotButton = gtk.Button("Draw Obs Plot")
        plotButton.connect("clicked", self.plot, "")
        table.attach(plotButton,2,3, 0,1)
        plotButton.show()

        labelconfig = gtk.Label("Configuration")
        table.attach(labelconfig, 0,1 ,1,2, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)

        self.cut = "Mono"
        self.buttonMono = gtk.RadioButton(None, "Mono Analyis")
        self.buttonMono.connect("toggled", self.Changecut, "Mono")
        table.attach(self.buttonMono, 1,2 ,1,2, gtk.EXPAND | gtk.FILL, 5, 5)

        self.buttonStereo = gtk.RadioButton(self.buttonMono, "Stereo Analyis")
        self.buttonStereo.connect("toggled", self.Changecut, "Stereo")
        table.attach(self.buttonStereo, 2, 3,  1,2, gtk.EXPAND | gtk.FILL, 5, 5)

        labelPeriod = gtk.Label("Period")
        table.attach(labelPeriod, 0,1 ,2,3, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)

        self.period = gtk.Entry()
        self.period.set_max_length(50)
        self.period.connect("activate", self.enter_callback,  self.period)
        self.period.set_text("ChooseCurrentPeriod")
        self.period.select_region(0, len( self.period.get_text()))
        table.attach(self.period, 1, 2,  2,3, gtk.EXPAND | gtk.FILL, 5, 5)
        self.period.show()




####################################################################
        frame11 = gtk.Frame("Run")
        table11 = gtk.Table(2,3, True)
        table11.set_size_request(200, 150)
        table11.show()
        frame11.add(table11)
        self.pack_start(frame11, expand=False, fill=False, padding=0)
#        v = gtk.VBox()
        table11.show()

        RunButton = gtk.Button("Run Analysis")
        RunButton.connect("clicked", self.runanalysis, "")
        table11.attach(RunButton,1,2,0,1)
        RunButton.show()

        checkButton = gtk.Button("Check Analysis")
        checkButton.connect("clicked", self.Checkanalysis, "")
        table11.attach(checkButton,1, 2,1,2)
        checkButton.show()

        checkrunButton = gtk.Button("Check and re-Run Analysis")
        checkrunButton.connect("clicked", self.CheckanalysisAndRerun, "")
        table11.attach(checkrunButton,2, 3,1,2)
        checkrunButton.show()


####################################################################
        frame2 = gtk.Frame("Cluster and Merge")


        table2 = gtk.Table(2, 3, True)
        table2.set_size_request(200, 100)
        table2.show()

        frame2.add(table2)
        
        self.pack_start(frame2, expand=False, fill=False, padding=0)
#        v = gtk.VBox()
        frame2.show()

        labelmergePeriod = gtk.Label("Merge Period or All")
        table2.attach(labelmergePeriod, 0,1, 1,2, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)



        self.OnlyPerido = True
        self.buttonOnlyPerido = gtk.RadioButton(None, "Merge Only Current Period")
        self.buttonOnlyPerido.connect("toggled", self.ChangeMerge, True)
        table2.attach(self.buttonOnlyPerido, 1,2 ,1,2, gtk.EXPAND | gtk.FILL, 5, 5)

        self.buttonAll = gtk.RadioButton(self.buttonOnlyPerido, "Merge All data ")
        self.buttonAll.connect("toggled", self.ChangeMerge, False)
        table2.attach(self.buttonAll, 2, 3,  1,2, gtk.EXPAND | gtk.FILL, 5, 5)

        self.Doall_button = SignalWidgets.CheckButton("Merge Past Obs too")
        self.Doall_button.set_active(False)
        table2.attach(self.Doall_button,1,2, 3,4, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.Doall_button.connect("toggled", self.changed_Doall)


        MergeButton = gtk.Button("Run Merging")
        MergeButton.connect("clicked", self.Merge, "")
        table2.attach(MergeButton,2, 3,3,4)
        MergeButton.show()

    def changed_Doall(self, widget, data = None):
        self.Doall_button.set_active(not(self.Doall_button.get_active()))

    def plot(self, widget, data = None):
        self.mainwindow.busy_cursor()
        MM = CM.MergeMaker(self.period.get_text(),cut=self.cut)
        MM.ClusterRun(self.runlist.get_filename())
        MM.CheckCluster()
        self.mainwindow.normal_cursor()
        MM.plot()

        
    def Merge(self, widget, data = None):
        self.mainwindow.busy_cursor()
        MM = CM.MergeMaker(self.period.get_text(),cut=self.cut)
        MM.MakeRootFileList(self.OnlyPerido,self.Doall_button.get_active())
        MM.RunMerge(self.OnlyPerido)
        self.mainwindow.normal_cursor()

    def ChangeMerge(self, widget, data=None):
       self.OnlyPerido = data

    def Changecut(self, widget, data=None):
       if data=="Mono":
          self.cut = "Mono"
       if data=="Stereo":
          self.cut = "Stereo"

    def enter_callback(self, widget, entry):
        entry_text = entry.get_text()

    def runanalysis(self,widget,data=None):
          self.mainwindow.busy_cursor()
          import RunListMaker
          try :
             runlist = RunListMaker.RunListMaker(runlist = self.runlist.get_filename())
          except:
             pass
          runlist.LoadRunlist()
          runlist.SelectRunFromRunlist()
          runlist.SaveRunFromRunlist()

          import ScriptMaker
          SMaker = ScriptMaker.ScriptMaker(runlist.runlist_info,self.period.get_text().replace(" ",""),cut=self.cut)
          SMaker.MakeAnalysisScript()
          self.mainwindow.normal_cursor()


    def Checkanalysis(self,widget,data=None):
          self.mainwindow.busy_cursor()
          import RunListMaker
          try :
             runlist = RunListMaker.RunListMaker(runlist = self.runlist.get_filename())
          except:
             pass
          runlist.LoadRunlist()
          runlist.SelectRunFromRunlist()
          runlist.SaveRunFromRunlist()

          import ScriptMaker
          SMaker = ScriptMaker.ScriptMaker(runlist.runlist_info,self.period.get_text().replace(" ",""),cut=self.cut)
          SMaker.CheckAna()
          self.mainwindow.normal_cursor()


    def CheckanalysisAndRerun(self,widget,data=None):
          self.mainwindow.busy_cursor()
          import RunListMaker
          try :
             runlist = RunListMaker.RunListMaker(runlist = self.runlist.get_filename())
          except:
             pass
          runlist.LoadRunlist()
          runlist.SelectRunFromRunlist()
          runlist.SaveRunFromRunlist()

          import ScriptMaker
          SMaker = ScriptMaker.ScriptMaker(runlist.runlist_info,self.period.get_text().replace(" ",""),cut=self.cut)
          SMaker.CheckAna(True)
          self.mainwindow.normal_cursor()

##############################################################################

#        self.fin = gtk.FileChooserButton("file")
#        self.fin.set_current_folder("$HOME")
#        self.fin.set_action(gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER)
#        self.fin.show()
#        f = gtk.Frame("Root files")
#        f.add(v)
#        t = gtk.Table(5, 2)
#        f.add(t)
#        t.attach(self.fin, 0, 2, 0, 1, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 5)

#        #Load file list button
#        self.load_button = gtk_supp.PixmapButton(gtk.STOCK_OPEN, "Load file list")
#        self.load_button.connect("clicked", self.load_run_list)
#        t.attach(self.load_button, 2, 3, 0, 1, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 5)

#        self.dec_label = gtk.Label("DEC :")
#        t.attach(self.dec_label, 1, 2, 4,5, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)



    def load_run_list(self, widget=None, selection=None, * args):
        '''
        load file list from the folder
        '''
        self.mainwindow.busy_cursor()
        self.list.clear()
