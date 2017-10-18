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


class CalibCheckWindow(gtk.HBox):
    '''
    '''
    def __init__(self, mainwindow):
        '''
        Initialize the window
        @param mainwindow: reference to the parent window
        '''
        gtk.HBox.__init__(self)
        self.mainwindow = mainwindow
        self.config = mainwindow.config
        self.section = type(self).__name__
        h1 = gtk.VBox()

        ##############################################################################


        self.fin = gtk.FileChooserButton("file")
#        self.fin.set_current_folder("")
#TODO
        self.fin.set_current_folder("$HESSUSER/RoundUp/Library")
        self.fin.set_action(gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER)
        self.fin.show()

        # Simulated spectrum
        f = gtk.Frame("Root files")
        h1.pack_start(f, expand=False, fill=False, padding=5)
        t = gtk.Table(5, 2)
        f.add(t)
        t.attach(self.fin, 0, 2, 0, 1, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 5)
#        t.attach(gtk.Label("Spectral Index"), 0, 1, 1, 2, gtk.EXPAND | gtk.FILL, 5, 5)

        #Load file list button
        self.load_button = gtk_supp.PixmapButton(gtk.STOCK_OPEN, "Load file list")
        self.load_button.connect("clicked", self.load_run_list)
#        h1.pack_start(self.load_button, expand=False, fill=False, padding=1)
        t.attach(self.load_button, 2, 3, 0, 1, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 5)



       ##############################################################################
        # ==== default values ====

        #Start button
        self.exec_buttonCOGNSB = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "Check CoG maps And NSB")
        self.exec_buttonCOGNSB.connect("clicked", self.start_COGNSB)
        h1.pack_start(self.exec_buttonCOGNSB, expand=False, fill=False, padding=1)



        self.exec_BP = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "Broken Pixel")
        self.exec_BP.connect("clicked", self.start_BP)
        h1.pack_start(self.exec_BP, expand=False, fill=False, padding=1)



#        ##############################################################################
#        #Selection window
        self.pack_start(h1, expand=False, fill=True, padding=5)
        h1 = gtk.VBox()
        self.pack_start(h1, expand=True, fill=True, padding=5)
#        self.list = gtk_supp.SortList(1, (""))
        self.list = gtk_supp.SortList(2, ("filename","size"))
        self.list.set_selection_mode(gtk.SELECTION_EXTENDED)
        self.list.set_size_request(400, -1)
        self.list.set_column_width(0, 400)
        self.list.set_column_width(1, 100)
        swin = gtk.ScrolledWindow()
        swin.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        swin.set_size_request(-1, -1)
        swin.add(self.list)
        h1.pack_start(swin, expand=True, fill=True, padding=5)


    def load_run_list(self, widget=None, selection=None, * args):
        '''
        load file list from the folder
        '''
        self.mainwindow.busy_cursor()
        self.list.clear()

        try: 
            os.chdir(self.fin.get_current_folder())
            for f in glob.glob("Result*.root"):
                siz = os.path.getsize(f)
                self.list.append([f,str(int(siz/1e6))+' MB'], info=f)
        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(x)

        self.mainwindow.normal_cursor()
  

    def start_COGNSB(self, widget=None, selection=None, * args):
        '''
        run open_analysis_results.C
        '''

        if not self.mainwindow.validate_batchsystem():
            return
        try:


            rootfiles = self.list.get_selection_info()
            if len(rootfiles) == 0:
                ErrorDialogs.Error("No ROOT files selected! Please load the list and select the ROOT files you want to process.")
            if len(rootfiles) > 1:
                ErrorDialogs.Error("Select only 1 ROOT file.")

            prefix = getPrefix(rootfiles[0])
#            folder = self.mainwindow.get_batchsystem().storage.find_directory(Storage.DIR_RESULTS).get_base_directory()

            folder = self.fin.get_current_folder()


            batchfilename = "open_calibration_results_"
            filename = "open_calibration_results"
            batchfile = self.mainwindow.get_batchfile(batchfilename,filename,
                                                          [Storage.DIR_DST,
                                                           Storage.DIR_TABLES,
                                                           Storage.DIR_RESULTS,
                                                           Storage.DIR_CALIBRATION],
                                                          "HEGS")

               #TODO official location of the code
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/open_calibration_results.C . ")
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/rootlogon.C . ")
#            batchfile.write_line("ln -s "+self.fin.get_current_folder()+"/"+prefix+"_MapGenerationInfos.txt . ")


            root_opts =  "-l"
            batchfile.write_line("root %s '%s.C+(\"%s\")'" % (root_opts, "open_calibration_results",folder+"/"+rootfiles[0]))
#               batchfile.save_result_file(".", "*.root", ".")
                
            r, jobname = batchfile.submit(terminal=True)

        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(str(x))
        self.mainwindow.normal_cursor()



    def start_BP(self, widget=None, selection=None, * args):
        '''
        run MakeSummary.C
        '''

        if not self.mainwindow.validate_batchsystem():
            return
        try:
            rootfiles = self.list.get_selection_info()
            if len(rootfiles) == 0:
                ErrorDialogs.Error("No ROOT files selected! Please load the list and select the ROOT files you want to process.")

            filename = open("rootfiles.list","w")

            for rf in rootfiles:
                filename.write(rf+"\n")

            filename.close()

            folder = self.fin.get_current_folder()


            batchfilename = "BpxViewer_"
            filename = "BpXViewer"
            batchfile = self.mainwindow.get_batchfile(batchfilename,filename,
                                                          [Storage.DIR_DST,
                                                           Storage.DIR_TABLES,
                                                           Storage.DIR_RESULTS,
                                                           Storage.DIR_CALIBRATION],
                                                          "HEGS")

               #TODO official location of the code
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/bpxviewer.C  . ")
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/rootlogon.C . ")

            prefix = getPrefix(rootfiles[0])
            runnumber = prefix.replace("run_","").replace("Mono_","").replace("Stereo_","")

            root_opts =  "-l"
            batchfile.write_line("root %s '%s.C+(%s)'" % (root_opts, "bpxviewer",runnumber))
                
            r, jobname = batchfile.submit(terminal=True)
        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(str(x))
        self.mainwindow.normal_cursor()





    ##
    ## Lookup target position in database
    ##
    def lookup_pos(self, *args):
#        try:
            print "debug"
            target = self.target_list.get_text()
            var = find_target.find_target(target)
            if var is not None:
                if var["Name"] != target:
                    w = ErrorDialogs.Warning("Found approximate match %s"%(var["Name"]),timeout = 5000)
                    self.target_list.set_text(var["Name"])
                self.ra_entry.set_value(float(var["RA"]))
                self.dec_entry.set_value(float(var["Dec"]))


def getPrefix(name):
    prefix= name.replace(".root","")
    prefix= prefix.replace("Results_","")
    prefix= prefix.replace("_merged","")
    return prefix	
