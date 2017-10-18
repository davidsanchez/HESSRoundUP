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


class ToolsWindow(gtk.HBox):
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
        self.fin.set_current_folder("$HESSUSER/RoundUp/Maps")
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
        self.ra_mode = True
        self.ra = 0.
        self.dec = 0.
        self.MinSig = 4.5
        self.UseConfigTarget = 1
  
        f = gtk.Frame("Target ")
        h1.pack_start(f, expand=False, fill=False, padding=5)
        t3= gtk.Table(3, 3)
        f.add(t3)

        self.UseConfigTarget_button = SignalWidgets.CheckButton("Config from File results")
        self.UseConfigTarget_button.set_active(True)
        t3.attach(self.UseConfigTarget_button,0, 3, 1, 2, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.UseConfigTarget_button.connect("toggled", self.changed_UseConfigTarget) 




        t3.attach(gtk.Label("Target Name"), 0, 1, 2, 3, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.target_list = gtk_supp.EntryCompletion(0)
        self.target_list.set_size_request(150, -1)
        t3.attach(self.target_list, 1, 2, 2, 3, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        

        self.filter_pos_button = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "Lookup")
        self.filter_pos_button.set_size_request(150, -1); 
        self.filter_pos_button.connect("clicked", self.lookup_pos)
        t3.attach(self.filter_pos_button, 2, 4, 2, 3, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)        

        self.ra_label = gtk.Label("RA :")
        t3.attach(self.ra_label, 1, 2, 3,4, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)


        self.ra_entry = gtk.SpinButton(gtk.Adjustment(self.ra, -360, 360, .1), .5, 4)
        self.ra_entry.set_numeric(True)
        self.ra_entry.set_size_request(30, -1);
        t3.attach(self.ra_entry, 2, 3, 3,4, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1) 


        self.dec_label = gtk.Label("DEC :")
        t3.attach(self.dec_label, 1, 2, 4,5, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)


        self.dec_entry = gtk.SpinButton(gtk.Adjustment(self.dec, -90, 90, .1), .5, 4)
        self.dec_entry.set_numeric(True)
        self.dec_entry.set_size_request(30, -1);
        t3.attach(self.dec_entry, 2, 3, 4,5, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1) 

        self.target_list.set_sensitive(False)
        self.filter_pos_button.set_sensitive(False)
        self.ra_label.set_sensitive(False)
        self.ra_entry.set_sensitive(False)
        self.dec_label.set_sensitive(False)
        self.dec_entry.set_sensitive(False)

        #Start button
        self.exec_button = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "ONOFF Test @ Pos")
        self.exec_button.connect("clicked", self.start_ONOFFTest)
        # h1.pack_start(self.exec_button, expand=False, fill=False, padding=1)
        t3.attach(self.exec_button, 0, 2,5,6, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)


        #Start button
        self.exec_buttonFP = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "Info @ Pos")
        self.exec_buttonFP.connect("clicked", self.start_findPos)
        # h1.pack_start(self.exec_buttonFP, expand=False, fill=False, padding=1)
        t3.attach(self.exec_buttonFP, 2,4,5,6, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)


        f2 = gtk.Frame("Tools ")
        h1.pack_start(f2, expand=False, fill=False, padding=5)
        ttools= gtk.Table(3, 2)
        f2.add(ttools)


        # #Start button
        # self.exec_buttonFP = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "Start FindPosition")
        # self.exec_buttonFP.connect("clicked", self.start_findPos)
        # # h1.pack_start(self.exec_buttonFP, expand=False, fill=False, padding=1)
        # ttools.attach(self.exec_buttonFP, 1, 2, 4,5, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)


        #Start button
        self.exec_buttonFindspot = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "Find Hot Spot (Significance) ")
        self.exec_buttonFindspot.connect("clicked", self.start_FindHotSpot)
        # h1.pack_start(self.exec_buttonFindspot, expand=False, fill=False, padding=1)
        ttools.attach(self.exec_buttonFindspot, 1,2, 2,3, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)


        self.minsig_label = gtk.Label("Min Significance :")
        ttools.attach(self.minsig_label, 2,3, 2,3, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)



        self.MinSig_entry = gtk.SpinButton(gtk.Adjustment(self.MinSig, 0, 90, .1), .1, 4)
        self.MinSig_entry.set_numeric(True)
        self.MinSig_entry.set_size_request(30, -1);
        ttools.attach(self.MinSig_entry, 3,4, 2,3, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1) 


        #Start button
        self.exec_buttonFindspot = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "Find Hot Spot (ONOFF) ")
        self.exec_buttonFindspot.connect("clicked", self.start_FindHotSpotONOFF)
        # h1.pack_start(self.exec_buttonFindspot, expand=False, fill=False, padding=1)
        ttools.attach(self.exec_buttonFindspot, 1,2, 3,4, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)


        #Start button
        self.exec_buttonDrawMap = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "Open Results Maps")
        self.exec_buttonDrawMap.connect("clicked", self.start_OpenResults)
        h1.pack_start(self.exec_buttonDrawMap, expand=False, fill=False, padding=1)


#        tsum= gtk.Table(3, 3)
#        f.add(tsum)

        self.exec_MakeSummary = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "Make Summary")
        self.exec_MakeSummary.connect("clicked", self.start_MakeSummary)
        h1.pack_start(self.exec_MakeSummary, expand=False, fill=False, padding=1)



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


    def changed_UseConfigTarget(self, *args):
        self.UseConfigTarget = 0;
        if self.UseConfigTarget_button.get_active() :
              self.UseConfigTarget = 1;
              self.target_list.set_sensitive(False)
              self.filter_pos_button.set_sensitive(False)
              self.ra_label.set_sensitive(False)
              self.ra_entry.set_sensitive(False)
              self.dec_label.set_sensitive(False)
              self.dec_entry.set_sensitive(False)
        else:
              self.target_list.set_sensitive(True)
              self.filter_pos_button.set_sensitive(True)
              self.ra_label.set_sensitive(True)
              self.ra_entry.set_sensitive(True)
              self.dec_label.set_sensitive(True)
              self.dec_entry.set_sensitive(True)


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
  

    def start_findPos(self, widget=None, selection=None, * args):
        '''
        run findPosition.C
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
            folder = self.fin.get_current_folder()

            batchfilename = "FindPosition_"
            filename = "FindPosition"
            batchfile = self.mainwindow.get_batchfile(batchfilename,filename,
                                                          [Storage.DIR_DST,
                                                           Storage.DIR_TABLES,
                                                           Storage.DIR_RESULTS,
                                                           Storage.DIR_CALIBRATION],
                                                          "HEGS")

               #TODO official location of the code
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/FindPosition.C . ")
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/rootlogon.C . ")
#            batchfile.write_line("ln -s "+self.fin.get_current_folder()+"/"+prefix+"_MapGenerationInfos.txt . ")


            root_opts =  "-q -b"
            ra = self.ra_entry.get_value();
            dec = self.dec_entry.get_value();
            batchfile.write_line("root %s '%s.C(%s,%s,\"%s\",%s,%s)'" % (root_opts, "FindPosition",ra,dec,folder+"/"+prefix,self.UseConfigTarget,rootfiles[0]))
#               batchfile.save_result_file(".", "*.root", ".")
                
            r, jobname = batchfile.submit(terminal=True)
        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(str(x))
        self.mainwindow.normal_cursor()


    def start_ONOFFTest(self, widget=None, selection=None, * args):
        '''
        run DrawONOFFTest_atPosition.C
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

            batchfilename = "ONOFFTest_"
            filename = "ONOFFTest"
            batchfile = self.mainwindow.get_batchfile(batchfilename,filename,
                                                          [Storage.DIR_DST,
                                                           Storage.DIR_TABLES,
                                                           Storage.DIR_RESULTS,
                                                           Storage.DIR_CALIBRATION],
                                                          "HEGS")

               #TODO official location of the code
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/DrawONOFFTest_atPosition.C . ")
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/rootlogon.C . ")
#            batchfile.write_line("ln -s "+self.fin.get_current_folder()+"/"+prefix+"_MapGenerationInfos.txt . ")


            root_opts =  "-l"
            onoff = folder+"/"+prefix+'_roundup_ONOFFTestMaps.root'

            ra = self.ra_entry.get_value();
            dec = self.dec_entry.get_value();
            batchfile.write_line("root %s '%s.C++(\"%s\",\"%s\",%s,%s,%s,\"%s\",%s)'" % (root_opts, "DrawONOFFTest_atPosition",folder+"/"+rootfiles[0],onoff,ra,dec,"true",folder+"/"+prefix,self.UseConfigTarget))
#               batchfile.save_result_file(".", "*.root", ".")
                
            r, jobname = batchfile.submit(terminal=True)

        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(str(x))
        self.mainwindow.normal_cursor()


    def start_FindHotSpot(self, widget=None, selection=None, * args):
        '''
        run hotspotposition.C
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


            batchfilename = "FindHotSpot_"
            filename = "FindHotSpot"
            batchfile = self.mainwindow.get_batchfile(batchfilename,filename,
                                                          [Storage.DIR_DST,
                                                           Storage.DIR_TABLES,
                                                           Storage.DIR_RESULTS,
                                                           Storage.DIR_CALIBRATION],
                                                          "HEGS")

               #TODO official location of the code
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode//hotspotposition.C . ")
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/rootlogon.C . ")
#            batchfile.write_line("ln -s "+self.fin.get_current_folder()+"/"+prefix+"_MapGenerationInfos.txt . ")


            root_opts =  "-l"
            MinSig = self.MinSig_entry.get_value();
            batchfile.write_line("root %s '%s.C+(\"%s\",%d,%i)'" % (root_opts, "hotspotposition",folder+"/"+prefix,MinSig,1))
#               batchfile.save_result_file(".", "*.root", ".")
                
            r, jobname = batchfile.submit(terminal=True)

        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(str(x))
        self.mainwindow.normal_cursor()


    def start_FindHotSpotONOFF(self, widget=None, selection=None, * args):
        '''
        run hotspotposition.C
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


            batchfilename = "FindHotSpot_"
            filename = "FindHotSpot"
            batchfile = self.mainwindow.get_batchfile(batchfilename,filename,
                                                          [Storage.DIR_DST,
                                                           Storage.DIR_TABLES,
                                                           Storage.DIR_RESULTS,
                                                           Storage.DIR_CALIBRATION],
                                                          "HEGS")

               #TODO official location of the code
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode//hotspotposition.C . ")
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/rootlogon.C . ")
#            batchfile.write_line("ln -s "+self.fin.get_current_folder()+"/"+prefix+"_MapGenerationInfos.txt . ")


            root_opts =  "-l"
            MinSig = self.MinSig_entry.get_value();
            batchfile.write_line("root %s '%s.C+(\"%s\",%d,%i)'" % (root_opts, "hotspotposition",folder+"/"+prefix,MinSig,0))
#               batchfile.save_result_file(".", "*.root", ".")
                
            r, jobname = batchfile.submit(terminal=True)

        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(str(x))
        self.mainwindow.normal_cursor()

    def start_OpenResults(self, widget=None, selection=None, * args):
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


            batchfilename = "open_analysis_results_"
            filename = "open_analysis_results"
            batchfile = self.mainwindow.get_batchfile(batchfilename,filename,
                                                          [Storage.DIR_DST,
                                                           Storage.DIR_TABLES,
                                                           Storage.DIR_RESULTS,
                                                           Storage.DIR_CALIBRATION],
                                                          "HEGS")

               #TODO official location of the code
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/open_analysis_results.C . ")
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/rootlogon.C . ")
#            batchfile.write_line("ln -s "+self.fin.get_current_folder()+"/"+prefix+"_MapGenerationInfos.txt . ")


            root_opts =  "-l"
            batchfile.write_line("root %s '%s.C+(\"%s\",\"%s\")'" % (root_opts, "open_analysis_results",folder+"/"+rootfiles[0],folder+"/"+prefix))
#               batchfile.save_result_file(".", "*.root", ".")
                
            r, jobname = batchfile.submit(terminal=True)

        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(str(x))
        self.mainwindow.normal_cursor()



    def start_MakeSummary(self, widget=None, selection=None, * args):
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


            batchfilename = "makesummary_"
            filename = "makesummary"
            batchfile = self.mainwindow.get_batchfile(batchfilename,filename,
                                                          [Storage.DIR_DST,
                                                           Storage.DIR_TABLES,
                                                           Storage.DIR_RESULTS,
                                                           Storage.DIR_CALIBRATION],
                                                          "HEGS")

               #TODO official location of the code

            for rf in rootfiles:
                     batchfile.write_line("ln -s "+folder+"/"+rf+" . ")

            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/MakeSummary.C . ")
            batchfile.write_line("ln -s $HESSUSER/RoundUpCode/rootlogon.C . ")
            batchfile.write_line("ln -s "+folder+"/rootfiles.list . ")


            root_opts =  "-l"
            batchfile.write_line("root %s 'MakeSummary.C+(\"rootfiles.list\")'")
                
            r, jobname = batchfile.submit(terminal=True)
            batchfile.write_line("rm rootfiles.list")
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
