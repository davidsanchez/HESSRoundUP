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

class MapsWindow(gtk.HBox):
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

        self.config = {}
        ##############################################################################

        self.fin = gtk.FileChooserButton("")
#        self.fin.set_current_folder("")
#TODO
        self.fin.set_current_folder("$HESSUSER/RoundUp/Maps/")
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
        self.target_position = True
        self.ra_mode = True
        self.ra = 329.7167
        self.dec = -30.2256
        self.UseConfigTarget = 1;

        f = gtk.Frame("Target ")
        h1.pack_start(f, expand=False, fill=False, padding=5)
        t3= gtk.Table(3, 3)
        f.add(t3)

        self.UseConfigTarget_button = SignalWidgets.CheckButton("Config from File results")
        self.UseConfigTarget_button.set_active(True)
        t3.attach(self.UseConfigTarget_button,0, 3, 1, 2, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.UseConfigTarget_button.connect("toggled", self.changed_UseConfigTarget) 




        t3.attach(gtk.Label("Target"), 0, 1, 2, 3, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.target_list = gtk_supp.EntryCompletion(0)
        self.target_list.set_size_request(150, -1)
        t3.attach(self.target_list, 1, 2, 2, 3, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.target_list.set_sensitive(False)

        

        self.filter_pos_button = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "Lookup")
        self.filter_pos_button.set_size_request(150, -1); 
        self.filter_pos_button.connect("clicked", self.lookup_pos)
        t3.attach(self.filter_pos_button, 2, 4, 2, 3, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)        
        self.filter_pos_button.set_sensitive(False)

        self.ra_label = gtk.Label("RA :")
        t3.attach(self.ra_label, 1, 2, 3,4, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.ra_label.set_sensitive(False)


        self.ra_entry = gtk.SpinButton(gtk.Adjustment(self.ra, -360, 360, .1), .5, 4)
        self.ra_entry.set_numeric(True)
        self.ra_entry.set_size_request(30, -1);
        t3.attach(self.ra_entry, 2, 3, 3,4, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1) 
        self.ra_entry.set_sensitive(False)


        self.dec_label = gtk.Label("DEC :")
        t3.attach(self.dec_label, 1, 2, 4,5, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.dec_label.set_sensitive(False)


        self.dec_entry = gtk.SpinButton(gtk.Adjustment(self.dec, -90, 90, .1), .5, 4)
        self.dec_entry.set_numeric(True)
        self.dec_entry.set_size_request(30, -1);
        t3.attach(self.dec_entry, 2, 3, 4,5, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1) 
        self.dec_entry.set_sensitive(False)

       
        self.ExtX_label = gtk.Label("Map Size X :")
        t3.attach(self.ExtX_label, 1, 2, 5,6, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)

        self.ExtX_entry = gtk.SpinButton(gtk.Adjustment(3, 0, 10, .1), .5, 4)
        self.ExtX_entry.set_numeric(True)
        self.ExtX_entry.set_size_request(30, -1);
        t3.attach(self.ExtX_entry, 2, 3, 5,6, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1) 


        self.ExtY_label = gtk.Label("Map Size Y :")
        t3.attach(self.ExtY_label, 1, 2, 6,7, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)

        self.ExtY_entry = gtk.SpinButton(gtk.Adjustment(3, 0, 10, .1), .5, 4)
        self.ExtY_entry.set_numeric(True)
        self.ExtY_entry.set_size_request(30, -1);
        t3.attach(self.ExtY_entry, 2, 3, 6,7, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1) 

        ##############################################################################


        #Start button
        self.exec_button = gtk_supp.PixmapButton(gtk.STOCK_EXECUTE, "Start Map Creation")
        self.exec_button.connect("clicked", self.start_map)
        h1.pack_start(self.exec_button, expand=False, fill=False, padding=1)


        ##############################################################################
        #Selection window
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



    ##
    ## Activation function, called when the GUI is fully loaded
    ##
    def activate(self):
        self.load_target_list()
#        self.filter_pos_button.connect("clicked", self.lookup_pos)
        self.target_list.connect("changed", self.changed_target);

    def load_run_list(self, widget=None, selection=None, * args):
        '''
        load file list from the folder
        '''
        self.mainwindow.busy_cursor()
        self.list.clear()

        try: 
            os.chdir(self.fin.get_current_folder())
            for f in glob.glob("Results_*.root"):
                siz = os.path.getsize(f)
                self.list.append([f,str(int(siz/1e6))+' MB'], info=f)
        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(x)

        self.mainwindow.normal_cursor()



    def start_map(self, *args):


        if not self.mainwindow.validate_batchsystem():
            return
        try:
            print self.list
            rootfiles = self.list.get_selection_info()
            if len(rootfiles) == 0:
                ErrorDialogs.Error("No ROOT files selected! Please load the list and select the ROOT files you want to process.")
            for rf in rootfiles:
               filename = "AnalyseRoundUpSource" 
               batchfilename = "start_%s"%(rf[:-5])

               self.config["char *Resfile"] = '"'+rf+'"'


               prefix = getPrefix(rf)

               self.config["char *OutfilePrefix"] = '"'+prefix+'"'


               batchfile = self.mainwindow.get_batchfile(batchfilename,filename,
                                                          [Storage.DIR_DST,
                                                           Storage.DIR_TABLES,
                                                           Storage.DIR_RESULTS,
                                                           Storage.DIR_CALIBRATION],
                                                          "HEGS")

               print batchfilename
               print filename
               #TODO official location of the code
               batchfile.write_line("ln -s $HESSUSER/RoundUpCode/RoundUpFancyMaps.C . ")
               batchfile.write_line("ln -s $HESSUSER/RoundUpCode/RoundUpFancyMaps.hh . ")
               batchfile.write_line("ln -s $HESSUSER/RoundUpCode/AnalyseRoundUpSource.C . ")
               batchfile.write_line("ln -s $HESSUSER/RoundUpCode/rootlogon.C . ")


               batchfile.write_line("ln -s "+self.fin.get_current_folder()+"/"+rf+" . ")
               root_opts =  "-q -b"
               if self.UseConfigTarget :
                    self.ra = 0
                    self.dec = 0
               batchfile.write_line("root %s '%s.C(\"%s\",\"%s\",%s,%s,%s)'" % (root_opts, filename,rf,prefix+'_roundup',self.UseConfigTarget,self.ra,self.dec))
               batchfile.write_line("rm "+rf)
               batchfile.write_line("cp *roundup*root "+self.fin.get_current_folder())
               batchfile.write_line("cp *roundup*png "+self.fin.get_current_folder())
#               batchfile.save_result_file(".", "*.root", "." , )
##               batchfile.save_result_file(".", "*.txt" , ".", "")
#               batchfile.save_result_file(".", "*.png" , ".",  "")
                
#               os.system("chmod +x "+batchfilename)
#               os.system("xterm  -e "+batchfilename)

               r, jobname = batchfile.submit(terminal=False)
        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(str(x))
        self.mainwindow.normal_cursor()


#    def Creatbatch(self):

    ##
    ## Sets the completion list of target
    ## @param targets list of available targets
    ##
    def set_target_list(self, targets):
        self.target_list.clear()

    ##
    ## Loads the list of observed targets from DB
    ##
    def load_target_list(self, *args):
#        try:
            targets = find_target.find_std_targets()
            self.set_target_list(targets)
            if not self.target_list.loaded_value_accessed:
                self.target_list.set_text(self.target_list.loaded_value)
                self.target_list.loaded_value_accessed = True
#        except Exception:
#            pass


    def set_target_list(self, targets):
        self.target_list.clear()
        self.target_list.set_completion_strings(targets)   

    ##
    ## Lookup target position in database
    ##
    def lookup_pos(self, *args):
#        try:
            target = self.target_list.get_text()
            var = find_target.find_target(target)
            if var is not None:
                if var["Name"] != target:
                    w = ErrorDialogs.Warning("Found approximate match %s"%(var["Name"]),timeout = 5000)
                    self.target_list.set_text(var["Name"])
                self.ra_entry.set_value(float(var["RA"]))
                self.dec_entry.set_value(float(var["Dec"]))
#                self.posmode_button.set_active(1)
#                self.parent.target_changed()
#            else:
#                w = ErrorDialogs.Error("Target %s not found in database" % (target))
#        except Exception:
#            traceback.print_exc()
#            w = ErrorDialogs.Error("Target %s not found in database"%(target))
   
def getPrefix(name):
    prefix= name.replace(".root","")
    prefix= prefix.replace("Results_","")
    prefix= prefix.replace("_merged","")
    return prefix
