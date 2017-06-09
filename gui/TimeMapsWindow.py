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

class TimeMapsWindow(gtk.HBox):
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

        self.fin = gtk.FileChooserButton("Config file")
#        self.fin.set_current_folder("")
#TODO
        self.fin.set_current_folder("$HOME")
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
  
        #  //Configure source position 
        self.config["bool UseConfigTarget"] = 1;
        self.config["double UserLambda"] = self.ra;
        self.config["double UserBeta"] = self.dec;

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
        f = gtk.Frame("Analyis ")
        h1.pack_start(f, expand=False, fill=False, padding=5)
        t6 = gtk.Table(3, 3)
        f.add(t6)
        self.config["char *Regionsfile"] = os.getenv("HESSUSER")+"/HEGS/ExcRegionsRadiusNew.txt";

        self.config["char *AnalysisConfig"] = "Loose";

        self.profile_combo = gtk_supp.Combo()
        self.profile_combo.set_popdown_strings(["Loose","Std","Faint"])
        t6.attach(gtk.Label("Profile"), 0, 1, 0,  1, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 25, 1)
        self.profile_combo.set_size_request(200, -1); 
        t6.attach(self.profile_combo, 1, 4, 0,  1, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.profile_combo.set_active_text_no_signal(self.config["char *AnalysisConfig"])

        self.Exclusion_file = SignalWidgets.FileChooserButton("Exclusion Regions file")
        self.Exclusion_file.set_local_only(True)
        self.Exclusion_file.set_filename(self.config["char *Regionsfile"])
        self.Exclusion_file.set_size_request(200, -1)
        # open existing file only
        self.Exclusion_file.set_action(gtk.FILE_CHOOSER_ACTION_OPEN)
        myfilter = gtk.FileFilter()
        myfilter.add_pattern("*.txt")
        self.Exclusion_file.set_filter(myfilter)
        t6.attach(gtk.Label("Exclusion Region"), 0, 1, 1, 2, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 25, 1)
        t6.attach(self.Exclusion_file, 1, 2, 1, 2, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)

        self.RunNumberMin_label = gtk.Label("RunNumber Min :")
        t6.attach(self.RunNumberMin_label, 0, 1, 5,6, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)

        self.RunNumberMin_entry = gtk.SpinButton(gtk.Adjustment(-1, -1, 1000000, 1),  1, 0)
        self.RunNumberMin_entry.set_numeric(True)
        self.RunNumberMin_entry.set_size_request(30, -1);
        t6.attach(self.RunNumberMin_entry, 1, 3, 5,6, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1) 


        self.RunNumberMax_label = gtk.Label("RunNumber Max")
        t6.attach(self.RunNumberMax_label, 0, 1, 6,7, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)

        self.RunNumberMax_entry = gtk.SpinButton(gtk.Adjustment(1000000,  -1, 1000000, 1),  1, 0)
        self.RunNumberMax_entry.set_numeric(True)
        self.RunNumberMax_entry.set_size_request(30, -1);
        t6.attach(self.RunNumberMax_entry, 1, 3, 6,7, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1) 
       ##############################################################################
        f = gtk.Frame("Acceptance ")
        h1.pack_start(f, expand=False, fill=False, padding=5)
        t = gtk.Table(3, 3)
        f.add(t)

        self.config["bool EventsAndAcceptanceFromFile"] = 0
        self.config["bool EventsAndAcceptanceFromRadial"] = 1
        self.config["bool EventsAndAcceptanceFromRadialFile"] = 0


        self.button2 = gtk.RadioButton(None, "EventsAndAcceptanceFromRadial")
        self.button2.connect("toggled", self.acceptance, "EventsAndAcceptanceFromRadial")
#        t.attach(gtk.Label("EventsAndAcceptanceFromFile"), 0, 2, 1,2, gtk.EXPAND | gtk.FILL, 5, 5)
        t.attach(self.button2, 2, 3, 0,1, gtk.EXPAND | gtk.FILL, 5, 5)

        self.button1 = gtk.RadioButton(self.button2, "EventsAndAcceptanceFromFile")
        self.button1.connect("toggled", self.acceptance, "EventsAndAcceptanceFromFile")
#        t.attach(gtk.Label("EventsAndAcceptanceFromFile"), 0, 2, 0,1, gtk.EXPAND | gtk.FILL, 5, 5)
        t.attach(self.button1, 2, 3,  1,2, gtk.EXPAND | gtk.FILL, 5, 5)


        self.button3 = gtk.RadioButton(self.button2, "EventsAndAcceptanceFromRadialFile")
        self.button3.connect("toggled", self.acceptance, "EventsAndAcceptanceFromRadialFile")
#        t.attach(gtk.Label("EventsAndAcceptanceFromFile"), 0, 2,2,3, gtk.EXPAND | gtk.FILL, 5, 5)
        t.attach(self.button3, 2, 3, 2,3, gtk.EXPAND | gtk.FILL, 5, 5)

        ##############################################################################
        f = gtk.Frame("Exclusion ")
        h1.pack_start(f, expand=False, fill=False, padding=5)
        t2 = gtk.Table(3, 3)
        f.add(t2)

        self.config["bool ExclusionFromFits"] = 0
        self.config["bool ExclusionFromRegionFile"] = 0
        self.button4 = gtk.RadioButton(None, "ExclusionFromFits")
        self.button4.connect("toggled", self.exclusion, "ExclusionFromFits")
        t2.attach(self.button4, 2, 3, 0,1, gtk.EXPAND | gtk.FILL, 5, 5)

        self.button5 = gtk.RadioButton(self.button4, "ExclusionFromRegionFile")
        self.button5.connect("toggled", self.exclusion, "ExclusionFromRegionFile")
        t2.attach(self.button5, 2, 3, 1,2, gtk.EXPAND | gtk.FILL, 5, 5)
        ##############################################################################
        #  //Configure Safe threshold cut
        self.config["bool ApplySafeThreshold"] = 1;
        self.config["bool SafeThresholdFromAcceptance"] = 1;
        self.config["bool SafeThresholdFromPsiCut"] = 1;
        self.config["double SafeThresholdFromPsiValue"] = 0;
        self.config["double SafeThresholdRatioParam"] = 0.15;

        #  //Configure Run range and runs to forbid or to authorize 
        self.config["int RunNumberMin"] = 0;
        self.config["int RunNumberMax"] = -1;
        self.config["bool UseRunListToMatch"] = 0;
        self.config["bool UseRunsToForbid"] = 0;
        self.config["char *RunListToMatch"] = ""; 
        self.config["char *RunsToForbid"] = '"'+'RunsToForbid.txt'+'"';
        #  //Configure map parameters
        self.config["bool UseConfigMapParams"] = 0;
        self.config["double BinSize"] = 0.02;
        self.config["double PsiCut"] = 2.0;
        self.config["bool SelectAllRunsContributingToTheMap"] = 1;

        #  //Configure ring background 
        self.config["bool AdaptFFT"] = 1; 
        self.config["bool AdaptCut_Alpha"] = 1;
        self.config["bool ConstantArea"] = 0;
        self.config["bool ConstantThickness"] = 1; 
        self.config["bool ConstantInnerRadius"] = 0; 
        self.config["bool SmoothRings"] = 1; 
        self.config["double RingMinimalRadius"] = 0.7;
        self.config["double InnerRingMax"] = 1.5;
        self.config["double OuterRingMax"] = 1.99;
        self.config["double RingStep"] = 0.08;
        self.config["double RingParam_AreaOverPi"] = 0.6;
        self.config["double RingParam_ExcFracMax"] = 0.5;
        self.config["double RingParam_Thickness"] = 0.44;
        self.config["double RingParam_AlphaMax"] = 0.25;
        self.config["double StandardRingRadius"] = 0.7;
        self.config["double StandardRingThickness"] = 0.2;
        self.config["bool AverageOff"] = 1;
        self.config["bool CorrectZenith"] = 0;
        self.config["double OSRadius"] = 0.1;


        #  //Configure energy range
        self.config["double EMin"]  = 0;
        self.config["double EMax"]  = -1;

        self.config["bool FitAcceptance"]  = 1;
 
        #  //Configure flux products
        self.config["bool ProduceFluxProducts"]  = 1;
        self.config["double SpectralIndex"]  = 3.0;
        self.config["bool PointLikeFluxMaps"]  = 1;
        self.config["bool ProduceSurfaceBrightnessMap"]  = 0;
        self.config["bool ExposureMapsFromFits"]  = 0;

        ##############################################################################
        f = gtk.Frame("Flux Products ")
        h1.pack_start(f, expand=False, fill=False, padding=5)
        t4 = gtk.Table(3, 3)
        f.add(t4)

        self.ProduceFluxProducts_button = SignalWidgets.CheckButton("Produce Flux Products")
        self.ProduceFluxProducts_button.set_active(True)
        t4.attach(self.ProduceFluxProducts_button, 0, 1, 0, 1, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.ProduceFluxProducts_button.connect("toggled", self.changed_ProduceFluxProducts) 

        self.PointLikeFluxMaps_button = SignalWidgets.CheckButton("Produce PointLike Flux Maps")
        self.PointLikeFluxMaps_button.set_active(True)
        t4.attach(self.PointLikeFluxMaps_button, 0, 1, 1,2, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.PointLikeFluxMaps_button.connect("toggled", self.changed_PointLikeFluxMaps) 

        self.ProduceSurfaceBrightnessMap_button = SignalWidgets.CheckButton("Produce Surface Brightness Map")
        self.ProduceSurfaceBrightnessMap_button.set_active(False)
        t4.attach(self.ProduceSurfaceBrightnessMap_button, 0, 1, 2,3, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)
        self.ProduceSurfaceBrightnessMap_button.connect("toggled", self.changed_SurfaceBrightnessMap) 



        self.index_label = gtk.Label("Assumed Source Index :")
        t4.attach(self.index_label,0, 1,3,4, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1)

        self.index_entry = gtk.SpinButton(gtk.Adjustment(self.config["double SpectralIndex"], 0, 10, .1), .5, 1)
#        self.index_entry.set_numeric(True)
        self.index_entry.set_size_request(30, -1);
        t4.attach(self.index_entry,  1,2, 3,4, gtk.EXPAND | gtk.FILL, gtk.EXPAND | gtk.FILL, 5, 1) 



       ##############################################################################
#        f = gtk.Frame("Job ")
#        h1.pack_start(f, expand=False, fill=False, padding=5)
#        t7 = gtk.Table(3, 3)
#        f.add(t7)

#        self.ExposureMaps = 1
#        self.CreateLightCurve = 0
#        self.CreateFluxMaps = 0
#        self.butt1 = gtk.RadioButton(None, "ExposureMaps")
#        self.butt1.connect("toggled", self.jobs, "ExposureMaps")
##        t.attach(gtk.Label("EventsAndAcceptanceFromFile"), 0, 2, 0,1, gtk.EXPAND | gtk.FILL, 5, 5)
#        t7.attach(self.butt1, 2, 3, 0,1, gtk.EXPAND | gtk.FILL, 5, 5)

#        self.butt2 = gtk.RadioButton(self.butt1, "CreateLightCurve")
#        self.butt2.connect("toggled", self.jobs, "CreateLightCurve")
#        t7.attach(self.butt2, 2, 3, 1,2, gtk.EXPAND | gtk.FILL, 5, 5)


#        self.butt3 = gtk.RadioButton(self.butt1, "CreateFluxMaps")
#        self.butt3.connect("toggled", self.jobs, "CreateFluxMaps")
#        t7.attach(self.butt3, 2, 3, 2,3, gtk.EXPAND | gtk.FILL, 5, 5)
       ######################################### #####################################

        self.config["bool SaveResults"]  = 1;
        self.config["bool Verbose"]  = 1;
 


        self.config["char *Evtfile"] = '""';
        self.config["char *Accfile"] = '""';
        self.config["char *ExposureMapFile"] = '""';
        self.config["char *ExtendedExposureMapFile"] = '""';
        self.config["char *ResFileForFlux"] = '""';
        
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
        self.list.set_size_request(500, -1)
        self.list.set_column_width(0, 400)
        self.list.set_column_width(1, 100)
        swin = gtk.ScrolledWindow()
        swin.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        swin.set_size_request(-1, -1)
        swin.add(self.list)
        h1.pack_start(swin, expand=True, fill=True, padding=5)

#        tool = gtk.HBox()
#        h1.pack_start(tool, expand=False, fill=False, padding=5)
#        tool = gtk.HBox()
#        self.selectentry = gtk_supp.HistoryCombo(app_name="smash",
#                                                 history_key="SQLFilter",
#                                                 set_text=False)
#        selectapply = gtk.Button("", gtk.STOCK_APPLY)
#        selectapply.set_size_request(80, -1)
#        selectcancel = gtk.Button("", gtk.STOCK_CANCEL)
#        selectcancel.set_size_request(80, -1)
#        selectcancel.connect("clicked", lambda b, e=self.selectentry : e.child.set_text(""))
#        selectapply.connect("clicked", self.do_apply)
#        tool.pack_start(self.selectentry, True, True, padding=5)
#        tool.pack_start(selectapply, 0, 0, 1)
#        tool.pack_start(selectcancel, 0, 0, 1)
#        h1.pack_start(tool, expand=False, fill=False, padding=5)

#            
#    def do_apply(self, *args):
#        '''
#        Apply settings to widgets
#        '''
#        t = self.selectentry.child.get_text()
#        self.selectentry.append_history()
#        self.apply_selection(t)


    def changed_UseConfigTarget(self, *args):
        self.config["bool UseConfigTarget"] = 0;
        if self.UseConfigTarget_button.get_active() :
              self.config["bool UseConfigTarget"]  = 1;
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


    def changed_SurfaceBrightnessMap(self, *args):
        self.config["bool ProduceSurfaceBrightnessMap"]  = 0;
        if self.PointLikeFluxMaps_button.get_active() :
              self.config["bool ProduceSurfaceBrightnessMap"]  = 1;


    def changed_PointLikeFluxMaps(self, *args):
        self.config["bool PointLikeFluxMaps"]  = 0;
        if self.PointLikeFluxMaps_button.get_active() :
              self.config["bool PointLikeFluxMaps"]  = 1;


    def changed_ProduceFluxProducts(self, *args):
        self.config["bool ProduceFluxProducts"]  = 0;
        if self.ProduceFluxProducts_button.get_active() :
              self.config["bool ProduceFluxProducts"]  = 1;

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
            for f in glob.glob("*.root"):
                siz = os.path.getsize(f)
                self.list.append([f,str(int(siz/1e6))+' MB'], info=f)
        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(x)

        self.mainwindow.normal_cursor()


    def update(self, *args):
        self.config["double UserLambda"] = self.ra_entry.get_value();
        self.config["double UserBeta"] = self.dec_entry.get_value();
        self.config["double SpectralIndex"] = self.index_entry.get_value();
        self.config["char *AnalysisConfig"] = '"'+self.profile_combo.get_active_text()+'"'
        self.config["char *Regionsfile"] = '"'+self.Exclusion_file.get_filename()+'"'
               
        self.config["double ExtX"] = self.ExtX_entry.get_value();
        self.config["double ExtY"] = self.ExtY_entry.get_value();

#int RunNumberMin, int RunNumberMax

    def start_map(self, *args):
        self.update()



        if not self.mainwindow.validate_batchsystem():
            return
        try:
            rootfiles = self.list.get_selection_info()
            if len(rootfiles) == 0:
                ErrorDialogs.Error("No ROOT files selected! Please load the list and select the ROOT files you want to process.")
            for rf in rootfiles:
               filename = "CreateTimeMaps" 
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

               #TODO official location of the code
               batchfile.write_line("ln -s /sps/hess/users/lapp/dsanchez/HEGS/Scripts/RunWise_Generation/SurveySuiteFancyTiming.C . ")
               batchfile.write_line("ln -s /sps/hess/users/lapp/dsanchez/HEGS/Scripts/RunWise_Generation/SurveySuiteFancyTiming.hh . ")
               batchfile.write_line("ln -s /sps/hess/users/lapp/dsanchez/HEGS/Scripts/rootlogon.C . ")
               batchfile.write_line("ln -s /sps/hess/users/lapp/dsanchez/HEGS/Scripts/RunWise_Generation/RadialTables/ . ")
               batchfile.write_line("ln -s "+self.fin.get_current_folder()+"/"+rf+" . ")
               batchfile.write_line("cat > "+filename+".C  << EOF")
               batchfile.write_line("void "+filename+"()\n{\ngROOT->LoadMacro(\"SurveySuiteFancyTiming.C+\");\n")


               batchfile.write_line("SurveySuite::MapMaker *m = new SurveySuite::MapMaker();\n")

               self.config["std::string table_path"] = '"'+self.mainwindow.get_batchsystem().storage.find_directory(Storage.DIR_TABLES).get_base_directory()+'"'

               for key in self.config:
                   batchfile.write_line(key+" = "+str(self.config[key])+";")
               batchfile.write_line("m->StartConfigure();")
               batchfile.write_line("bool ConfigExclusionSuccess = m->ConfigureExclusions(ExclusionFromFits, ExclusionFromRegionFile);\n")
               batchfile.write_line("bool ConfigAccSuccess = m->ConfigureAcceptance(EventsAndAcceptanceFromFile, EventsAndAcceptanceFromRadial, EventsAndAcceptanceFromRadialFile, PsiCut, FitAcceptance, CorrectZenith, SelectAllRunsContributingToTheMap, ApplySafeThreshold, SafeThresholdFromAcceptance, SafeThresholdFromPsiCut, SafeThresholdFromPsiValue, SafeThresholdRatioParam);   \n")
               batchfile.write_line("bool ConfigMapSuccess = m->ConfigureMaps(UseConfigTarget, UseConfigMapParams, UserLambda, UserBeta, BinSize, ExtX, ExtY, OSRadius, EMin, EMax);\n")
               batchfile.write_line("bool ConfigRingSuccess = m->ConfigureRingMethod(AdaptFFT, AdaptCut_Alpha, ConstantArea, ConstantThickness, ConstantInnerRadius, SmoothRings, RingMinimalRadius, InnerRingMax, OuterRingMax, RingStep, RingParam_AreaOverPi, RingParam_ExcFracMax, RingParam_Thickness, RingParam_AlphaMax, StandardRingRadius, StandardRingThickness, AverageOff);\n")
               batchfile.write_line("bool ConfigFluxSuccess = m->ConfigureFluxProducts(ProduceFluxProducts, SpectralIndex, PointLikeFluxMaps, ProduceSurfaceBrightnessMap, ExposureMapsFromFits, AnalysisConfig);\n")
               batchfile.write_line("bool ConfigOutSuccess = m->ConfigureOutputs(SaveResults, Verbose);\n")
               batchfile.write_line("bool ConfigEndSuccess = m->EndConfigure(AnalysisConfig, Resfile, RunNumberMin, RunNumberMax, UseRunListToMatch, UseRunsToForbid, RunListToMatch, RunsToForbid, table_path);\n")
               batchfile.write_line("bool ConfigSuccess = 0;\n")
               batchfile.write_line("if(ConfigExclusionSuccess && ConfigAccSuccess && ConfigMapSuccess && ConfigRingSuccess && ConfigFluxSuccess && ConfigOutSuccess && ConfigEndSuccess){")
               batchfile.write_line("ConfigSuccess = 1;")
               batchfile.write_line("m->PrintConfigSummary();")
               batchfile.write_line("}\n")

#               batchfile.write_line("m->CreateExposureMaps_PerRun(AnalysisConfig, Resfile, Regionsfile, table_path, OutfilePrefix);\n")
#               batchfile.write_line("m->CreateLightCurve_PerRun(AnalysisConfig, Resfile, Evtfile, Accfile, ExposureMapFile, ExtendedExposureMapFile, Regionsfile, table_path, OutfilePrefix,UserLambda,UserBeta);\n")
               batchfile.write_line("m->CreateFluxMaps_PerRun(AnalysisConfig, Resfile, Evtfile, Accfile, ExposureMapFile,ExtendedExposureMapFile, Regionsfile,  table_path,  OutfilePrefix, \"withGui\");\n")

               batchfile.write_line("}")
               batchfile.write_line("EOF")
               root_opts =  "-q -b"
               batchfile.write_line("root %s '%s.C'" % (root_opts, filename))
               batchfile.save_result_file(".", "*.root" , ".", "")
               batchfile.save_result_file(".", "*.txt" , ".", "")
                
               r, jobname = batchfile.submit(terminal=False)
        except Exception, x:
            traceback.print_exc()
            ErrorDialogs.show_error(str(x))
        self.mainwindow.normal_cursor()


#    def Creatbatch(self):

    def acceptance(self, widget, data=None):
        self.config["bool EventsAndAcceptanceFromFile"] = 0
        self.config["bool EventsAndAcceptanceFromRadial"] = 0
        self.config["bool EventsAndAcceptanceFromRadialFile"] = 0
        
        self.config["bool "+data] = 1

#    def jobs(self,widget,data=None):
#        self.CreateExposureMaps = 0
#        self.CreateLightCurve = 0
#        self.CreateFluxMaps = 0
#        if data=="CreateLightCurve":
#           self.CreateLightCurve = 1
#        if data=="CreateExposureMaps":
#           self.CreateExposureMaps = 1
#        if data=="CreateFluxMaps":
#           self.CreateFluxMaps = 1

    def exclusion(self, widget, data=None):
        self.config["bool ExclusionFromFits"] = 0
        self.config["bool ExclusionFromRegionFile"] = 0
        
        self.config["bool "+data] = 1

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
        self.target_list.set_popdown_strings(targets)   

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
