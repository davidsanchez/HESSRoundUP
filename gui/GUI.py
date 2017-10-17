#!/usr/bin/python
# author David Sanchez david.sanchez@lapp.in2p3.fr
'''
Main GUI for the HEGS module
'''
import sys
import os
import gtk, gobject

try:
    sys.path.append(os.environ["HESSUSER"])
except KeyError:
    pass
sys.path.append(os.environ["HESSROOT"])

from display.python import MainGUI, Storage
from MapsWindow import MapsWindow
from TimeMapsWindow import TimeMapsWindow
from ToolsWindow import ToolsWindow
from CalibCheckWindow import CalibCheckWindow
from AnalyseWindow import AnalyseWindow
#from Environment import EnvironmentGUI



def main():
    '''
    Main method which creates the GUI
    '''
    gtk.rc_parse(MainGUI.RC_FILE)
    directories = [Storage.DIR_CALIBRATION,
                   Storage.DIR_DST,
                   Storage.DIR_TABLES,
                   Storage.DIR_OUTPUT_TABLES,
                   Storage.DIR_RESULTS]
    window = MainGUI.MainGUI("HEGS_Manager",
                             "HEGSGUI",
                             directories,
                             "HEGS Manager",
                             xrootd=False)

    window.add_page("Analysis", AnalyseWindow(window))
    window.add_page("Maps Generation", MapsWindow(window))
#    window.add_page("Time Maps Generation", TimeMapsWindow(window))
    window.add_page("Tools", ToolsWindow(window))
    window.add_page("Calibration", CalibCheckWindow(window))
#    window.init(envtab=EnvironmentGUI(window))
    window.init()
    gobject.threads_init()
    gtk.main()




if __name__ == '__main__': main()
