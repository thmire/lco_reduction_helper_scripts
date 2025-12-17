#!/home/treynolds/mambaforge/envs/pypeit_dev/bin/python3

import os
import glob
from argparse import ArgumentParser
import subprocess
from astropy.io import fits
import json
import shutil
import logging


pypeit_path = "/home/treynolds/mambaforge/envs/pypeit_dev"

# Need to allow this to be changed to make this script portable.
path = "/home/treynolds/data/LCO_spectra_reduction/scripts"

parser = ArgumentParser(
    description="Reduces all the data, after downloading and sorting"
)

parser.add_argument(
    "--input_folder",
    default=None,
    help="Path where the data was downloaded and sorted",
)

parser.add_argument(
    "--sensfunc",
    default=True,
    help="Make a sensfunc or use the Master sensfunc",
)

parser.add_argument(
    "--retry",
    default=None,
    help="If True, retry previously failed reductions",
)

parser.add_argument(
    "--extraction_only",
    default=None,
    help="If True, redo all reductions, but dont remake calibs.",
)

args = parser.parse_args()
input_folder = args.input_folder
sensfunc = args.sensfunc
retry = args.retry
extraction_only = args.extraction_only

HERE = os.getcwd()
if os.path.isabs(input_folder):
    input_folder_abs_path = input_folder
else:
    input_folder_abs_path = os.path.join(HERE, input_folder)
all_folders = glob.glob(input_folder_abs_path+"/*")

# Should do some check in the diags that the folder is ready to go.
with open(input_folder_abs_path + "/data_summary_all.json", 'r') as f:
    data_summary_all = json.load(f)

# Set up the logger
# should probably do more logging, and produce a proper log file.
logger = logging.getLogger(__name__)
logging.basicConfig(filename=input_folder_abs_path + "/logging_reduction.log", encoding='utf-8', level=logging.DEBUG)

#print(data_summary_all)
reduction_status_dictionary = {}
for folder in all_folders :
    date = folder.split("/")[-1]
    # Currently, I am excluding any random files in the folder by hand
    # do something smarter than this.
    if date == "data_summary_all.json":
        continue
    if date == "logging_reduction.log":
        continue
    if date == "reduction_diags_all.json":
        continue
    try :
        with open(folder + '/reduction_diags_'+ date + '.json', 'r') as f: 
            reduction_status_dictionary[date] = json.load(f)
        if reduction_status_dictionary[date]["Reduction_complete"] == True:
            if not extraction_only :
                print(date + " reduction already complete",
                  "skipping...")
                continue                                                             
        else :
            print(date + " reduction already attempted, but failed.\n")
            if not retry :
                print("skipping...")
                continue
            else :
                print("retrying science reduction...")
    except FileNotFoundError :
        reduction_status_dictionary[date] = {"Reduction_complete":False}
        with open(folder + '/reduction_diags_' + date + ".json", 'w') as f:
            json.dump(reduction_status_dictionary[date], f, indent=4)    
    # Is the folder ready?
    if not data_summary_all[date]["File_complete"] == True :
        with open(folder + '/data_diags_'+ date + '.json', 'r') as f:
            data_diags_dictionary = json.load(f)
            if data_diags_dictionary["science"] and data_diags_dictionary["sci_arc"] and data_diags_dictionary["sci_flat"] and sensfunc == "Master":
                print("We have all the science files, and you requested the Master sensfunc. Continuing...")
            else:
                print(date + " is not ready, skipping...")
                continue
    else :
        print(date + " has all the required files, reducing...")

    telescope = data_summary_all[date]["telescope"]
    if telescope == "LCO_north":
        reduction_status_dictionary[date]["telescope"] = telescope
        pypeit_spectrograph = "lco_floyds_north"
    elif telescope == "LCO_south" :
        reduction_status_dictionary[date]["telescope"] = telescope
        pypeit_spectrograph = "lco_floyds_south"

## I used this code block to re-reduce with a change to extraction. Probably obselete 
##    skip=None    
##    reduction_pypeit_file_sci_list = glob.glob(folder + "/sci/*"+ pypeit_spectrograph +"*/*.pypeit")
##    with open(reduction_pypeit_file_sci_list[0],"r") as file:
##        pypeit_file_data = file.readlines()
##        for i in range(0,len(pypeit_file_data)) :
##            if i == 6 :
##                if pypeit_file_data[i] == "[reduce]\n":
##                    print(date + " re-reduction already complete",
##                              "skipping...")
##                    skip = True
##                    break
##        
##    if skip :
##        skip = None
##        continue
    
    # Now reduce the data.
    # First reduce the std data
    # Probably want the option to completely turn off the std reduction.
    # Pypeit needs us to be working in the right directory.
    os.chdir(folder+"/std")

    if sensfunc == "Master":
        if pypeit_spectrograph == "lco_floyds_south":
            print("Master sensfunc requested, using archive sensfunc for LCO south...")
            sensfunc_file = path + "/archive_sensfunc_LCO_south.fits"
            reduction_status_dictionary[date]["sensfunc"] = "archive_sensfunc_LCO_south.fits"
        elif pypeit_spectrograph == "lco_floyds_north":
            print("Master sensfunc requested, using archive sensfunc for LCO north...")
            sensfunc_file = path + "/archive_sensfunc_LCO_north.fits"
            reduction_status_dictionary[date]["sensfunc"] = "archive_sensfunc_LCO_north.fits"
    else:
        # check if sensfunc already exists:
        sensfunc_file = folder + "/std/" + "sensfunc_"+date+".fits"
        if not os.path.isfile(sensfunc_file) and sensfunc :
            print("no sensfunc ready")
            if not retry :
                # Create the pypeit reduction file
                subprocess.call([pypeit_path + "/bin/pypeit_setup",
                             "-s",pypeit_spectrograph,
                             "-r", folder + "/std/raw_files/",
                             "-c" "all"]) 
                # Should check the setup is good, make any edits we want.
                # Function goes here?
                
                reduction_pypeit_file_std = glob.glob(folder + "/std/*" + pypeit_spectrograph + "*/*.pypeit")
                # Always redo the calibs and overwrite old reductions?
                print("running std reduction now...")
                subprocess.call([pypeit_path + "/bin/run_pypeit",
                             reduction_pypeit_file_std[0],
                             "-o","-m",
                             "-v","0", # turn the verbosity to none.
                             #"-c" "all"
                                 ])

                reduced_std = glob.glob(folder + "/std/Science/spec1d*.fits")
                # Produce a sensitivity function
                print("Producing sensitivity function...")
                # Copy a parameter file with preset good parameters.

                if sensfunc == True :
                    #try :
                    subprocess.call([pypeit_path + "/bin/pypeit_sensfunc",
                             reduced_std[0],
                             #"--algorithm","UVIS", # use the UVIS algorithm
                             "-v","0", # turn the verbosity to none.
                             "--par_outfile","sensfunc_pars.par", # Write the parameters used to this file
                             "-o","sensfunc_"+date+".fits", # Write the sensfunc to this file
                             "-s", path + "/sensfunc_good_params.par" #link to best param file
                             #"-c" "all"
                                 ])
                    sensfunc_file = folder + "/std/" + "sensfunc_"+date+".fits"
                    if not os.path.exists(sensfunc_file):
                        if pypeit_spectrograph == "lco_floyds_south":
                            print("Sensfunc failed, using archive sensfunc for LCO south...")
                            sensfunc_file = path + "/archive_sensfunc_LCO_south.fits"
                            reduction_status_dictionary[date]["sensfunc"] = "archive_sensfunc_LCO_south.fits"
                        elif pypeit_spectrograph == "lco_floyds_north":
                            print("Sensfunc failed, using archive sensfunc for LCO north...")
                            sensfunc_file = path + "/archive_sensfunc_LCO_north.fits"
                            reduction_status_dictionary[date]["sensfunc"] = "archive_sensfunc_LCO_north.fits"
                    else :
                        reduction_status_dictionary[date]["sensfunc"] = "sensfunc_"+date+".fits"
                else :
                    print("Sensfunc exists and you don't want to retry, skipping std reduction")
            

                
    #print(sensfunc_file)
    #break
    # Reduce the science
    # Pypeit needs us to be working in the right directory.
    os.chdir(folder+"/sci")
    # Create the pypeit reduction file
    subprocess.call([pypeit_path + "/bin/pypeit_setup",
                 "-s",pypeit_spectrograph,
                 "-r", folder + "/sci/raw_files/",
                 "-c" "all",
                 "-v","0", # turn verbosity down
                     ])
    # Should check the setup is good, make any edits we want.
    # Function goes here?

    # Now reduce the science
    reduction_pypeit_file_sci_list = glob.glob(folder + "/sci/*"+ pypeit_spectrograph +"*/*.pypeit")
    # assuming here that there is one skyframe problem file and one not
    # If we get a nice set of skyframes, this code block will no longer be needed.
    if len(reduction_pypeit_file_sci_list) != 1 and pypeit_spectrograph == "lco_floyds_north":
        print("We have multiple setups, need to consolidate...")
        with open(reduction_pypeit_file_sci_list[0],"r") as file:
            reduction_pypeit_file_sci_0 = file.readlines()
        with open(reduction_pypeit_file_sci_list[1],"r") as file:
            reduction_pypeit_file_sci_1 = file.readlines()
        for line in reduction_pypeit_file_sci_0:
            if "trace" in line:
                print(reduction_pypeit_file_sci_list[0], " has the skyframe")
                skyframe_file = reduction_pypeit_file_sci_list[0]
                sci_file = reduction_pypeit_file_sci_list[1]
                break
        if not sci_file:
            skyframe_file = reduction_pypeit_file_sci_list[1]
            sci_file = reduction_pypeit_file_sci_list[0]
        with open(sci_file,"r") as file:
            reduction_pypeit_file_sci = file.readlines()
        with open(skyframe_file,"r") as file:
            reduction_pypeit_file_skyframe = file.readlines()
        reduction_pypeit_file_sci[-2] = reduction_pypeit_file_skyframe[-3][:-3] +"1\n"
        reduction_pypeit_file_sci.append("data end")
   #     if pypeit_spectrograph == "lco_floyds_south":
   #         with open(folder + "/sci/" + pypeit_spectrograph + "_A/lco_floyds_south_A.pypeit","w") as file:
   #             print(reduction_pypeit_file_sci)
   #             file.writelines(reduction_pypeit_file_sci)          
   #         reduction_pypeit_file_sci[0] = folder + "/sci/" + pypeit_spectrograph + "_A/lco_floyds_south_A.pypeit"
        with open(folder + "/sci/" + pypeit_spectrograph + "_A/lco_floyds_north_A.pypeit","w") as file:
            file.writelines(reduction_pypeit_file_sci)          
        reduction_pypeit_file_sci[0] = folder + "/sci/" + pypeit_spectrograph + "_A/lco_floyds_north_A.pypeit"

    # If you need to edit the .pypeit file, do so here.
    raw_files = glob.glob(folder + "/sci/raw_files/*")
    for key in data_summary_all[date]["files"].keys():
        if data_summary_all[date]["files"][key]["filetype"] == "SPECTRUM":
            if data_summary_all[date]["files"][key]["target"]=="science" :
                slit_width = data_summary_all[date]["files"][key]["slitwidth"]
    print("slit_width: ",slit_width)            
    if folder + "/sci/raw_files/master_arc_LCO_north_20240701.fits" or folder + "/sci/raw_files/master_arc_LCO_south_20240706.fits" in raw_files :
        if slit_width != 1.6:
            print("Adding the master arc to the _A pypeit_file")
            if len(reduction_pypeit_file_sci_list) != 1:
                for pypeit_file in reduction_pypeit_file_sci_list :
                    with open(pypeit_file,"r") as file:
                        reduction_pypeit_file_sci_0 = file.readlines()
                        for line in reduction_pypeit_file_sci_0:
                            if "master_arc_" in line:
                                print(pypeit_file, " has the master_arc")
                                master_arc_file = pypeit_file
                                arc_line_to_add = line
                                break
            # add the master arc to the _A

            # edit the master arc line:
            if pypeit_spectrograph == "lco_floyds_north":
                arc_line_to_add = arc_line_to_add[:-3] + "1\n"
            elif pypeit_spectrograph == "lco_floyds_south":
                arc_line_to_add = arc_line_to_add[:-3] + "0\n"
                
            
            with open(reduction_pypeit_file_sci_list[0],"r") as file:
                pypeit_file_data = file.readlines()
            pypeit_file_data_new = []
            print(len(pypeit_file_data))
            for i in range(0,len(pypeit_file_data)) :
                if i == 20 :
                    #print("adding the line",arc_line_to_add)
                    pypeit_file_data_new.append(arc_line_to_add)
                    pypeit_file_data_new.append(pypeit_file_data[i])     
                else :
                    pypeit_file_data_new.append(pypeit_file_data[i])
            #print(pypeit_file_data_new)
            with open(reduction_pypeit_file_sci_list[0],"w") as file:
                file.writelines(pypeit_file_data_new)
        else :
            print('Data taken with 1.6" slit, can happily use the master arc with no edits')

    # Add the no-extract-mask setting
    with open(reduction_pypeit_file_sci_list[0],"r") as file:
        pypeit_file_data = file.readlines()
    pypeit_file_data_new = []
    for i in range(0,len(pypeit_file_data)) :
        if i == 6 :
            #print("adding the line",arc_line_to_add)
            pypeit_file_data_new.append("[reduce]\n    [[extraction]]\n        use_2dmodel_mask = False\n")
            pypeit_file_data_new.append(pypeit_file_data[i])     
        else:
            pypeit_file_data_new.append(pypeit_file_data[i])
            #print(pypeit_file_data_new)
    with open(reduction_pypeit_file_sci_list[0],"w") as file:
        file.writelines(pypeit_file_data_new)

              
    # If the arc lamp is 6.0", need to change the fitting order.
    # This should not be needed any more?
    # Added a line to stop these being added when we are using the master arc
    # Should combine all the if statements below!
    
    #print(data_summary_all[date][key])
##    for key in data_summary_all[date]["files"].keys():
##        #print(data_summary_all[date][key])
##        if data_summary_all[date]["files"][key]["filetype"] == "SPECTRUM":
##            if data_summary_all[date]["files"][key]["target"]=="science" :
##                if data_summary_all[date]["files"][key]["slitwidth"] == 6.0:
##                        if not folder + "/sci/raw_files/master_arc_LCO_north_20240701.fits" or folder + "/sci/raw_files/master_arc_LCO_south_20240706.fits" in raw_files: 
##                            with open(reduction_pypeit_file_sci_list[0],"r") as file:
##                                pypeit_file_data = file.readlines()
##                            pypeit_file_data_new = []
##                            for i in range(0,len(pypeit_file_data)) :
##                                if i == 6:
##                                    print("adding the line")
##                                    pypeit_file_data_new.append("[calibrations]\n [[wavelengths]]\n  n_first = 2\n  ech_nspec_coeff = 3\n")
##                                    pypeit_file_data_new.append(pypeit_file_data[i])     
##                                else :
##                                    pypeit_file_data_new.append(pypeit_file_data[i])
##                            with open(reduction_pypeit_file_sci_list[0],"w") as file:
##                                file.writelines(pypeit_file_data_new)
##                            break
                
        
    # Always redo the calibs and overwrite old reductions?
    print("running sci reduction now...")
    print(reduction_pypeit_file_sci_list[0])
    if extraction_only :
        print("Only redo extraction, using old calibrations")
        subprocess.call([pypeit_path + "/bin/run_pypeit",
                 reduction_pypeit_file_sci_list[0],
                 "-o",
                 "-v","0", # turn the verbosity to none.
                     ])
    else :
        subprocess.call([pypeit_path + "/bin/run_pypeit",
                 reduction_pypeit_file_sci_list[0],
                 "-o","-m",
                 "-v","0", # turn the verbosity to none.
                     ])
    # copy the sensfunc to the sci folder
    shutil.copy(sensfunc_file,
                folder + "/sci/Science/" + sensfunc_file.split("/")[-1])
    
    print("Setup flux calibration and coadds...")
    subprocess.call([pypeit_path + "/bin/pypeit_flux_setup",
                 folder + "/sci/Science"  ,
                 "--name",date,
                     ])
    # Some edits to the flux calibration file.
    with open(folder + "/sci/" + date + ".flux","r") as file:
        flux_calib_file_data = file.readlines()
    # As we use the UVIS algo, need to update the .flux file to extinct_correct
    flux_calib_file_data[5] = "  extinct_correct = True\n"
    # We also need to add the sensitivity function to the file.
    flux_calib_file_data[12] = flux_calib_file_data[12][:-2] + sensfunc_file.split("/")[-1] + "\n"
    with open(folder + "/sci/" + date + ".flux","w") as file:
        file.writelines(flux_calib_file_data)
    
    print("Apply the sensitivity function...")
    subprocess.call([pypeit_path + "/bin/pypeit_flux_calib",
                 folder + "/sci/" + date + ".flux",
                 "-v","0", # turn the verbosity to none.
                     ])
    
    # Need to make a number of edits to the coadd file:
    # Could be a big problem here, if there are multiple science files that we don't want to coadd.
    with open(folder + "/sci/" + date + ".coadd1d","r") as file:
        flux_coadd1d_file_data = file.readlines()
    flux_coadd1d_file_data_new = []
    # Provide a name for the final file
    for i in range(0,len(flux_coadd1d_file_data)-1) :
        if i < 5 :
            flux_coadd1d_file_data_new.append(flux_coadd1d_file_data[i])
        elif i == 5 :
            # Provide a name for the final file
            flux_coadd1d_file_data_new.append("  coaddfile = "+ date +"_coadded.fits  # Please set your output file name\n")
        elif i == 6 :
            flux_coadd1d_file_data_new.append("  wave_method = velocity  # creates a uniformly space grid in log10(lambda)\n")
        elif i == 7 :
            flux_coadd1d_file_data_new.append("  dv = 180  # size for velocity grid\n")
            flux_coadd1d_file_data_new.append(flux_coadd1d_file_data[i])
        elif i == 15 :
            print(flux_coadd1d_file_data[i])
            line = " | ".join(flux_coadd1d_file_data[i].split("|")[0:2]) + " | " +\
                   sensfunc_file.split("/")[-1] + " | " + flux_coadd1d_file_data[i].split("|")[3] + "\n"
            print("Tricky line: ",line)
            flux_coadd1d_file_data_new.append(line)
        else :
            flux_coadd1d_file_data_new.append(flux_coadd1d_file_data[i])
    with open(folder + "/sci/" + date + ".coadd1d","w") as file:
        file.writelines(flux_coadd1d_file_data_new)
            
    print("Coadding the spectra...")
    
    subprocess.call([pypeit_path + "/bin/pypeit_coadd_1dspec",
                 folder + "/sci/" + date + ".coadd1d",
                 "-v","0", # turn the verbosity to none.
                     ])
    
    if os.path.isfile(folder + "/sci/" + date + "_coadded.fits"):
        print(date + " should be done! On to the next...")
        reduction_status_dictionary[date]["Reduction_complete"] = True
    else :
        print("Something went wrong, no coadded file found")
        reduction_status_dictionary[date]["Reduction_complete"] = False
    with open(folder + '/reduction_diags_' + date + ".json", 'w') as f:
        json.dump(reduction_status_dictionary[date], f, indent=4)
    with open(input_folder_abs_path + "/reduction_diags_all.json","w") as f:
        json.dump(reduction_status_dictionary, f, indent=4)
        


