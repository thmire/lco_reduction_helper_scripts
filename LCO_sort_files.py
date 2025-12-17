#!/home/treynolds/mambaforge/envs/astro/bin/python3

import os
import glob
from argparse import ArgumentParser
import subprocess
from astropy.io import fits
import json
import shutil
import numpy as np

def arc_check(arc_fits_file):
    """
    Returns True if the arc is messed up (see 20220714 std for example)
    Takes a median (for cosmics) of some rows that should be low counts.
    """
    with fits.open(arc_fits_file) as arc_data:
        if np.max(np.median(arc_data[0].data[290:310],axis=0)) > 700:
            print("Bad arc, pixel value max is:,", np.max(arc_data[0].data[300]))
            return True
        else :
            return False


parser = ArgumentParser(
    description="Sorts the data from LCO for reduction"
)

parser.add_argument(
    "--input_folder",
    default=None,
    help="Path where the data was downloaded",
)

parser.add_argument(
    "--science_target",
    default=None,
    help="Name of the science target. Currently unused! (2025-04-29)",
)


parser.add_argument(
    "--skyframe_path",
    default="/home/treynolds/data/LCO_spectra_reduction/SKYFRAMES",
    help="Path to location where the skyframes are stored.",
)

parser.add_argument(
    "--use_master_arc",
    default=True,
    help="Flag to use the master arcs for all wavelength calibrations. Default is to use masters.",

    )


args = parser.parse_args()
input_folder = args.input_folder
skyframe_path = args.skyframe_path
use_master_arc = args.use_master_arc

# Update later for different science targets
science_target = "Mrk 590"

# this should not be hard-coded, change later.
master_arc_path = "/home/treynolds/data/LCO_spectra_reduction/MASTER_ARC"


HERE = os.getcwd()
if os.path.isabs(input_folder):
    input_folder_abs_path = input_folder
else:
    input_folder_abs_path = os.path.join(HERE, input_folder)

index_letter_dictionary = {
    "e":"SPECTRUM",
    "a":"ARC",
    "w":"FLAT",
    }

skyframes_dict= {"LCO_south":{1.2:"SKYFLAT_LCO_south_1.2arcsec_Jan2025.fits",
                              1.6:"SKYFLAT_LCO_south_1.6arcsec_Jan2025.fits",
                              2.0:"SKYFLAT_LCO_south_2.0arcsec_Jan2025.fits",
                              6.0:"SKYFLAT_LCO_south_6.0arcsec_Jan2025.fits"},
                 "LCO_north":{}}

# Get all the folders
all_folders = glob.glob(input_folder_abs_path+"/*")
all_files_dictionary = {}

incomplete_folders = []
    
# Iterate through the folders
# Want to:
# 1) Produce a dictionary which summarises all the data
# 2) create sci and std folders, and mv the data into the correct places for the
#    reduction

# Note that this is not set up to work on already sorted fodlers yet.


for folder in all_folders :
    # Should check for an output file from a previous run here, and skip if
    # everything is already sorted.
    break_code = 0
    try :
        with open(folder + '/data_diags.json', 'r') as f: 
            my_dict = json.load(f)
        print("Data already downloaded and sorted by this script. ",
              "Check the jsons for info. Skipping...")
        
    except FileNotFoundError :
        pass
    
    #populate the data dictionary
    try :
        os.mkdir(folder+"/sci")
        os.mkdir(folder+"/sci/raw_files")
        os.mkdir(folder+"/std")
        os.mkdir(folder+"/std/raw_files")
    except FileExistsError:
        pass
    date = folder.split("/")[-1]
    all_files_dictionary[date]={}
    file_list = glob.glob(folder + "/*.fits*")
    print(folder)
    telescope_code = file_list[0].split("/")[-1][0:3]
    if telescope_code == "coj":
        telescope = "LCO_south"
    elif telescope_code == "ogg":
        telescope = "LCO_north"
    else :
        print("Strange fits file, I don't know what telescope it comes from!")
    all_files_dictionary[date]["telescope"] = telescope
    all_files_dictionary[date]["files"]={}
    for file in file_list :
        # unpack the files
        if ".fits.fz" in file:
            #print("Unpacking file: ",file)
            subprocess.call(["funpack", "-F",file])
            new_filename = file[:-3]
        else :
            new_filename = file
        new_filename_short = new_filename.split("/")[-1]
        all_files_dictionary[date]["files"][new_filename_short]={}
        
        # check there are not two different telescopes in this folder
        telescope_code = new_filename_short[0:3]
        if telescope_code == "coj":
            if all_files_dictionary[date]["telescope"] == "LCO_north":
                print("There are files from both north and south here!",folder)
                break_code = 1
                break
        elif telescope_code == "ogg":
            if all_files_dictionary[date]["telescope"] == "LCO_south":
                print("There are files from both north and south here!",folder)
                break_code = 1
                break
            
        # categorise the files into a type: arc,flat or spectrum
        index_letter = file.split("-")[-1][0]
        filetype = index_letter_dictionary[index_letter]
        all_files_dictionary[date]["files"][new_filename_short]["filetype"] = filetype
        if filetype == "ARC" :
            bad_arc_check = arc_check(new_filename)
            if bad_arc_check == True :
                all_files_dictionary[date]["files"][new_filename_short]["arc_status"] = "bad"
            else :
                all_files_dictionary[date]["files"][new_filename_short]["arc_status"] = "good"
            
        with fits.open(new_filename) as hdul:
            header = hdul[0].header
            date_obs = header["DATE"]
            slitwidth = header["APERWID"]    
            if header["OBJECT"] == science_target:
                filetype_2 = "science"
            else :
                filetype_2 = "standard"
        all_files_dictionary[date]["files"][new_filename_short]["slitwidth"] = slitwidth
        all_files_dictionary[date]["files"][new_filename_short]["target"] = filetype_2
        all_files_dictionary[date]["files"][new_filename_short]["date"] = date_obs
    if break_code == 1:
        print("Cancelling: ", date)
        diag_dict = {"Failure":"There are files from multiple telescopes here"}
        with open(folder + '/data_diags_' + date + '.json', 'w') as f:
            json.dump(diag_dict, f, indent=4)  # indent for readability
        break
        
    # Now sort the files into reduction ready folders. Should also check here
    # that you have all the necessary files.
    #print(all_files_dictionary)       

    diag_dict = {"science":False, "standard":False,
                 "std_arc":False, "sci_arc":False,
                 "sci_flat":False, "std_flat":False}

    # Do we have: science frames, with arc and flat, and a std, with arc and flat.
    # Also, shuffle the files into the correct folders.
    for key in all_files_dictionary[date]["files"].keys():
        file_dict = all_files_dictionary[date]["files"][key]
        if file_dict["target"] == "science":
            if file_dict["filetype"] == "FLAT":
                diag_dict["sci_flat"] = True
            if file_dict["filetype"] == "ARC":
                if file_dict["arc_status"] == "bad" :
                    diag_dict["sci_arc"] = False
                else :
                    diag_dict["sci_arc"] = True
            elif file_dict["filetype"] == "SPECTRUM":
                diag_dict["science"] = True
            os.rename(folder+"/"+key, folder + "/sci/raw_files/" + key)
        elif file_dict["target"] == "standard":
            if file_dict["filetype"] == "FLAT":
                diag_dict["std_flat"] = True
            if file_dict["filetype"] == "ARC":
                if file_dict["arc_status"] == "bad" :
                    diag_dict["std_arc"] = False
                else :
                    diag_dict["std_arc"] = True
            elif file_dict["filetype"] == "SPECTRUM":
                diag_dict["standard"] = True
            os.rename(folder+"/"+key, folder + "/std/raw_files/" + key)
    files_complete = True        
    for key in diag_dict.keys():
        if diag_dict[key] == False:
            print("Missing: ", key)
            files_complete = False
    if files_complete == False :
        incomplete_folders.append(date)
    all_files_dictionary[date]["File_complete"] = files_complete
            
    all_files_dictionary[date]["skyframes"]={}

    # Put the skyframes in the folders ready to be used, record in dictionary
    # currently only have one LCO_north skyframe, update later
    if all_files_dictionary[date]["telescope"] == "LCO_north":
        shutil.copy(skyframe_path+"/LCO_north/LCO_floyds_north_skyframe.fits",
                    folder + "/sci/raw_files/LCO_floyds_north_skyframe.fits")
        all_files_dictionary[date]["skyframes"]["science"]="LCO_floyds_north_skyframe.fits"
        shutil.copy(skyframe_path+"/LCO_north/LCO_floyds_north_skyframe.fits",
                    folder + "/std/raw_files/LCO_floyds_north_skyframe.fits")
        all_files_dictionary[date]["skyframes"]["standard"]="LCO_floyds_north_skyframe.fits"

    # Have a nice new recent set of skyframes for the south
    elif all_files_dictionary[date]["telescope"] == "LCO_south":
        for key in all_files_dictionary[date]["files"].keys():

            file_dict = all_files_dictionary[date]["files"][key]
            if file_dict["filetype"] == "SPECTRUM" and file_dict["target"] == "science":
                slitwidth = file_dict["slitwidth"]
                shutil.copy(skyframe_path+"/LCO_south/" + skyframes_dict["LCO_south"][slitwidth],
                            folder + "/sci/raw_files/" + skyframes_dict["LCO_south"][slitwidth])
                all_files_dictionary[date]["skyframes"]["science"]=skyframes_dict["LCO_south"][slitwidth]
            if file_dict["filetype"] == "SPECTRUM" and file_dict["target"] == "standard":
                slitwidth = file_dict["slitwidth"]
                shutil.copy(skyframe_path+"/LCO_south/" + skyframes_dict["LCO_south"][slitwidth],
                            folder + "/std/raw_files/" + skyframes_dict["LCO_south"][slitwidth])                        
                all_files_dictionary[date]["skyframes"]["standard"]=skyframes_dict["LCO_south"][slitwidth]


    if use_master_arc :
        # Swap to use the master arc for science reductions
        if all_files_dictionary[date]["telescope"] == "LCO_north":
            shutil.copy(master_arc_path+"/LCO_north/master_arc_LCO_north_20240701.fits",
                    folder + "/sci/raw_files/master_arc_LCO_north_20240701.fits")
        elif all_files_dictionary[date]["telescope"] == "LCO_south":
            shutil.copy(master_arc_path+"/LCO_south/master_arc_LCO_south_20240109.fits",
                    folder + "/sci/raw_files/master_arc_LCO_south_20240706.fits")
            
        for file in all_files_dictionary[date]["files"].keys():
            if all_files_dictionary[date]["files"][file]["filetype"] == "ARC" and all_files_dictionary[date]["files"][file]["target"] == "science":
                shutil.move(folder+"/sci/raw_files/" + file,
                    folder + "/sci/" + file)                
                
    if files_complete == True:
        print(folder, " is complete and ready for reduction")

    # Should write an output file with the dictionary here.

    # Write jsons with both dictionaries:
    with open(folder + "/data_summary_" + date + ".json", 'w') as f:
        json.dump(all_files_dictionary[date], f, indent=4)    
    with open(folder + '/data_diags_' + date + '.json', 'w') as f:
        json.dump(diag_dict, f, indent=4)  

with open(input_folder_abs_path + "/data_summary_all.json", 'w') as f:
    json.dump(all_files_dictionary, f, indent=4)  # indent for readability    
            
print("Folders with missing files: ", incomplete_folders)

