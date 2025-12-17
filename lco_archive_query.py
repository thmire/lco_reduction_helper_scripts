#!/home/treynolds/mambaforge/envs/astro/bin/python3

# Script to download lco files from the archive

# Script is based on:
# https://github.com/cylammarco/FLOYDS_pipeline/blob/main/query_lco_archive.py
# That script is based on:
# https://github.com/LCOGT/lcogtsnpipe/blob/master/trunk/bin/LCOGTingest.py
# Also: https://lco.global/documentation/archive-documentation/

from argparse import ArgumentParser
from collections import Counter
from datetime import timedelta, datetime
import getpass
import glob
import os
import sys
import requests
import re

from astropy.io import fits
import numpy as np
import ruamel.yaml

# Move functions elsewhere?
def get_metadata(authtoken={}, limit=None, logger=None, **kwargs):
    """
    Get the list of files meeting criteria in kwargs

    """

    url = "https://archive-api.lco.global/frames/?" + "&".join(
        [
            key + "=" + str(val)
            for key, val in kwargs.items()
            if val is not None
        ]
    )
    url = url.replace("False", "false")
    url = url.replace("True", "true")
    if logger is not None:
        logger.info(url)
    else:
        print(url)

    response = requests.get(url, headers=authtoken, stream=True).json()
    frames = response["results"]
    while response["next"] and (limit is None or len(frames) < limit):
        if logger is not None:
            logger.info(response["next"])
        else:
            print(response["next"])
        response = requests.get(
            response["next"], headers=authtoken, stream=True
        ).json()
        frames += response["results"]
    return frames[:limit]

def download_frame(frame, base_directory, no_date=False, overwrite=False):
    """
    Download a single image from the LCOGT archive

    Parameters
    ==========
    frame: dict
        Dictionary containing the response from requests.
    base_directory: str
        The relative/full path of where the files are to be stored.
    no_date: boolean
        Set to True to store everything in base_directory, else, the files
        will be stored in folders named by the date of the night of the
        observation.

    """

    filename = frame["filename"]
    dayobs = re.search("(20\d\d)(0\d|1[0-2])([0-2]\d|3[01])", filename).group()
    if no_date:
        filepath = os.path.join(base_directory)
    else:
        filepath = os.path.join(base_directory, dayobs)

    if not os.path.isdir(filepath):
        os.makedirs(filepath)

    filename = frame["filename"]

    if not os.path.exists(os.path.join(filepath, filename)) or overwrite:
        with open(os.path.join(filepath, filename), "wb") as f:
            f.write(requests.get(frame["url"]).content)

    return filepath, filename, dayobs

parser = ArgumentParser(
    description="Downloads and Reduce data from archive.lco.global."
)

##parser.add_argument(
##    "--target_name",
##    default=None,
##    help="The target name to be queried on the TNS.",
##)

parser.add_argument(
    "--ra",
    default=None,
    help="Right Ascension in decimal. Only used if target_name is None.",
)

parser.add_argument(
    "--dec",
    default=None,
    help="Declination in decimal. Only used if target_name is None.",
)

parser.add_argument(
    "--input_folder",
    default=None,
    help="Path to install the the data",
)

parser.add_argument(
    "--telescope",
    default=None,
    help="Choose either the 'north' or 'south' telescope.",
)

parser.add_argument(
    "--lco_token",
    #default=None,
    help=(
        "LCO token. Only used if --login is None. Will ask for one if"
        "needed and nothing is provided."
    ),
)

parser.add_argument(
    "--date_start",
    default="1900-01-01",
    help="The date of the beginning of the night of the observation.",
)

parser.add_argument(
    "--date_end",
    default="2100-12-31",
    help="The date of the beginning of the night of the observation.",
)

parser.add_argument(
    "--propID",
    default=None,
    help="The proposal ID.",
)

#parser.add_argument(
#    "--propID",
#    default=None,
#    help="target_name",
#)

args = parser.parse_args()

input_folder = args.input_folder
#target_name = args.target_name
#ra = float(args.ra)
#dec = float(args.dec)
lco_token = args.lco_token
telescope = args.telescope
authtoken = {"Authorization": "Token {}".format(lco_token)}
propID = args.propID

# Get the absolute path to the input folder
HERE = os.getcwd()
if os.path.isabs(input_folder):
    input_folder_abs_path = input_folder
else:
    input_folder_abs_path = os.path.join(HERE, input_folder)
    
date_start = args.date_start
date_end = args.date_end

science_metadata = []
standard_metadata = []
request_id = []
obstype = []
# The "day" of the beginning of the observation
day_obs = []
# The exact time of the start of the observation
date_obs = []
instrume_list = []

# first north (en06) then south (en12)
# Changing this to either north or south, based on the user choice, for ease

if telescope == "south":
    instrume = "en12"
elif telescope == "north":
    instrume = "en06"
else :
    print("Choose a telescope! Quitting...")
    #return

#for instrume in ["en06", "en12"]:
# Get the metadata
_science_metadata = get_metadata(
        authtoken=authtoken,
        limit=1000,
        INSTRUME=instrume,
        start=date_start,
        end=date_end,
        RLEVEL=0,
        #covers="POINT({} {})".format(ra, dec),
        propID=propID,
        target_name="590",
    )
#print(_science_metadata[0]["request_id"])
# There is only one spectrum from a single group of request_id
# check that each SPECTRUM is associated with an ARC and a LAMPFLAT
# mismatched can happen if acquisition fail or being overrid by TOO
_request_id = []
for metadata in _science_metadata:
    _request_id.append(metadata["request_id"])
#print(_request_id)
c = Counter(_request_id)
mask = np.ones_like(_request_id, dtype=bool)
for i in set(_request_id):
        if c[i] < 3:
            mask[np.argwhere(np.array(_request_id) == i)] = False

_science_metadata = list(np.array(_science_metadata)[mask])
#print(_science_metadata)
for metadata in _science_metadata:
        request_id.append(metadata["request_id"])
        obstype.append(metadata["OBSTYPE"])
        day_obs.append(metadata["DAY_OBS"])
        date_obs.append(metadata["DATE_OBS"][:-1])
        instrume_list.append(instrume)

science_metadata += _science_metadata

request_id = np.array(request_id)
obstype = np.array(obstype)
day_obs = np.array(day_obs)
date_obs = np.array(date_obs)
instrume_list = np.array(instrume_list)

# Get the request_id that contains useful spectral data
request_id_science = [
        i for i, j in zip(request_id, obstype) if j == "SPECTRUM"
    ]

science_metadata = [
        v for v in science_metadata if v["request_id"] in request_id_science
    ]

# Look for the standard frames taken closest to the time of observation
for date, instrume, day in [
        (i, j, l)
        for i, j, k, l in zip(date_obs, instrume_list, day_obs, obstype)
        if l == "SPECTRUM"
    ]:
    _date = datetime.fromisoformat(date)
    day_range = 0.5
    if _date > datetime.fromisoformat("2021-09-01T00:00:00.000"):
        prop_id = "FLOYDS standards"
    else:
        if instrume == "en06":
            prop_id = "OGG_calib"
        else:
            prop_id = "COJ_calib"
    _standard_metadata = []
    while _standard_metadata == []:
        # Get the standard star
        _standard_metadata = get_metadata(
                authtoken=authtoken,
                limit=1000,
                INSTRUME=instrume,
                start=(_date - timedelta(days=day_range)).isoformat(),
                end=(_date + timedelta(days=day_range)).isoformat(),
                PROPID=prop_id,
                OBSTYPE="SPECTRUM",
                RLEVEL=0,
        )
        day_range += 1
        
    # If there are more than 1 frame returned, check for matching day_obs
    # first and then check for the one with the least time difference
    _standard_metadata_day_obs = []
    _time = []
    if len(_standard_metadata) > 1:
        # first check for the day_obs
        for meta in _standard_metadata:
            print(meta["DAY_OBS"])
            if np.in1d(meta["DAY_OBS"], day_obs):
                _standard_metadata_day_obs.append(meta)
            if len(_standard_metadata_day_obs) > 1:
                # Get all the timestamps
                for meta in _standard_metadata_day_obs:
                    _time.append(
                        datetime.fromisoformat(
                            meta["observation_date"][:-1]
                        )
                    )
                _time_diff = [i - _date for i in _time]
                min_time_idx = np.argmin(np.abs(np.array(_time_diff)))
            elif len(_standard_metadata_day_obs) == 1:
                min_time_idx = 0
                _time.append(
                    datetime.fromisoformat(meta["observation_date"][:-1])
                )
            else:
                # Get all the timestamps
                for meta in _standard_metadata:
                    print(meta["DAY_OBS"])
                    _time.append(
                        datetime.fromisoformat(
                            meta["observation_date"][:-1]
                        )
                    )
                    _time_diff = [i - _date for i in _time]
                    min_time_idx = np.argmin(np.abs(np.array(_time_diff)))
    else:
        min_time_idx = 0
        _time.append(
            datetime.fromisoformat(
                _standard_metadata[0]["observation_date"][:-1]
            )
        )
        
    # make sure there is a light frame...
    standard_spectrum_metadata = []
    while standard_spectrum_metadata == []:
        standard_spectrum_metadata = get_metadata(
            authtoken=authtoken,
            limit=1000,
            INSTRUME=instrume,
            start=(
                _time[min_time_idx] - timedelta(minutes=60)
            ).isoformat(),
            end=(_time[min_time_idx] + timedelta(minutes=60)).isoformat(),
            PROPID=prop_id,
            OBSTYPE="SPECTRUM",
            RLEVEL=0,
        )
        if len(standard_spectrum_metadata) > 1:
            standard_spectrum_metadata = [
                d for d in standard_spectrum_metadata
            ][0]
            
    # make sure there is an arc...
    standard_arc_metadata = []
    day_range = 0.0
    while standard_arc_metadata == []:
        standard_arc_metadata = get_metadata(
            authtoken=authtoken,
            limit=1000,
            INSTRUME=instrume,
            start=(
                _time[min_time_idx]
                - timedelta(minutes=60)
                - timedelta(days=day_range)
            ).isoformat(),
            end=(
                _time[min_time_idx]
                + timedelta(minutes=60)
                + timedelta(days=day_range)
            ).isoformat(),
            PROPID=prop_id,
            OBSTYPE="ARC",
            RLEVEL=0,
        )
        day_range += 1
        if len(standard_arc_metadata) > 1:
            standard_arc_metadata = [d for d in standard_arc_metadata][0]
    # make sure there is a flat...
    standard_flat_metadata = []
    day_range = 0.0
    while standard_flat_metadata == []:
        standard_flat_metadata = get_metadata(
            authtoken=authtoken,
            limit=1000,
            INSTRUME=instrume,
            start=(
                _time[min_time_idx]
                - timedelta(minutes=60)
                - timedelta(days=day_range)
            ).isoformat(),
            end=(
                _time[min_time_idx]
                + timedelta(minutes=60)
                + timedelta(days=day_range)
            ).isoformat(),
            PROPID=prop_id,
            OBSTYPE="LAMPFLAT",
            RLEVEL=0,
        )
        day_range += 1
        if len(standard_flat_metadata) > 1:
            standard_flat_metadata = [d for d in standard_flat_metadata][0]

    standard_metadata += standard_spectrum_metadata
    standard_metadata += standard_arc_metadata
    standard_metadata += standard_flat_metadata



# Pack the light, flat & arc fro science and standard as a dictionary item
target_list = {}

# Download the data and distinguish north-south for the rectification
science_filepath_list = []
science_filename_list = []
science_hemisphere_list = []
standard_filepath_list = []
standard_filename_list = []
standard_hemisphere_list = []

# Download the science and the best matched standard frames
for science_frame, standard_frame in zip(
    science_metadata, standard_metadata
):
    if (
        (".tar.gz" not in science_frame["filename"])
        or ("fits.fz" not in science_frame["filename"])
    ) & (science_frame["request_id"] in request_id_science):
        [
            x.append(y)
            for x, y in zip(
                [science_filepath_list, science_filename_list],
                download_frame(
                    frame=science_frame,
                    base_directory=input_folder,
                    no_date=False,
                ),
            )
        ]
        # Arrange into the target list
        if science_frame["request_id"] not in target_list.keys():
            target_list[science_frame["request_id"]] = {
                "science": {},
                "standard": {},
            }
        if science_frame["instrument_id"] == "en06":
            science_hemisphere_list.append("north")
            target_list[science_frame["request_id"]]["science"][
                "hemisphere"
            ] = "north"
        elif science_frame["instrument_id"] == "en12":
            science_hemisphere_list.append("south")
            target_list[science_frame["request_id"]]["science"][
                "hemisphere"
            ] = "south"
        else:
            print("This frame is not generated by FLOYDS.")
        target_list[science_frame["request_id"]]["science"][
            science_frame["OBSTYPE"]
        ] = science_frame["filename"]
        target_list[science_frame["request_id"]]["science"][
            "DAY_OBS"
        ] = science_frame["DAY_OBS"]
        target_list[science_frame["request_id"]]["science"][
            "OBJECT"
        ] = science_frame["OBJECT"]
        #print(target_list)
        # Download the standard frames to the SCIENCE FOLDER
        #print(standard_frame)
        try :
            if (".tar.gz" not in standard_frame["filename"]) or (
                ".fits.fz" not in standard_frame["filename"]
            ):
                [
                    x.append(y)
                    for x, y in zip(
                        [standard_filepath_list, standard_filename_list],
                        download_frame(
                            frame=standard_frame,
                            base_directory=os.path.join(
                                input_folder,
                                science_frame["DAY_OBS"].replace("-", ""),
                            ),
                            no_date=True,
                        ),
                    )
                ]
                if standard_frame["instrument_id"] == "en06":
                    standard_hemisphere_list.append("north")
                    target_list[science_frame["request_id"]]["standard"][
                        "hemisphere"
                    ] = "north"
                elif standard_frame["instrument_id"] == "en12":
                    standard_hemisphere_list.append("south")
                    target_list[science_frame["request_id"]]["standard"][
                        "hemisphere"
                    ] = "south"
                else:
                    print("This frame is not generated by FLOYDS.")
                target_list[science_frame["request_id"]]["standard"][
                    standard_frame["OBSTYPE"]
                ] = standard_frame["filename"]
                target_list[science_frame["request_id"]]["standard"][
                    "DAY_OBS"
                ] = standard_frame["DAY_OBS"]
                target_list[science_frame["request_id"]]["standard"][
                    "OBJECT"
                ] = standard_frame["OBJECT"]
        except TypeError:
            continue

for x, y in zip(science_filepath_list, science_filename_list):
    print("Science frame: {} is downloaded.".format(x + os.sep + y))

for x, y in zip(standard_filepath_list, standard_filename_list):
    print("Standard frame: {} is downloaded.".format(x + os.sep + y))

print(target_list)
