import os
import json
import glob
from pathlib import Path
from appdirs import user_data_dir

from ccinput.__init__ import __version__
from ccinput.calculation import Calculation, Parameters
from ccinput.exceptions import InvalidParameter
from ccinput.utilities import warn

data_dir = user_data_dir("ccinput", "CYLlab")


def warn_invalid_json(preset_name):
    warn(f"Invalid JSON found for existing preset {preset_name} (in {data_dir})")


def save_preset(args, default_args):
    preset_name = args.save

    if len(preset_name.strip()) == 0:
        raise InvalidParameter(
            "The preset name must contain at least one valid character"
        )

    params = vars(args)
    default_params = vars(default_args)

    del params["save"]
    del default_params["save"]

    for k, v in default_params.items():
        if v == params[k] or params[k] is None:
            del params[k]

    if "xyz" in params:
        warn("Ignoring the xyz structure")
        del params["xyz"]

    if "file" in params:
        warn("Ignoring the input structure")
        del params["file"]

    if "output" in params:
        warn("Ignoring the output name")
        del params["output"]

    if "name" in params:
        warn("Ignoring the name")
        del params["name"]

    ## Test the validity of the parameters?

    if not os.path.isdir(data_dir):
        Path(data_dir).mkdir(parents=True, exist_ok=True)

    preset_json = {}

    path = os.path.join(data_dir, f"{preset_name}.preset")
    if os.path.isfile(path):
        try:
            with open(path) as f:
                preset_json = json.load(f)
        except json.decoder.JSONDecodeError:
            warn_invalid_json()

    for k, v in params.items():
        preset_json[k] = v

    preset_json["version"] = __version__

    with open(path, "w") as out:
        json.dump(preset_json, out, indent=4)

    return preset_name


def load_preset(preset_name):
    if len(preset_name.strip()) == 0:
        raise InvalidParameter(
            "The preset name must contain at least one valid character"
        )

    path = os.path.join(data_dir, f"{preset_name}.preset")

    if os.path.isfile(path):
        try:
            with open(path) as f:
                preset_json = json.load(f)
        except json.decoder.JSONDecodeError:
            warn_invalid_json()
    else:
        raise InvalidParameter(
            "No preset found with the name {preset_name} (in {data_dir})"
        )

    return preset_json


def get_preset_names():
    _presets = list(glob.glob(os.path.join(data_dir, "*.preset")))
    return [os.path.splitext(os.path.basename(p))[0] for p in _presets]


def is_preset(preset_name):
    presets = get_preset_names()
    if f"{preset_name}" in presets:
        return True

    print(f"Unknown preset: '{preset_name}'")
    print("")
    return False


def list_presets():
    presets = get_preset_names()
    print("--- Available presets:")
    for p in presets:
        print(p)


def print_preset(preset_name):
    path = os.path.join(data_dir, f"{preset_name}.preset")

    if not os.path.isfile(path):
        raise InvalidParameter(
            f"No preset called {preset_name} found in your user directory ({data_dir})"
        )

    try:
        with open(path) as f:
            preset_json = json.load(f)
    except json.decoder.JSONDecodeError:
        warn_invalid_json(preset_name)
        return

    print(f"--- Saved preset '{preset_name}'")

    for k, v in preset_json.items():
        print(f"{k:<30}{v:<30}")
