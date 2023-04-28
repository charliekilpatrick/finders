# finders

Finder-generation code, originally developed by Kaew Tinyanont and extended by Charlie Kilpatrick

Usage: python prepare_obs_run.py filename telescope

Filename is your target list and must be formatted with four columns with name, ra, dec, magnitude.  See target_list.txt for an example input file formatted so the code can parse.

Telescope can be Keck, Lick, or SOAR.

Outputs will be a reformatted target list (in Keck format) that includes the target and offset stars as well as finder charts in a "finders/" directory.
