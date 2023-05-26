import sys
if len(sys.argv) < 3:
	print('Usage: python LRIS_exp_time.py individual_blue_exp_time #red_desired #blue_desired')
	print('If #blue_desired not provided, defaults to 1.\n\n')
	print('Example: python LRIS_exp_time.py 1200 2 \n--> will give you the individual exposure time for the red side if the blue side is 1x1200 sec, \nand you want to take 2 red exposures.')
else:
	blue_exp = int(sys.argv[1])
	num_red = int(sys.argv[2])
	if len(sys.argv) == 4:
		num_blue = int(sys.argv[3])
	else:
		num_blue = 1

	#Readout times from https://www2.keck.hawaii.edu/inst/lris/win_bin.html
	blue_readout = 54 #s, 1x1 binning, full frame 4 amps
	red_readout = 65 #x, Spec 2x1 2 amp

	### We want to make sure that the last frame of blue and red finish exposing at the same time so that we can slew. 
	###  (blue_exp + blue_readout) * num_blue - blue_readout = (red_exp + red_readout) * num_red - red_readout
	red_exp = (((blue_exp + blue_readout)*num_blue - blue_readout) + red_readout)/num_red -red_readout

	print("Appropriate exposure time for each red frame is:\n%d"%red_exp)