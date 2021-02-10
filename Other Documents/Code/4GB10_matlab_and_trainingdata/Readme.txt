Data folder contains example measurement

matlab folder contains three files
	ConvertFiles.m
	NasaUseExample.m
	ReadCaseDataExample.m
and subfolder General and General/Nasa

ConvertFiles.m 			is needed to convert txt files where a ',' is used as decimal separator to ones where '.' is used. 
					Use it with caution since it will overwrite the original files. Make backups first.

NasaUseExample. m 		is an example how to use the Nasa thermal database.

ReadCaseDataExample.m 	processes a set of measurement txt files in a given directory and plots the encoder signal 
					and the pressure signal as function of time in a single figure (using subplot commmand).
