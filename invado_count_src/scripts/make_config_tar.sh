cd ../../data; 
tar cvf all_configs.tar * --exclude *.avi* --exclude *.tiff --exclude *.TIF --exclude *.tif --exclude *.png --exclude config;
if [ -n "$1" ]; then 
	scp all_configs.tar $1:
fi
