cd ../../data; 
tar cvf all_configs.tar * --exclude *.mp4 --exclude *.mov --exclude *.tiff --exclude *.TIF --exclude *.tif --exclude *.png --exclude config;
scp all_configs.tar $1:
