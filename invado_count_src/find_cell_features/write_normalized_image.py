#!/usr/bin/python

import Image
import sys

#open the image name provided in command line argument one
im = Image.open(sys.argv[1])

#we will cut out the bottom fifty pixels as they are often contain noise
box = (0,0,im.size[0], im.size[1] - 50)
region = im.crop(box)

#now we convert the file to a format that PNG is happy with
out = region.convert('I')

#and save the final image into command line argument two
out.save(sys.argv[2])
