#!/usr/bin/python

import Image
import sys

#open the image name provided in command line argument one
im = Image.open(sys.argv[1])

#we will cut out the pixels at the very edge as they often noisy
box = (1,1,im.size[0] - 1, im.size[1] - 1)
region = im.crop(box)

#now we convert the file to a format that PNG is happy with
out = region.convert('I')

#and save the final image into command line argument two
out.save(sys.argv[2])
