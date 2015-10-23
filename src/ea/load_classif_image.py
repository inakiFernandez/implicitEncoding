# -*- coding: utf-8 -*-
"""
Read and load a black and white image into a table associating
pairs of x and y coordinates (starting from top-leftmost pixel) to the
corresponding color/label (white=0, black=1)
"""
import png
import sys


def read_image(name, folder="../../../image_db"):
    r = png.Reader(folder + '/' + name)
    f = r.asDirect()

    width = f[0]
    height = f[1]

    pixels = f[2]
    img_prop = f[3]

    pixel_byte_width = 3

    grey_pixels = [-1] * (width * height)
    aux_pixels = list(pixels)

    i = 0
    for y in xrange(height):
        for x in xrange(width):
            grey_pixels[i] = aux_pixels[y][x * pixel_byte_width]
            i = i + 1

    return width, height, grey_pixels, img_prop


if __name__ == "__main__":

    filename = ""
    if len(sys.argv) > 1:
        filename = sys.argv[1]

    w, h, pix, metadat = read_image(filename)

    print w
    print h
    print pix
    print metadat
