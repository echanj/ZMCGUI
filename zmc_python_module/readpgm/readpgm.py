import re
import numpy

# 
# got this form stack echange but there was a bug
# how to get 16 bit to work ?? use dt = numpy.uint16
#                                  dt = dt.newbyteorder('>')
#
# I have edited this and seems to work for the .pgm we use for diffuse work  
#
# note values go from black to white but in YELL the .h5 format may differ 
#

def read_pgm(filename, byteorder='>'):
    """Return image data from a raw PGM file as numpy array.

    Format specification: http://netpbm.sourceforge.net/doc/pgm.html

    """
    with open(filename, 'rb') as f:
        buffer = f.read()
    try:
        header, width, height, maxval = re.search(
            b"(^P5\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n])*"
            b"(\d+)\s(?:\s*#.*[\r\n]\s)*)", buffer).groups()
    except AttributeError:
        raise ValueError("Not a raw PGM file: '%s'" % filename)
    # dt = numpy.dtype(numpy.uint16)
    dt = numpy.dtype(numpy.uint16)
    dt = dt.newbyteorder('>')
    return numpy.frombuffer(buffer,
                            # dtype='u1' if int(maxval) < 256 else byteorder+'u2',
                            dtype='u1' if int(maxval) < 256 else dt,
                            # dtype=dt,
                            count=int(width)*int(height),
                            offset=len(header)
                            ).reshape((int(height), int(width)))

if __name__ == "__main__":
    from matplotlib import pyplot
    image = read_pgm("h_0.66_l_16bit_box1.pgm", byteorder='>')
    print image[1,:]
    pyplot.imshow(image, pyplot.cm.gray)
    pyplot.show()




