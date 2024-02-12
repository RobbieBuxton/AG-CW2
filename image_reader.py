import numpy as np
import os
import sys

def read_pfm(filepath: str)->np.ndarray:
    with open(filepath, 'rb') as file:
        file_contents = file.read()
        header, other_contents = file_contents.split(chr(10).encode("ascii"), 1)
        header = header.decode("ascii")
        if header[0] != "P" or header[1].lower() != "f":
            raise TypeError(f"File header is malformed, expected one of PF or Pf, got {header}")
        image_dims, scale, _ = other_contents.split(chr(10).encode("ascii"), 2)
        file.seek(len(header) + len(image_dims) + len(scale) + 3)
        image_dims = tuple(map(int, image_dims.decode("ascii").split(" ")))[::-1]
        scale = float(scale)
        little_endian = scale < 0
        color = (header[1] == "F")
        np_dims = (*image_dims, 3 if color else 1)
        np_dtype = "<f4" if little_endian else ">f4"
        np_buffer = bytearray(file.read(np_dims[0] * np_dims[1] * np_dims[2] * 4))
        result = np.ndarray(shape=np_dims, dtype=np_dtype, buffer=np_buffer)
        result = result[::-1, :, :]
        return result

def write_pfm(image: np.ndarray, filepath: str):
    if image.dtype not in [np.float32]:
        raise TypeError("Expected type of image to be one of np.float32, got {}".format(image.dtype))
    with open(filepath, 'wb') as file:
        image_float32 = image.astype(np.float32)
        scale = image_float32.max()
        if scale == 0:
            scale = 1.0
        header = ("Pf" if image_float32.shape[2] == 1 else "PF") + chr(10)
        file.write(header.encode('ascii'))
        file.write((" ".join(map(str, image_float32.shape[:2][::-1]))+chr(10)).encode("ascii"))
        if image_float32.dtype.byteorder == "<":
            file_byteorder = "-1.0"
        elif image_float32.dtype.byteorder == ">":
            file_byteorder = "1.0"
        else:
            if sys.byteorder.lower() == "little":
                file_byteorder = "-1.0"
            else:
                file_byteorder = "1.0"
        file_byteorder = file_byteorder + chr(10)
        file.write(file_byteorder.encode("ascii"))
        file.write(image_float32[::-1, :, :].tobytes())
		
def __skip_whitespace_and_comments(file):
    c = file.read(1).decode("ascii")
    currently_a_comment = False
    while c.isspace() or currently_a_comment:
        c = file.read(1).decode("ascii")
        if not currently_a_comment and c == "#":
            currently_a_comment = True
        elif c == "\n":
            currently_a_comment = False
    return c
        
def read_ppm(filepath: str)->np.ndarray:
    with open(filepath, 'rb') as file:
        header = file.read(2)
        header = header.decode("ascii")
        data_format = int(header[1])
        c = __skip_whitespace_and_comments(file)
        width = c
        while not c.isspace():
            c = file.read(1).decode("ascii")
            width = width + c
        c = __skip_whitespace_and_comments(file)
        height = c
        while not c.isspace():
            c = file.read(1).decode("ascii")
            height = height + c
        c = __skip_whitespace_and_comments(file)
        width = int(width)
        height = int(height)
        maxval = c
        while not c.isspace():
            c = file.read(1).decode("ascii")
            maxval = maxval + c
        maxval = int(maxval)
        c = __skip_whitespace_and_comments(file)
        file.seek(-1, os.SEEK_CUR)
        np_dims = (height, width)
        if maxval > 255:
            np_dtype = np.uint16
        else:
            np_dtype = np.uint8
        if data_format == 6:
            # This is a ppm formatted file
            np_dims = (*np_dims, 3)
            np_buffer = bytearray(file.read(width * height * 3))
            result = np.ndarray(shape=np_dims, dtype=np_dtype, buffer=np_buffer)
        elif data_format == 5:
            # This is a pgm formatted file
            np_dims = (*np_dims, 1)
            np_buffer = bytearray(file.read(width * height * 1))
            result = np.ndarray(shape=np_dims, dtype=np_dtype, buffer=np_buffer)
        elif data_format == 2:
            # This is a pgm formatted file
            np_dims = (*np_dims, 1)
            values = file.read()
            values = values.decode("ascii").split()
            values = list(map(float, values))
            result = np.array(values).reshape(*np_dims)
        elif data_format == 3:
            # This is a ppm formatted file
            np_dims = (*np_dims, 3)
            values = file.read()
            values = values.decode("ascii").split()
            values = list(map(float, values))
            result = np.array(values).reshape(*np_dims)
        return result
        
def write_ppm(image: np.ndarray, filepath: str):
    if image.dtype not in [np.uint8, np.uint16]:
        raise TypeError("Expected type of image to be one of np.uint8 or np.uint16, got {}".format(image.dtype))
    with open(filepath, 'wb') as file:
        data_format = ""
        if len(image.shape) == 2:
            image = image[:,:,np.newaxis]
        if image.shape[-1] == 1:
            data_format = "5"
        else:
            data_format = "6"
        file.write(("P" + data_format + chr(10)).encode("ascii"))
        file.write((str(image.shape[1]) + " " + str(image.shape[0]) + chr(10)).encode("ascii"))
        file.write((str(image.max())+chr(10)).encode("ascii"))
        if image.max() < 256:
            image_uint8 = image.astype(np.uint8)
            file.write(image_uint8.tobytes())
        else:
            image_uint16 = np.clip(image, 0, 65535).astype(np.uint16)
            file.write(image_uint16.tobytes())
        