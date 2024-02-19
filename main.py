import image_reader
import copy
import numpy as np

class KDTree:
    def __init__(self, data: np.ndarray, overallShape, offset=(0,0),):
        self.depth = 0
        self.data = data
        self.offset = offset
        self.overallShape = overallShape
        self.left = None
        self.right = None
        self.splitAxis = None

    def split(self):
        if (self.depth == 0):
            a, aOffset, b, bOffset, self.splitAxis = medianCut(self.data, self.offset, self.overallShape)
            print(a.shape)
            print(b.shape)
            self.left = KDTree(a, self.overallShape, aOffset)
            self.right = KDTree(b, self.overallShape,bOffset)
        else: 
            self.left.split()
            self.right.split()
        self.depth += 1
        
    def size(self):
        return self.depth

    def getlightSources(self):
        if self.depth == 0:
            split_y, split_x = findSplitPoint(self.data, getWeightFunction(self.offset, self.overallShape))
            return [((self.offset[0] + split_y) ,(self.offset[1] + split_x))]
        else:
            return self.left.getlightSources() + self.right.getlightSources()
        
    def getLines(self):
        if self.depth == 0:
            return []
        else:
            if self.splitAxis == 0:
                length = self.data.shape[1] 
            if self.splitAxis == 1:
                length = self.data.shape[0]
            line = (self.right.offset,self.splitAxis,length)
            return  [line] + self.left.getLines() + self.right.getLines()



def getWeightFunction(offset, shape):
    weightFunc = lambda y, x, value: ((offset[0] + y)/shape[0])*np.pi * value
    return weightFunc
        
        

def findSplitPoint(data: np.ndarray, weightFunc):
    # Calculate weight function values for each point in the data
    weight_values = np.zeros(data.shape[:2])
    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            weight_values[y, x] = weightFunc(y, x, np.sum(data[y, x]))

    # Calculate cumulative sums along x and y axes
    cum_sum_x = np.cumsum(weight_values.sum(axis=0))
    cum_sum_y = np.cumsum(weight_values.sum(axis=1))

    # Find indices where cumulative sums are closest to half of total sum
    total_sum_x = cum_sum_x[-1]
    total_sum_y = cum_sum_y[-1]
    half_sum_x = total_sum_x / 2.0
    half_sum_y = total_sum_y / 2.0
    split_x = np.argmin(np.abs(cum_sum_x - half_sum_x))
    split_y = np.argmin(np.abs(cum_sum_y - half_sum_y))
    print(split_y, split_x)

    return split_y, split_x

def medianCut(data: np.ndarray, offset, overallShape):

    split_y, split_x, = findSplitPoint(data, getWeightFunction(offset, overallShape))

    splitAxis = None

    # Check if height > width - Split on y axis
    if data.shape[0] > data.shape[1]:
        splitIndex = split_y
        a = data[:splitIndex,:,:]
        b = data[splitIndex:,:,:]

        aOffset = (offset[0], offset[1])
        bOffset = (offset[0] + splitIndex, offset[1])
        splitAxis = 0
    # Otherwise, width > height - Split on x
    else:
        splitIndex = split_x
        
        a = data[:,:splitIndex,:]
        b = data[:,splitIndex:,:]
        aOffset = (offset[0], offset[1])
        bOffset = (offset[0] ,offset[1] + splitIndex)
        splitAxis = 1

    return a, aOffset, b, bOffset, splitAxis

def draw(image, lines, points,thickness=1):
    # Draw lines
    for line in lines:
        start, axis, length = line
        if axis == 1:  # Horizontal line
            end = (start[0] + length, start[1])
            image[int(start[0]):int(end[0]),int(start[1]-thickness):int(start[1]+thickness)] = [255, 255, 255]  # White
        else:  # Assuming '0' represents vertical
            end = (start[0], start[1] + length)
            image[int(start[0]-thickness):int(start[0]+thickness),int(start[1]):int(end[1])] = [255, 255, 255]  # White

    # Draw points
    for point in points:
        # Ensure the slice indices are integers
        x, y = int(point[1]), int(point[0])
        image[y-2:y+3, x-2:x+3] = [0, 0, 255]  # Blue

    return image


def apply_gamma_correction(image: np.ndarray, gamma_correction: float) -> np.ndarray:
    """
    Applies gamma correction to the image.
    """
    # Create a copy of the image to avoid modifying the original data
    normalized_image = np.copy(image)
        
    # Apply gamma correction
    corrected_image = np.power(normalized_image, gamma_correction)
    return corrected_image

def convert_pfm_to_ppm(pfm_image, gamma_correction: float):
    """
    Reads a PFM image, applies gamma correction, and prepares it for PPM format.
    """
    # Apply gamma correction
    gamma_corrected_image = apply_gamma_correction(pfm_image, 1.0/gamma_correction)
    
    # Convert to the appropriate range for PPM (0-255)
    ppm_ready_image = np.clip(gamma_corrected_image * 255, 0, 255).astype(np.uint8)
    
    return ppm_ready_image

def main():
    # Assuming image_reader.read_pfm and image_reader.write_ppm are correctly defined elsewhere
    image = image_reader.read_pfm("GraceCathedral/grace_latlong.pfm")
    
    cutTree = KDTree(image,image.shape)
    print(image.shape)
    loops = 2
    for i in range(loops):
        cutTree.split()
        splits = pow(2,i+1)

        # points = cutTree.getlightSources()
        lines = cutTree.getLines()
        points = []
        
        out_image = copy.deepcopy(image)
        drawn_image = draw(out_image, lines, points)

        gamma_corrected_image = convert_pfm_to_ppm(drawn_image, 2.2)
        name = "Out/test_copy_" + str(splits) +  ".ppm"
        image_reader.write_ppm(gamma_corrected_image, name)


main()


def viewPFM():
    image = image_reader.read_pfm("pbrt/simple_sphere.pfm")
    gamma_corrected_image = convert_pfm_to_ppm(image, 1.0)
    image_reader.write_ppm(gamma_corrected_image, "Out/simple_sphere.ppm")
    
# viewPFM()