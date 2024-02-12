import image_reader

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
            return [((self.offset[0] + self.data.shape[0]/2) ,(self.offset[1] + self.data.shape[1]/2))]
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

def medianCut(data: np.ndarray, offset, shape):
    
    splitAxis = None

    # Check if width > height
    if data.shape[1] > data.shape[0]:
        # Split on width
        sumArray = np.array(np.zeros(data.shape[1])) 
        for a in range(data.shape[1]):
            for b in range(data.shape[0]):
                scaleFactor = ((offset[1] + b)/shape[0])*np.pi
                sumArray[a] += (scaleFactor * sum(data[b,a])/3)

        splitIndex = find_split_index(sumArray)
        
        a = data[:,:splitIndex,:]
        b = data[:,splitIndex:,:]
        aOffset = (offset[0], offset[1])
        bOffset = (offset[0], offset[1] + splitIndex)
        splitAxis = 1

    else:
        # Split on height
        sumArray = np.array(np.zeros(data.shape[0]))
        for b in range(data.shape[0]):
            for a in range(data.shape[1]):
                scaleFactor = ((offset[0] + a)/shape[1])*np.pi
                sumArray[b] += (scaleFactor * sum(data[b,a])/3)

        splitIndex = find_split_index(sumArray)
        
        a = data[:splitIndex,:,:]
        b = data[splitIndex:,:,:]

        aOffset = (offset[0], offset[1])
        bOffset = (offset[0] + splitIndex, offset[1])
        splitAxis = 0

    # print(data.shape)  
    # print(a.shape)
    # print(b.shape)
    return a, aOffset, b, bOffset, splitAxis


# def draw(image, lines, points):
#     #Image is a numpy array of shape (height, width, 3) where 3 is the RGB channels
#     #Lines is a list of tuples (start, axis, length)
#     #Points is a list of tuples (x,y)

#     #Write me a function that draws the lines and points on the image of the numpy array
#     #The lines should be drawn in white (1,1,1) and be 5 pixels wide
#     #The points should be in blue (0,0,1) and be 5*5 pixels 
#     pass


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




def main():
    # Assuming image_reader.read_pfm and image_reader.write_ppm are correctly defined elsewhere
    image = image_reader.read_pfm("GraceCathedral/grace_latlong.pfm")

    cutTree = KDTree(image,image.shape)

    loops = 4
    for _ in range(loops):
        cutTree.split()

    print(pow(2,loops))

    points = cutTree.getlightSources()
    lines = cutTree.getLines()
    print(lines)
    print(points)
    drawn_image = draw(image, lines, points)

    gamma_corrected_image = convert_pfm_to_ppm(drawn_image, 0.2)
    image_reader.write_ppm(gamma_corrected_image, "Out/test_copy.ppm")



def find_split_index(arr):
    """
    Finds the index to split the array into two halves such that the sum of the halves is as equal as possible.
    
    Parameters:
    arr (np.array): The input numpy array to split.
    
    Returns:
    int: The index after which the array should be split.
    """
    # Calculate the cumulative sum of the array
    cumulative_sum = np.cumsum(arr)

    # Find the total sum and calculate half of it
    total_sum = cumulative_sum[-1]
    half_sum = total_sum / 2

    # Find the index where cumulative sum is closest to half of the total sum
    index = np.abs(cumulative_sum - half_sum).argmin()

    # Adjust index if necessary to ensure the split is as equal as possible
    if cumulative_sum[index] < half_sum:
        index += 1

    return index

def apply_gamma_correction(image: np.ndarray, gamma_correction: float) -> np.ndarray:
    """
    Applies gamma correction to the image.
    """
    # Normalize the image to the range 0-1
    normalized_image = image - image.min()
    normalized_image /= normalized_image.max()
    # Apply gamma correction directly as a float value, not as a callable
    corrected_image = np.power(normalized_image, gamma_correction)
    return corrected_image

def convert_pfm_to_ppm(pfm_image, gamma_correction: float):
    """
    Reads a PFM image, applies gamma correction, and prepares it for PPM format.
    """
    # Apply gamma correction
    gamma_corrected_image = apply_gamma_correction(pfm_image, gamma_correction)
    
    # Convert to the appropriate range for PPM (0-255)
    ppm_ready_image = np.clip(gamma_corrected_image * 255, 0, 255).astype(np.uint8)
    
    return ppm_ready_image


def viewPFM():
    image = image_reader.read_pfm("pbrt/simple_sphere.pfm")
    gamma_corrected_image = convert_pfm_to_ppm(image, 1.0)
    image_reader.write_ppm(gamma_corrected_image, "Out/simple_sphere.ppm")

viewPFM()
# main()