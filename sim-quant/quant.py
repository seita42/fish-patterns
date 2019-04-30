import sys

import numpy as np
import cv2
from PIL import Image


def quant(img_mono):
    """Quantify the lightness and complexity of an image

    Quantify the lightness and complexity of an image

    Args:
        img_mono: Numpy array of an image (8bit, binarized)

    Returns:
        lightness: The lightness of the image (white pixel/total pixel)
        pcs: The pattern complexity score (1 - area weighted mean circularity)
        ele_num: The number of pattern elements in the image
        ele_area_mean: The mean area of pattern elements
        ele_area_std: The standard deviation of area of pattern elements

    """

    # Lightness = <light part>/<total area> (0 to 255 -> 0 to 1)
    lightness = np.sum(img_mono)/(img_mono.shape[0]*img_mono.shape[1])/255.0

    if (lightness == 1.0 or lightness == 0.0):
        pcs = 0.0
        ele_num = 0
        ele_area_mean = 0.0
        ele_area_std = 0.0
    else:
        vals0 = _circula(img_mono)
        vals1 = _circula(cv2.bitwise_not(img_mono))
        if vals0[0] >= vals1[0]:
            w_circ, ele_num, ele_area_mean, ele_area_std = vals0
        else:
            w_circ, ele_num, ele_area_mean, ele_area_std = vals1

        # PatternComplexityScore = 1-(area weighted mean circularity) (0 to 1)
        pcs = 1-w_circ

    return lightness, pcs, ele_num, ele_area_mean, ele_area_std


def _circula(img_mono):

    # findContours
    _, contours, hierarchy = cv2.findContours(img_mono.copy(),
                                              cv2.RETR_CCOMP,
                                              cv2.CHAIN_APPROX_SIMPLE)
    # perimeter
    perimeters = np.array([cv2.arcLength(c, True) for c in contours])

    # contourArea
    areas = np.array([cv2.contourArea(c) for c in contours])

    for i, h in enumerate(hierarchy[0]):
        if h[3] < 0:        # [i] has no parent
            pass
        if h[3] >= 0:       # [i] has parent (parent = h[3])
            perimeters[h[3]] += perimeters[i]     # add child's perimeter
            areas[h[3]] -= areas[i]               # subtract child's area

    # parents only
    elements_perimeter = perimeters[hierarchy[0][:, 3] < 0]
    elements_area = areas[hierarchy[0][:, 3] < 0]

    # discard if area <= 0
    elements_perimeter = elements_perimeter[elements_area > 0]
    elements_area = elements_area[elements_area > 0]

    elements_circularity = 4.0*np.pi*elements_area/(elements_perimeter**2)

    w_circ = np.sum(elements_circularity*elements_area)/np.sum(elements_area)
    ele_num = len(elements_area)
    ele_area_mean = np.mean(elements_area)
    ele_area_std = np.std(elements_area)

    return w_circ, ele_num, ele_area_mean, ele_area_std


def binarize(img_gray,
             diameter=10,
             sigmaColor=150,
             sigmaSpace=150,
             a_thresh_type=2,
             a_thresh_blocksize=20,
             a_thresh_offset=0,
             a_thresh_offset2=0):
    """Binarize an image

    Binarize an image using OpenCV library

    Args:
        img_gray: Numpy array of an image (8bit, grayscale)
        diameter: param for bilateralFilter
        sigmaColor: param for bilateralFilter
        sigmaSpace: param for bilateralFilter
        a_thresh_type: 0, 1: Adaptive thresholding, 2: Otsu's thresholding
        a_thresh_blocksize: used for blocksize of adaptiveThreshold
        a_thresh_offset: constant for subtraction
        a_thresh_offset2: constant for negative subtraction (addition)

    Returns:
        img_mono: Numpy array of binarized image

    """

    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))

    img_blur = cv2.bilateralFilter(img_gray,
                                   diameter,
                                   sigmaColor,
                                   sigmaSpace)

    ret, img_mono0 = cv2.threshold(img_blur,
                                   127,
                                   255,
                                   cv2.THRESH_BINARY)

    if a_thresh_type == 0:
        img_mono = cv2.adaptiveThreshold(img_blur,
                                         255,
                                         cv2.ADAPTIVE_THRESH_MEAN_C,
                                         cv2.THRESH_BINARY,
                                         (a_thresh_blocksize+1)*2+1,
                                         a_thresh_offset - a_thresh_offset2)
    elif a_thresh_type == 1:
        img_mono = cv2.adaptiveThreshold(img_blur,
                                         255,
                                         cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                                         cv2.THRESH_BINARY,
                                         (a_thresh_blocksize+1)*2+1,
                                         a_thresh_offset - a_thresh_offset2)

    elif a_thresh_type == 2:
        ret, img_mono = cv2.threshold(img_blur,
                                      0,
                                      255,
                                      cv2.THRESH_BINARY+cv2.THRESH_OTSU)

    if np.sum(img_mono0) == 0:
        img_mono = img_mono0
    elif np.sum(img_mono0) == (255*img_mono0.shape[0]*img_mono0.shape[1]):
        img_mono = img_mono0

    img_closing = cv2.morphologyEx(img_mono, cv2.MORPH_CLOSE, kernel)
    img_mono = img_closing

    img_opening = cv2.morphologyEx(img_mono, cv2.MORPH_OPEN, kernel)
    img_mono = img_opening

    return img_mono


def crop_center(pil_img, crop_width, crop_height):

    img_width, img_height = pil_img.size
    
    left = (img_width - crop_width) // 2
    top = (img_height - crop_height) // 2
    right = (img_width + crop_width) // 2
    bottom = (img_height + crop_width) // 2
    
    return pil_img.crop((left, top, right, bottom))


def pil_img(img):

    return Image.fromarray(img)


if __name__ == '__main__':
    argv = sys.argv
    in_file = argv[1]
    img_mono = cv2.imread(in_file, cv2.IMREAD_GRAYSCALE)
    print(in_file, quant(img_mono))
