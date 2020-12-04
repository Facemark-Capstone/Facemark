#!/usr/bin/env python
# coding: utf-8

import dlib
import cv2
import numpy as np
from skimage import filters
from skimage.color import rgb2gray
from scipy import ndimage as nd
import math

detector = dlib.get_frontal_face_detector()
predictor = dlib.shape_predictor('shape_predictor_68_face_landmarks.dat')

def extractLandmark(image) :
    face = detector(image, 1)[0]
    landmarks = predictor(image, face)
    points = landmarks.parts()
    return points

def readImage(imageName) :
    #read image
    image = cv2.imread(imageName) 
    image = cv2.resize(image,(768,1024))
    image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    
    #extract landmark
    points = extractLandmark(image)
    return image, points


def calAngle(points) :
    midx = (points[39].x + points[42].x) / 2
    midy = (points[39].y + points[42].y) / 2
    midx = midx - points[33].x
    midy = midy - points[33].y
    return (90-math.acos(midx*1 / math.sqrt(midx*midx + midy*midy))*180/math.pi)

def rotate(origin, point, angle):
    angle = angle * -1
    ox, oy = origin
    px, py = point
    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy


def rotateImage(image, points, angle) :
    h,w = image.shape[:2]
    rotatedImage = np.copy(image)
    m1 = cv2.getRotationMatrix2D((w/2,h/2), angle,1)
    rotatedImage = cv2.warpAffine(image,m1,(w,h))

    for p in points :
        p.x, p.y = list(map(int, rotate((w/2,h/2),(p.x,p.y),math.radians(angle))))
    return rotatedImage, points


def cutImage(image, points):
    centerPoint = points[29]
    dx = math.ceil((points[14].x - points[27].x))
    dy = math.ceil((points[9].y - points[29].y))
    image = image[centerPoint.y - dy:centerPoint.y + dy + 1, centerPoint.x - dx:centerPoint.x + dx + 1]
    return image, dx


def sobleFilter(image) :
    image = filters.gaussian(image, sigma=1.5)
    gray_image = rgb2gray(image)
    grad_x = nd.convolve(gray_image, np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]])) 
    grad_y = nd.convolve(gray_image, np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]]))
    mag = np.sqrt(np.power(grad_x,2)+np.power(grad_y,2))
    theta = np.arctan2(np.power(grad_y,2), np.power(grad_x,2))
    return mag, theta


def to_string(landmark, DD, EOS) :
    data = ""
    for p in landmark :
        data = data + str(p.x) + " "
        data = data + str(p.y) + " "
    data = data + str(DD) + " "
    data = data + str(EOS) + " "
    return data

def cal_DD_EOS(image) : 
    mag, theta = sobleFilter(image)
    y_size = image.shape[0]
    x_size = image.shape[1]
    cnt = y_size * x_size / 2
    DD = 0
    for y in range(0,y_size) : 
        for x in range(0,dx) : 
            DD = DD + abs(mag[y][x] - mag[y][x_size-1-x])
    DD = round(DD / cnt,11)

    EOS = 0
    for y in range(0,y_size) : 
        for x in range(0,dx) :
            EOS = EOS + abs(theta[y][x] - theta[y][x_size-1-x])
    EOS = round(math.cos(EOS / cnt),11)
    return DD, EOS

image, points = readImage("0 (14).jpg") #image address
Landmark = points
angle = calAngle(points) 
image, points = rotateImage(image,points, angle) 
image, dx = cutImage(image, points) 
image = cv2.resize(image,(600,800)) 
DD, EOS =  cal_DD_EOS(image)


CPP_Input = to_string(Landmark, DD, EOS)


#output : CPP_Input(type: string)
print(CPP_Input)

