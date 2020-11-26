#!/usr/bin/env python
# coding: utf-8

import tensorflow as tf
from tensorflow import keras
from keras.models import load_model
import numpy as np
import pandas as pd


#input : FromCPP (type : string)
floatStr = list(map(float, FromCPP))

#input string From CPP

df = pd.DataFrame(floatStr).T
X = df.drop([29], axis=1)
Y = df[29]
keras.backend.set_floatx('float64')
LoadModel = load_model('project/model.h5') #model dic
prediction = LoadModel.predict(X)
result = int(np.argmax(prediction))

rate = [df[24],df[25],df[26],df[27],df[28]]
sum = 0
for i in range(0,5) :
    sum = sum + floatStr[i+24]
score = (sum / 5) * 100

a = str(result) + ' '
a = a + str(score) + ' '
for i in range(0,5) :
    a = a + str(floatStr[i+24]) + ' '
print(a)
# output a (asy(1/2), score, eye, nose, mouth, jaw, face)
# if asy == 1 ---> Asymmetry 
#else asy == 2 ---> symmetry

