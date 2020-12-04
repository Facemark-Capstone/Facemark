#!/usr/bin/env python
# coding: utf-8

import tensorflow as tf
from tensorflow import keras
from keras.models import load_model
import numpy as np
import pandas as pd


#input : FromCPP (type : string)

FromCpp = FromCpp.split() #for split string to list
floatStr = list(map(float, FromCPP))


df = pd.DataFrame(floatStr[0:32]).T
X = df
keras.backend.set_floatx('float64')
LoadModel = load_model('model.h5') #model dic
prediction = LoadModel.predict(X)
result = int(np.argmax(prediction))

a = ''
a = str(result) + ' '
a = a + str(floatStr[37]) + ' '
for i in range(0,5) :
    a = a + str(floatStr[i+32]) + ' '
print(a)

# output a (asy(1/2), score, eye, nose, mouth, jaw, face)
# if asy == 1 ---> Asymmetry 
#else asy == 2 ---> symmetry

