import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Input, Layer, Dense, Flatten, LayerNormalization, MultiHeadAttention, LSTM, GRU, Embedding, Add, Dropout, Multiply, Concatenate
import os
import time
from keras import backend
from IPython import display
import math


def mlp(x, hidden_units, dropout_rate=0.25, activation=tf.nn.relu):
    for units in hidden_units:
        x = Dense(units, activation=activation)(x)
        x = Dropout(dropout_rate)(x)
    return x
    

def MuRNN(lstm=True):
        
    input_layers = Input(shape=(12, 12))
    x = input_layers
    
    if lstm:
        x = LSTM(units=500, return_sequences=False, dropout=0.0, recurrent_dropout=0.0)(x)
    else:
        x = GRU(units=500, return_sequences=False, dropout=0.0, recurrent_dropout=0.0)(x)
            
    x = mlp(x, [250, 150, 12])

    return tf.keras.models.Model(input_layers, x)


if __name__ == "__main__":
    MuLSTM = MuRNN(lstm=True)
    MuLSTM.summary()

    MuGRU = MuRNN(lstm=False)
    MuGRU.summary()
