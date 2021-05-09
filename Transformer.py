import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Input, Layer, Dense, Conv1D, Flatten, TimeDistributed, LayerNormalization, MultiHeadAttention, LSTM, GRU, Embedding, Add, Dropout, Multiply, Concatenate
import os
import time
import keras
from keras import backend as K
from IPython import display
import math

def mlp(x, hidden_units, dropout_rate=0.25, activation=tf.nn.relu):
    for units in hidden_units:
        x = Dense(units, activation=activation)(x)
        x = Dropout(dropout_rate)(x)
    return x

class Time2Vec(Layer):
    def __init__(self, kernel_size=1):
        super(Time2Vec, self).__init__(trainable=True, name='Time2VecLayer')
        self.k = kernel_size
    
    def build(self, input_shape):
        # trend
        self.wb = self.add_weight(name='wb',shape=(input_shape[1],),initializer='uniform',trainable=True)
        self.bb = self.add_weight(name='bb',shape=(input_shape[1],),initializer='uniform',trainable=True)
        # periodic
        self.wa = self.add_weight(name='wa',shape=(1, input_shape[1], self.k),initializer='uniform',trainable=True)
        self.ba = self.add_weight(name='ba',shape=(1, input_shape[1], self.k),initializer='uniform',trainable=True)
        super(Time2Vec, self).build(input_shape)
    
    def call(self, inputs, **kwargs):
        bias = self.wb * inputs + self.bb
        dp = K.dot(inputs, self.wa) + self.ba
        wgts = K.sin(dp) # or K.cos(.)
        ret = K.concatenate([K.expand_dims(bias, -1), wgts], -1)
        ret = K.reshape(ret, (-1, inputs.shape[1]*(self.k+1)))
        return ret
    
    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[1]*(self.k + 1))
    
class AttentionBlock(keras.Model):
    def __init__(self, name='AttentionBlock', num_heads=2, head_size=128, ff_dim=None, dropout=0, **kwargs):
        super().__init__(name=name, **kwargs)

        if ff_dim is None:
            ff_dim = head_size

        self.attention = MultiHeadAttention(num_heads=num_heads, key_dim=head_size, dropout=dropout)
        self.attention_dropout = Dropout(dropout)
        self.attention_norm = LayerNormalization(epsilon=1e-6)

        self.ff_conv1 = Conv1D(filters=ff_dim, kernel_size=1, activation='relu')
        # self.ff_conv2 at build()
        self.ff_dropout = Dropout(dropout)
        self.ff_norm = LayerNormalization(epsilon=1e-6)

    def build(self, input_shape):
        self.ff_conv2 = Conv1D(filters=input_shape[-1], kernel_size=1) 

    def call(self, inputs):
        x = self.attention(inputs, inputs)
        x = self.attention_dropout(x)
        x = self.attention_norm(inputs + x)

        x = self.ff_conv1(x)
        x = self.ff_conv2(x)
        x = self.ff_dropout(x)

        x = self.ff_norm(inputs + x)
        return x

def MuTran(num_layers=4):
        
    input_layers = Input(shape=(12, 12))
    
    time_embedding = TimeDistributed(Time2Vec(kernel_size=1))(input_layers)
    output = Concatenate(axis=-1)([input_layers, time_embedding])
    
    for i in range(num_layers):
        output = AttentionBlock(name='AttBlock' + str(i), num_heads=2, head_size=128, ff_dim=None, dropout=0.2)(output)
    
    output = Flatten()(output)
    output = mlp(output, [250, 150, 12])
    
    return tf.keras.models.Model(input_layers, output)

if __name__ == "__main__":
    #MuTrans = ModelTrans()
    MuTrans = MuTran()
    #MuTrans.build((1, 12, 12))
    MuTrans.summary()
