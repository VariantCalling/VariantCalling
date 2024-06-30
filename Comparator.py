# Import keras and other libraries
import keras
from keras.layers import Input, Embedding, Dense, Lambda, LayerNormalization, Dropout
from keras.models import Model
from keras.optimizers import Adam
import numpy as np
import tensorflow as tf
import pandas as pd
import random
import PickleUtil as PU
import VariantCalling as vc

# Define the vocabulary size and the embedding dimension
vocab_size = 4 # A, T, C, G
embed_dim = 128

# Define a custom Transformer layer
class Transformer(keras.layers.Layer):
    def __init__(self, num_heads, dim_feedforward, dropout_rate, name, **kwargs):
        super(Transformer, self).__init__(name=name, **kwargs)
        self.num_heads = num_heads
        self.dim_feedforward = dim_feedforward
        self.dropout_rate = dropout_rate
        self.attention = keras.layers.MultiHeadAttention(num_heads=num_heads, key_dim=embed_dim, dropout=dropout_rate, name='attention')
        self.attention_dropout = keras.layers.Dropout(dropout_rate, name='attention_dropout')
        self.attention_norm = keras.layers.LayerNormalization(epsilon=1e-6, name='attention_norm')
        self.ffn = keras.Sequential([
            keras.layers.Dense(dim_feedforward, activation='relu', name='ffn_1'),
            keras.layers.Dense(embed_dim, name='ffn_2')
        ], name='ffn')
        self.ffn_dropout = keras.layers.Dropout(dropout_rate, name='ffn_dropout')
        self.ffn_norm = keras.layers.LayerNormalization(epsilon=1e-6, name='ffn_norm')

    def call(self, inputs, training):
        # Apply the self-attention layer
        attention_output = self.attention(inputs, inputs, training=training)
        # Apply the dropout and residual connection
        attention_output = self.attention_dropout(attention_output, training=training)
        attention_output = self.attention_norm(inputs + attention_output)
        # Apply the feed-forward network
        ffn_output = self.ffn(attention_output)
        # Apply the dropout and residual connection
        ffn_output = self.ffn_dropout(ffn_output, training=training)
        ffn_output = self.ffn_norm(attention_output + ffn_output)
        return ffn_output