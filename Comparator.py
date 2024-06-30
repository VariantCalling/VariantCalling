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

def load_comparator(model_path):
    """
    Wrapper to load the comparator model with self-handled initialization of custom Transformer model
    """
    return tf.keras.models.load_model(model_path, custom_objects={"Transformer": Transformer},safe_mode=False)

def truncate_seq(truncation_range,seq_in):
    truncated_chunks = []
    for truncation_chunk in truncation_range:
        truncated_chunks.append(seq_in[:, truncation_chunk[0]:truncation_chunk[1]])
    for i in range(len(truncated_chunks)):
        print(i)
        if i == 0:
            seq_in = truncated_chunks[i]
            print(np.shape(np.array(seq_in)))
        else:
            seq_in = np.concatenate((seq_in, truncated_chunks[i]), axis=1)
            print(np.shape(np.array(seq_in)))
    return seq_in

def run_comparator_inference(model, pickle_file, truncate_range=[]):
    """
    Wrapper to run inference using comparator model and a user-input pickle file
    """
    pickle_file = pickle_file
    pickle_seq = []

    with (open(pickle_file, "rb")) as openfile:
        while True:
            try:
                pickle_seq.append(pickle.load(openfile))
            except EOFError:
                break
    X_inp_char = []
    for element in pickle_seq[0]:
        X_inp_char.append(list(element[2]))
    X_inp_char = np.array(X_inp_char)    

    transdict = {"A":0, "C": 1, "G":2, "T":3}
    X_inp_val = np.vectorize(transdict.get)(X_inp_char)

    if truncate_range != []:
        X_inp_val = truncate_seq(X_inp_valr)

    clones = []
    with open("clones_lib/" + clone_file, "r") as f:
        for clone in f:
            alignment = []
            for char in clone.strip():
                alignment.append(char)
            clones.append(alignment)
    proportion_tracker = [0 for i in len(clones) + 1]
    clones_int = np.vectorize(transdict.get)(clones)
    pred_list = []
    for i in range(len(clones_int)):
        ref_list = []
        inp_list = []
        for j in range(np.shape(X_inp_val)[0]):
            ref_list.append(clones_int[i])
            inp_list.append(X_inp_val[j].tolist())
        pred_list.append(model.predict([pd.DataFrame(ref_list).values, pd.DataFrame(inp_list).values]))
    return pred_list