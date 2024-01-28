"""
TransformerComparator.py

Constructor of a Transformer-based comparator for long-sequence read genetic data.

The model is a generic model inspired by both Siamese network and the original paper on Transformer (Attention is all you need).

Most of the variables in this file are hardcoded, consider convert to input arguments for flexibility
"""
# Import keras and other libraries
import keras
from keras.layers import Input, Embedding, Dense, Lambda, LayerNormalization, Dropout, Conv2D, MaxPooling2D, Activation, Flatten
from keras.models import Model
from keras.optimizers import Adam
import numpy as np
import tensorflow as tf

class Transformer(tf.keras.Model):
  def __init__(self, *, num_layers, d_model, num_heads, dff,
               input_vocab_size, target_vocab_size, dropout_rate=0.1):
    super().__init__()
    self.encoder = Encoder(num_layers=num_layers, d_model=d_model,
                           num_heads=num_heads, dff=dff,
                           vocab_size=input_vocab_size,
                           dropout_rate=dropout_rate)

    self.decoder = Decoder(num_layers=num_layers, d_model=d_model,
                           num_heads=num_heads, dff=dff,
                           vocab_size=target_vocab_size,
                           dropout_rate=dropout_rate)

    self.final_layer = tf.keras.layers.Dense(target_vocab_size)

  def call(self, inputs):
    # To use a Keras model with `.fit` you must pass all your inputs in the
    # first argument.
    context, x  = inputs

    context = self.encoder(context)  # (batch_size, context_len, d_model)

    x = self.decoder(x, context)  # (batch_size, target_len, d_model)

    # Final linear layer output.
    logits = self.final_layer(x)  # (batch_size, target_len, target_vocab_size)

    try:
      # Drop the keras mask, so it doesn't scale the losses/metrics.
      # b/250038731
      del logits._keras_mask
    except AttributeError:
      pass

    # Return the final output and the attention weights.
    return logits

# Define a custom Transformer layer
class Transformer(keras.layers.Layer):
    def __init__(self, num_heads, dim_feedforward, dropout_rate,embed_dim, name, **kwargs):
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

    def get_config(self):
        # Here we generate the config for loading
        config = super().get_config()
        config.update({
            'num_heads': self.num_heads,
            'dim_feedforward': self.dim_feedforward,
            'dropout_rate': self.dropout_rate,
            'name': self.name
            })

class SaveBestModel(tf.keras.callbacks.Callback):
    """
    Keras Callback to save the best performing model
    """
    def __init__(self, save_best_metric='val_accuracy', this_max=False):
        self.save_best_metric = save_best_metric
        self.max = this_max
        if this_max:
            self.best = float('-inf')
        else:
            self.best = float('inf')

    def on_epoch_end(self, epoch, logs=None):
        metric_value = logs[self.save_best_metric]
        if self.max:
            if metric_value > self.best:
                self.best = metric_value
                self.best_weights = self.model.get_weights()

        else:
            if metric_value < self.best:
                self.best = metric_value
                self.best_weights= self.model.get_weights()

def transformer_model():
    """
    Returns
    -------
    Keras Model :
        Keras model of the Transformer-based comparator
    """
    # Define the vocabulary size and the embedding dimension
    vocab_size = 4 # A, T, C, G
    embed_dim = 16

    # Define the input sequences
    ref_input = Input(shape=(178,), dtype='int32', name='ref_input')
    inp_input = Input(shape=(178,), dtype='int32', name='inp_input')    

    # Embed the input sequences
    ref_embed = Embedding(vocab_size, embed_dim, name='ref_embed')(ref_input)
    inp_embed = Embedding(vocab_size, embed_dim, name='inp_embed')(inp_input)

    # Apply the transformer layer to encode the input sequences
    transformer = Transformer(num_heads=4, dim_feedforward=16, dropout_rate=0.1,embed_dim=embed_dim, name='transformer')
    ref_encoded = transformer(ref_embed)
    inp_encoded = transformer(inp_embed)
    # Compute the cosine similarity between the encoded sequences
    cosine_similarity = Lambda(lambda x: tf.keras.losses.cosine_similarity(x[0], x[1], axis=-1), name='cosine_similarity')([ref_encoded, inp_encoded])    

    # Add some FF laters    
    deep1 = Dense(32, activation='relu', name='deep1')(cosine_similarity)
    #deep2 = Dense(32, activation='relu', name='deep2')(deep1)
    #deep3 = Dense(32, activation='relu', name='deep3')(deep2)  
    #deep4 = Dense(32, activation='relu', name='deep4')(deep3)    

    # Define the output layer
    output = Dense(1, activation='sigmoid', name='output')(deep1)

    # Define the model
    model = Model(inputs=(ref_input, inp_input), outputs=output, name='comparator')    

    # Compile the model
    model.compile(optimizer=Adam(learning_rate=0.002), loss='binary_crossentropy', metrics=['accuracy'])
    return model