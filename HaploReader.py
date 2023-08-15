"""
HaploReader.py

This file contains the inference running and training utilities for Haplo-Reader
"""

import tensorflow as tf
import numpy as np
from keras.datasets import cifar10
from keras.models import Sequential
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.layers.core import Dense, Dropout, Activation, Flatten, Lambda
from keras.utils.np_utils import to_categorical
import VariantCalling as vc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from sklearn.metrics import r2_score

class HaploReaderTFTrain:
    def __init__(self, alignments, prob_lists, nb_epoch, batch_size, nb_filters, nb_conv, nb_layer_conv, nb_pool, nb_node, nb_layer_ff, train_loop) -> None:
        """Class definition for the training algorithm class.

        Args:
            alignments (np_array): Alignments data that follows the following format: [nb_alignment, nb_channel, img_row, img_col]
            nb_epoch (int): Number of training epoch
            batch_size (int): Mini-batch size 
            nb_filters (int): Kernel size, square
            nb_conv (int): Number of kernels
            nb_layer_conv (int): Number of Conv layers
            nb_pool (int): MaxPool kernel size
            nb_node (int): Number of nodes in FF
            nb_layer_ff (int): Number of feed-forward layer
            train_loop (int): Number of training loop
        """
        self.alignments = alignments
        self.prob_lists = prob_lists
        self.nb_epoch = nb_epoch
        self.batch_size = batch_size
        self.nb_filters = nb_filters
        self.nb_conv = nb_conv
        self.nb_layer_conv = nb_layer_conv
        self.nb_pool = nb_pool
        self.nb_node = nb_node
        self.nb_layer_ff = nb_layer_ff
        self.train_loop = train_loop
        self.img_row = np.shape(self.alignments)[2]
        self.img_col = np.shape(self.alignments)[3]
        self.nb_channel = np.shape(self.alignments)[1]
        self.nb_mutations = np.shape(self.prob_lists)[1]
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
        hfont = {'fontname':'serif'}

    def model_builder(self):
        """Model-builder, generate Keras model based on the parameters
        
        Returns:
            Keras Model: Model generated based on the following architecture
 
            Please note that the model follows the standard architecture of CNN as shown below:
                Input -> (Conv / MaxPool / ReLU) * nb_layer_conv -> (FF / ReLU / DropOut) * nb_layer_ff -> Softmax(nb_clones)
        
        """
        model = Sequential()
        for _ in range(0, self.nb_layer_conv):
            model.add(Conv2D(self.nb_filters, self.nb_conv, self.nb_conv, padding="same"))
            model.add(MaxPooling2D(pool_size = (self.nb_pool, self.nb_pool),padding="same"))
            model.add(Dropout(0.1))
            model.add(Activation('relu'))
        model.add(Flatten())
        for _ in range(0, self.nb_layer_ff):
            model.add(Dense(self.nb_node))
            model.add(Dropout(0.1))
            model.add(Activation('relu'))
        model.add(Dense(self.nb_mutations))
        model.add(Activation('softmax'))
        model.add(Lambda(lambda x: x * 100))
        model.compile(metrics=['mse'],loss='categorical_crossentropy',optimizer='adam')
        return model

    def train_val_split(self):
        """Summary
        """
        rng = np.random.default_rng(seed=42) # use a fixed random generator so runs are consistent
        idxs = np.arange(self.alignments.shape[0])    
        rng.shuffle(idxs)    
        split_idx = int(self.alignments.shape[0]*0.8)
        self.train_alignments, self.valid_alignments = self.alignments[idxs[:split_idx]], self.alignments[idxs[split_idx:]]
        self.train_prob_lists, self.valid_prob_lists = np.array(self.prob_lists)[idxs[:split_idx]], np.array(self.prob_lists)[idxs[split_idx:]]
        self.train_alignments = self.train_alignments.reshape(self.train_alignments.shape[0], self.img_row, self.img_col, self.nb_channel)
        self.valid_alignments = self.valid_alignments.reshape(self.valid_alignments.shape[0], self.img_row, self.img_col, self.nb_channel)

        self.train_alignments = self.train_alignments.astype('float32')
        self.valid_alignments = self.valid_alignments.astype('float32')
        
        # Normalise against 3 (we only have 4 basic nucleotides)
        self.train_alignments /= 3
        self.valid_alignments /= 3

    def model_train(self):
        """Build the model and run the training
        """
        self.model = self.model_builder()
        save_best_model = SaveBestModel()
        self.trained_model = self.model.fit(self.train_alignments, self.train_prob_lists, batch_size = self.batch_size, epochs = self.nb_epoch,  verbose = 1, validation_data = (self.valid_alignments, self.valid_prob_lists),callbacks=[save_best_model])
        self.model.set_weights(save_best_model.best_weights)

    def print_train_trend(self):
        """Print the last model training trend
        """
        print(self.trained_model.history.keys())
        # summarize history for accuracy
        plt.plot(self.trained_model.history['mse'])
        plt.plot(self.trained_model.history['val_mse'])
        plt.title('model mse')
        plt.ylabel('mse')
        plt.xlabel('epoch')
        plt.legend(['train', 'test'], loc='upper left')
        plt.show()
        # summarize history for loss
        plt.plot(self.trained_model.history['loss'])
        plt.plot(self.trained_model.history['val_loss'])
        plt.title('model loss')
        plt.ylabel('loss')
        plt.xlabel('epoch')
        plt.legend(['train', 'test'], loc='upper left')
        plt.show()

    def eval_co_def(self):
        """Calculate coefficient of determination using the validation set
        """
        self.predict_array = self.model.predict(self.valid_alignments).tolist()
        self.df_predict = pd.DataFrame(self.predict_array)
        self.df_actual = pd.DataFrame(self.valid_prob_lists)
        
        # Setup the headers
        self.df_predict.columns = self.clone_names
        self.df_actual.columns = self.clone_names

        plt.figure(figsize=(4 + (8*len(clone_names)), 8), dpi=80)
        for i in range(0,len(self.clone_names)):
            plt.subplot(1, len(self.clone_names), i+1)
            r2 = r2_score(self.df_actual[self.clone_names[i]], self.df_predict[self.clone_names[i]])
            yp_1 = np.polyval([1, 0], self.df_actual[self.clone_names[i]])
            plt.plot(self.df_actual[self.clone_names[i]], yp_1, color='orange')
            plt.title('{seq_name} - R2 = {r2_score}'.format(seq_name=self.clone_names[i],r2_score=round(r2,2)))
            plt.ylabel('Predicted')
            plt.xlabel('Actual')
            plt.scatter(self.df_predict[clone_names[i]], self.df_actual[clone_names[i]],color='red')

    def eval_mse(self):
        """Calculate MSE using the validation set
        
        Returns:
            list: MSE of each clones
        """
        mse = []
        for i in range(0,len(clone_names)):
            mse.append(((df_predict[clone_names[i]] - df_actual[clone_names[i]])**2).mean(axis=ax))
        return mse

class SaveBestModel(tf.keras.callbacks.Callback):

    """Class definition for callbacks to save the best mse score
    """
    
    def __init__(self, save_best_metric='val_mse', this_max=False):
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
