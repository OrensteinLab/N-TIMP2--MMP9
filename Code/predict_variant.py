import os
import pandas as pd
import numpy as np
import tensorflow as tf
from utils import oneHot
from tensorflow.keras.layers import *



if __name__ == '__main__':
    num_shuffles = 10
    predictions_ala = []
    predictions_no_ala = []
    WTaa = 'SINSVHT'
    aaList = 'ACDEFGHIKLMNPQRSTVWY'
    seq = input(f"Please insert the variant 7 positions sequence:")
    seq = seq.upper()
    if len(seq)!= 7 or not set(seq).issubset(set(aaList)):
        raise ValueError('Invalid sequence')


    seq_one_hot_test = np.reshape(np.array(oneHot(seq), dtype='int8'),(1,140))
    
    for shuffle in range(num_shuffles):
    
        ala_model = tf.keras.models.load_model(f'ala_all_data_model_final{shuffle}')
        predictions_ala.append(ala_model.predict(seq_one_hot_test, verbose=0))
        no_ala_model = tf.keras.models.load_model(f'no_ala_all_data_model_final{shuffle}')
        predictions_no_ala.append(no_ala_model.predict(seq_one_hot_test, verbose=0))
    
    
    pred_ala = np.asarray(predictions_ala)
    avg_pred_ala = np.mean(pred_ala, axis=0)
    pred_no_ala = np.asarray(predictions_no_ala)
    avg_pred_no_ala = np.mean(pred_no_ala, axis=0)
    print (f'The predicted log2 ER of variant {seq} is {(avg_pred_ala[0] + avg_pred_no_ala[0])[0]}')
