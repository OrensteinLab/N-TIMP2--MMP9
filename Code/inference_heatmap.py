import os
import pandas as pd
import numpy as np
import tensorflow as tf
from utils import oneHot
from tensorflow.keras.layers import *



if __name__ == '__main__':

    path = '../Data'
    num_shuffles = 10
    predictions_ala = []
    predictions_no_ala = []
    WTaa = 'SINSVHT'
    aaList = 'ACDEFGHIKLMNPQRSTVWY'
    mutlist = []
    for i in range(len(WTaa)):
        for aa in aaList:
            mutlist.append(WTaa[:i] + aa + WTaa[i+1:])
    df = pd.DataFrame(mutlist)
    seq_one_hot_test = np.reshape(np.array(list(df[0].apply(oneHot)), dtype='int8'), (len(df[0]), 140))
    
    for shuffle in range(num_shuffles):
    
        ala_model = tf.keras.models.load_model(f'ala_all_data_model_final{shuffle}')
        predictions_ala.append(ala_model.predict(seq_one_hot_test))
        no_ala_model = tf.keras.models.load_model(f'no_ala_all_data_model_final{shuffle}')
        predictions_no_ala.append(no_ala_model.predict(seq_one_hot_test))
    
    
    pred_ala = np.asarray(predictions_ala)
    avg_pred_ala = np.mean(pred_ala, axis=0)
    pred_no_ala = np.asarray(predictions_no_ala)
    avg_pred_no_ala = np.mean(pred_no_ala, axis=0)
    
    reg_dataset = pd.DataFrame({'seq': np.array(list(df[0])),
                                'ala_test': avg_pred_ala.reshape(len(avg_pred_ala)),
                                'no_ala_test': avg_pred_no_ala.reshape(len(avg_pred_no_ala))})
    reg_dataset.to_csv(f'{path}/heatmap_pred.csv')



