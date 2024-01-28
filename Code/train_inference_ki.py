from scipy.stats.stats import pearsonr
import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"  # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"] = ""
from utils import oneHot
from train_model import train_evaluate_model
import pandas as pd
import numpy as np
from tensorflow.keras.layers import *


kwargs_ala = {'dense_size_1': 8, 'drop_1': 0, 'drop_2': 0.3, 'dense_size_2': 2, 'epochs': 30, 'batch': 2, 'lr': 0.0005}
kwargs_no_ala = {'dense_size_1': 32, 'drop_1': 0.2, 'drop_2': 0.1, 'dense_size_2': 4, 'epochs': 50, 'batch': 2, 'lr': 0.0005}

def mut_to_seq(mut, wt):
    pos = int(mut[1:-1])
    seq = wt[:pos - 1] + mut[-1] + wt[pos:]
    return seq

if __name__ == '__main__':

    path = '../Data'
    WTaa = 'CSCSPVHPQQAFCNADVVIRAKAVSEKEVDSGNDIYGNPIKRIQYEIKQIKMFKGPEKDIEFIYTAPSSAVCGVSLDVGGKKEYLIAGKAEGDGKMHIT'
    num_shuffles = 10

    # read data
    ala_total_data = pd.read_csv(f'{path}/All_variant_ala.csv')
    no_ala_total_data = pd.read_csv(f'{path}/All_variant_no_ala.csv')
    # read ki data
    ki_df = pd.read_csv(f'{path}/KI_VALUES.csv')
    # change A2C format to seq
    ki_df['Variant'][0] = WTaa
    ki_df['Variant'][1:] = ki_df['Variant'][1:].apply(mut_to_seq, wt=WTaa)
    ki_df['seq'] = ki_df['Variant'].str[3:4] + ki_df['Variant'].str[34:35] +ki_df['Variant'].str[37:38] + \
                 ki_df['Variant'].str[67:68] + ki_df['Variant'].str[70:71] + ki_df['Variant'].str[96:97] + ki_df['Variant'].str[98:99]

    # Remove shared data points from data
    ala_total_data = ala_total_data[~ala_total_data['seq'].isin(list(ki_df['seq']))]
    ala_train = ala_total_data
    no_ala_total_data = no_ala_total_data[~no_ala_total_data['seq'].isin(list(ki_df['seq']))]
    no_ala_train = no_ala_total_data

    # Train ala model (without ki variants) and predict
    predictions_ala = train_evaluate_model(num_shuffles, ala_train, ki_df, save_model=0, model_name='ala',
                                           kwargs=kwargs_ala, verbose=0)
    # Train no ala model (without ki variants) and predict
    predictions_no_ala = train_evaluate_model(num_shuffles, no_ala_train, ki_df, save_model=0,
                                              model_name='no_ala', kwargs=kwargs_no_ala, verbose=0)

    pred_ala = np.asarray(predictions_ala)
    avg_pred_ala = np.mean(pred_ala, axis=0)
    pred_no_ala = np.asarray(predictions_no_ala)
    avg_pred_no_ala = np.mean(pred_no_ala, axis=0)
    pred_ER = avg_pred_ala + avg_pred_no_ala
    reg_dataset = pd.DataFrame({'seq': np.array(list(ki_df['seq'])), 'ki_obs': np.array(list(ki_df['log_ki_ratio'])), 'ala_test': avg_pred_ala.reshape(len(avg_pred_ala)),
                                'no_ala_test': avg_pred_no_ala.reshape(len(avg_pred_no_ala)),
                                'ER_pred': pred_ER.reshape(len(pred_ER))})
    reg_dataset.to_csv(f'{path}/variants_ki_prediction.csv')
    r, p = pearsonr(reg_dataset['ER_pred'], reg_dataset['ki_obs'])
    print(f'ki model prediction: R={r}, P-VALUE={p}')


