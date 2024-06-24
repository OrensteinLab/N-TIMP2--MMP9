import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"] = ""
import pandas as pd
import numpy as np
from utils import oneHot
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import *
from tensorflow import random
import time
import sys


# Models parameters
num_shuffles = 10
minRepetitionThrAla = 100
minRepetitionThrNoAla = 40

kwargs_ala = {'dense_size_1': 8, 'drop_1': 0, 'drop_2': 0.3, 'dense_size_2': 2, 'epochs': 30, 'batch': 2, 'lr': 0.0005}
kwargs_no_ala = {'dense_size_1': 32, 'drop_1': 0.2, 'drop_2': 0.1, 'dense_size_2': 4, 'epochs': 50, 'batch': 2, 'lr': 0.0005}

def split_data(ala_total_data, minRepetitionThrAla, no_ala_total_data, minRepetitionThrNoAla, split):
    '''
    Create train, validation and test set.
    Test, validation sets - keep only data with more than minRepetitionThr - 10% each
    rain set - remove test and validation sets from data
    :param ala_total_data: all ala variants data frame
    :param no_ala_total_data: all no ala variants data frame
    :param split: indicates the user preference for data partitioning and the desired prediction task
    :return: train, validation and test set.
    '''

    ala_test = ala_total_data[ala_total_data['read_count'] >= minRepetitionThrAla]
    rest, ala_test = train_test_split(ala_test, test_size=0.303, random_state=1)
    _, ala_valid = train_test_split(rest, test_size=0.435, random_state=1)

    no_ala_test = no_ala_total_data[no_ala_total_data['read_count'] >= minRepetitionThrNoAla]
    rest, no_ala_test = train_test_split(no_ala_test, test_size=0.361, random_state=1)
    _, no_ala_valid = train_test_split(rest, test_size=0.566, random_state=1)

    if split == '1':
        ala_train = ala_total_data
        no_ala_train = no_ala_total_data
    elif split == '2':
        ala_train = ala_total_data.drop(ala_test.index.union(ala_valid.index), axis=0)
        no_ala_train = no_ala_total_data.drop(no_ala_test.index.union(no_ala_valid.index), axis=0)
    elif split == '3':
        ala_train = ala_total_data.drop(ala_test.index, axis=0)
        no_ala_train = no_ala_total_data.drop(no_ala_test.index, axis=0)
    return ala_train, ala_valid, ala_test, no_ala_train, no_ala_valid, no_ala_test

def train_evaluate_model(num_shuffles, data_train, data_test, save_model, model_name, kwargs, verbose=0):
    '''
    train a model
    :param num_shuffles: int, number of different shuffled models
    :param data_train: data frame, training variants data
    :param data_test: data frame, test variants data
    :param save_model: int, '1'- save the model, '0'- do not save
    :param model_name: str, when save_model=1, name of the model
    :param kwargs: dict, all hyperparameters of the model
    :param verbose: int, '1'- print, otherwise- do not print
    :return: list of predictions on the test set
    '''

    predictions = []
    # One hot encoding- inputs
    train_one_hot = np.reshape(np.array(list(data_train['seq'].apply(oneHot)), dtype='int8'), (len(data_train['seq']), 140))
    test_one_hot = np.reshape(np.array(list(data_test['seq'].apply(oneHot)), dtype='int8'), (len(data_test['seq']), 140))

    # Labels and weights
    labels_ER_train = np.array(list(data_train['log2_ER']))
    weights_train = np.array(list(np.log2(data_train['log2_read_count'])))

    for shuffle in range(num_shuffles):
        np.random.seed(shuffle)
        idx = np.random.permutation(len(train_one_hot))
        train_inputs_shuff = train_one_hot[idx]
        ER_train_shuff = labels_ER_train[idx]
        weights_shuff = weights_train[idx]

        random.set_seed(shuffle)
        model = Sequential()
        model.add(Dense(kwargs['dense_size_1'], activation='relu'))
        model.add(Dropout(kwargs['drop_1']))
        model.add(Dense(kwargs['dense_size_2'], activation='relu'))
        model.add(Dropout(kwargs['drop_2']))
        model.add(Dense(1, activation='linear'))
        model.compile(optimizer=tf.optimizers.Adam(learning_rate=kwargs['lr']), loss=tf.keras.losses.mean_squared_error)
        model.fit(train_inputs_shuff, ER_train_shuff, epochs=kwargs['epochs'], batch_size=kwargs['batch'],
                      sample_weight=weights_shuff, verbose=verbose)
        if save_model == '1':
            model.save(f'{model_name}_all_data_model_final{shuffle}')
        else:
            predictions.append(model.predict(test_one_hot))

    return predictions


if __name__ == '__main__':
    path = '../Data'

    split = input("Please indicate your preference for data partitioning and the desired prediction task. "
                  "You have the following options:\n"
                  "'1': Use all the data to train the model without making any predictions.\n"
                  "'2': Split the data into three sets: a training set (80%), a validation set (10%), and a test set (10%). "
                  "Train the model using the training set and make predictions on the validation set.\n"
                  "'3': Split the data into a training set (90%) and a test set (10%). Train the model using the "
                  "training set and make predictions on the validation set.\n>")


    # Import datasets
    ala_total_data = pd.read_csv(f'{path}/All_variant_ala.csv')
    no_ala_total_data = pd.read_csv(f'{path}/All_variant_no_ala.csv')

    # Create train, validation and test set.
    ala_train, ala_valid, ala_test, no_ala_train, no_ala_valid, no_ala_test = \
        split_data(ala_total_data, minRepetitionThrAla, no_ala_total_data, minRepetitionThrNoAla, split)

    # Train ala model and predict
    predictions_ala = train_evaluate_model(num_shuffles, ala_train, ala_test, save_model=split, model_name='ala', kwargs=kwargs_ala, verbose=1)
    # Train no ala model and predict
    predictions_no_ala = train_evaluate_model(num_shuffles, no_ala_train, no_ala_test, save_model=split, model_name='no_ala', kwargs=kwargs_no_ala, verbose=1)

    if split == '1':
        sys.exit()

    pred_ala = np.asarray(predictions_ala)
    avg_pred_ala = np.mean(pred_ala, axis=0)
    reg_dataset_ala = pd.DataFrame(
        {'seq': np.array(list(ala_test['seq'])), 'obs_test': np.array(list(ala_test['log2_ER'])), 'pred_test': avg_pred_ala.reshape(len(avg_pred_ala))})
    reg_dataset_ala.to_csv(f'{path}/ala_predData_train90_test10.csv')
    r_ala, p_ala = pearsonr(reg_dataset_ala['obs_test'], reg_dataset_ala['pred_test'])
    print(f'ala prediction: R={r_ala}, P-VALUE={p_ala}')

    pred_no_ala = np.asarray(predictions_no_ala)
    avg_pred_no_ala = np.mean(pred_no_ala, axis=0)
    reg_dataset_no_ala = pd.DataFrame(
        {'seq': np.array(list(no_ala_test['seq'])), 'obs_test': np.array(list(no_ala_test['log2_ER'])), 'pred_test': avg_pred_no_ala.reshape(len(avg_pred_no_ala))})
    reg_dataset_no_ala.to_csv(f'{path}/no_ala_predData_train90_test10.csv')
    r_no_ala, p_no_ala = pearsonr(reg_dataset_no_ala['obs_test'], reg_dataset_no_ala['pred_test'])
    print(f'no ala prediction: R={r_no_ala}, P-VALUE={p_no_ala}')



