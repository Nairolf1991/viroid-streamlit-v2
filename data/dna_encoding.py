import pandas as pd
import numpy as np
from tensorflow.keras.preprocessing.sequence import pad_sequences
import math
import itertools
import matplotlib.pyplot as plt

################################################################################################################
#ONE-HOT-ENCODING DNA
################################################################################################################
def onehote(sequence):
    '''Takes in entry a DNA sequence,
    Returns a One-Hot-Encoded matrix'''
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3,"-":4}
    seq2 = [mapping[i] for i in sequence]
    matrix=np.array([[1,0,0,0],
           [0,1,0,0],
            [0,0,1,0],
            [0,0,0,1],
            [0,0,0,0]])
    return matrix[seq2]

def preprocessing_ohe(data,maxlen):
    '''Takes in entry a dataframe containing a column named Sequence, containing the DNA sequences,
    Returns X = array of one-hot-encoded matrixes corresponding'''
    #preprocessing X, One hot encoding, padding
    new_features =[]
    for seq in data['Sequence'] :
        new_features.append(pad_sequences(onehote(seq).T, maxlen=maxlen, padding='post').T)
    X = np.array(new_features)
    print(f"    Shape of feature : {X.shape}")
    return X

################################################################################################################
#ENCODING DNA IN FREQUENCY MATRIXES
################################################################################################################

def list_kmers(sequence,k):
    words=[]
    for i in range(0,len(sequence)-(k-1)):
        words.append(sequence[i:i+k])

    return words

def count_kmers(sequence,k):
    return pd.Series(list_kmers(sequence,k)).value_counts()

def frequencies_kmers(sequence,k):
    counts = count_kmers(sequence,k)
    dictionnary=dict(counts)
    for key,value in dictionnary.items():
        dictionnary[key]=value/(counts.sum())
    return dictionnary

"""
[AAAA, AAAT, AAAC, AAAG, AATA, AATC, AATG, AACA, AACT, AACC, AACG"""

def initialize_word_map(k,map_type):
    alphabet = ['A','T','C','G']
    keywords = [''.join(i) for i in itertools.product(alphabet, repeat = k)]
    if map_type=='chaos':
        dimension = int(math.sqrt(4**k))
        xs=[]
        ys=[]
        for word in keywords:
            x=0.5
            y=0.5
            #C(0,1), A(0,0), T(1,0), G(1,1)
            xC, yC = 0, 1
            xA, yA = 0, 0
            xT, yT = 1, 0
            xG, yG = 1, 1
            for letter in word:
                if letter == 'A':
                    x=(x+xA)/2
                    y=(y+yA)/2
                if letter == 'T':
                    x=(x+xT)/2
                    y=(y+yT)/2
                if letter == 'C':
                    x=(x+xC)/2
                    y=(y+yC)/2
                if letter == 'G':
                    x=(x+xG)/2
                    y=(y+yG)/2
            xs.append(x)
            ys.append(y)

        df_words = pd.DataFrame(keywords,columns=['word'])
        df_words['x']=xs
        df_words['y']=ys
        #Position d'un mot dans la matrice [dimension,dimension]:
        # Pos = [rang_x,rang_y]
        xs_to_rank = {np.unique(xs)[i] : i  for i in range(len(np.unique(xs)))}
        ys_to_rank = {np.unique(ys)[i] : i  for i in range(len(np.unique(ys)))}

        word_map = np.zeros((dimension,dimension))
        word_map = word_map.astype('object')

        for i in range(len(df_words)) :
        #Coordinates (x,y) ==> word_map[15-y][x]
            x=xs_to_rank[df_words.loc[i, "x"]]
            y=ys_to_rank[df_words.loc[i, "y"]]
            word_map[(dimension-1)-y][x]=df_words.loc[i, "word"]

    elif map_type=='dumb':
        word_map = np.array(keywords).reshape(int(math.sqrt(4**k)),int(math.sqrt(4**k)))

    return word_map

def matrix_frequencies(sequence,k,word_map):
    frequencies_dict=frequencies_kmers(sequence,k)
    matrix_frequencies=[]
    for i in range(int(math.sqrt(4**k))):
        row=[]
        for word in word_map[i]:
            if word in frequencies_dict.keys():
                row.append(frequencies_dict[word])
            else:
                row.append(0)
        matrix_frequencies.append(row)
    return np.array(matrix_frequencies).reshape(int(math.sqrt(4**k)),int(math.sqrt(4**k)))

def list_frequency_matrix(dataframe,k,map_type):
    '''
    Renvoie la liste des images (np array de shape : (len(df), racine(4**k), racine(4**k), 1)
    '''
    word_map=initialize_word_map(k=k,map_type=map_type)
    images_list = []
    for sequence in dataframe['Sequence']:
        images_list.append(matrix_frequencies(sequence,k=k,word_map=word_map))
    #images_array = np.array(images_list).reshape(len(images_list),int(math.sqrt(4**K)),int(math.sqrt(4**K)),1)
    return np.array(images_list).reshape(len(images_list),int(math.sqrt(4**k)),int(math.sqrt(4**k)),1)

def dna_to_image_preprocessing(df,k,map_type):
    #clean le dataset
    df = df[df['Sequence'].str.contains('-')==False] # On enlève les lignes avec séquences contenant des "-"
    #transformer dna en images
    X = list_frequency_matrix(dataframe=df,k=k,map_type=map_type)
    return X

def remove_unknown_letters(df):
    df = df[df['Sequence'].str.contains('-')==False] # On enlève les lignes avec séquences contenant des "-"
    return df

#########################################################################################################################
#GRAPHICAL FUNCTIONS
#########################################################################################################################

def plot_image_of_sequence(sequence,k,map_type,cmap='viridis'):
    matrix=matrix_frequencies(sequence,k,map_type)
    plt.imshow(matrix, cmap=cmap)


def print_10_random_viroids_and_viruses(df):
    plots_per_row=10
    k=4
    chaos_word_map=initialize_chaos_word_map(k=k)
    dumb_word_map=initialize_dumb_word_map(k=k)
    fig, axes = plt.subplots(2,plots_per_row, figsize=(15,5))
    for i in range(plots_per_row):
        j=np.random.randint(len(viroids))
        axes[0,i].imshow(matrix_frequencies(viroids['Sequence'][j],k,matrix_words=chaos_word_map))
        axes[1,i].imshow(matrix_frequencies(viroids['Sequence'][j],k,matrix_words=dumb_word_map))

    fig.suptitle(f'{plots_per_row} Random Viroids - Chaos (Up) vs Dumb (Down)', fontsize=32, weight='bold')
    plt.tight_layout()

    fig2, axes2 = plt.subplots(2, plots_per_row, figsize=(15,5))
    for i in range(plots_per_row):
        j=np.random.randint(len(viruses))
        axes2[0,i].imshow(matrix_frequencies(viruses['Sequence'][j],k,matrix_words=chaos_word_map))
        axes2[1,i].imshow(matrix_frequencies(viruses['Sequence'][j],k,matrix_words=dumb_word_map))
    fig2.suptitle(f'{plots_per_row} Random Viruses - Chaos (Up) vs Dumb (Down)', fontsize=32, weight='bold')
    plt.tight_layout()
