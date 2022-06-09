### imports
import streamlit as st
from PIL import Image
import pandas as pd
import numpy as np
import requests
import data.dna_encoding as enc
import matplotlib.pyplot as plt


# fonction qui import le dataframe
def get_dataframe_data():
        df=pd.read_csv('data/df_final_v2.csv')
        df = df[['Sequence', 'IsViroid', 'Length']]
        return df

### header -  title and subtitle
st.markdown("""# *Bio informatique : Les séquences ARN des viroids vs celles  des virus* """)


### SIDEBAR LOGO Le wagon ####
# title sidebar
with st.sidebar:
    original_title = '<p style="font-family:sans-serif; color:white; font-size: 30px;">Data Science  Viroids Project</p>'
    st.markdown(original_title, unsafe_allow_html=True)

# image wagon
with st.sidebar:
    image = Image.open('images-streamlit/wagon.png')
    st.image(image, caption='Le Wagon', use_column_width=False)


### SIDEBAR BUTTON to navigate through pages ####
with st.sidebar:
    item_select = st.radio('Navigation', ('Transformation', 'Demo'))
    st.write(item_select)


# ### CODE FOR PAGE CONTEXTE ###
# if item_select == 'Introduction':
#     # contexte, definition et images
#     st.markdown("""
#         ## Introduction""")

#     # photo viroid sur plante
#     st.markdown(""" 👉 A quoi ressemble un viroid sur une plante? """)
#     image = Image.open('images-streamlit/viroid-plant.jpeg')
#     st.image(image, caption='viroid-plant', use_column_width=False)

#     # image viroid
#     st.markdown(""" 👉 A quoi ressemble un viroid ?""")
#     image = Image.open('images-streamlit/viroid-image.png')
#     st.image(image, caption='viroid-image', use_column_width=True)

#     # image séquence ARN viroid
#     st.markdown(""" 👉 La traduction d'une séquence ARN sur un viroid?""")
#     image = Image.open('images-streamlit/GCAT viroid - forme.png')
#     st.image(image, caption='séquence ARN sur viroid', use_column_width=True)


### CODE FOR PAGE DATA
# elif item_select == 'Dataset':
#     st.markdown("""## Dataset utilisé : viroids & virus""")

#     # 3 petites infos sur les tableaux
#     col1, col2 = st.columns(2)
#     col1.metric("tableau viroid 👇", "9397 lignes")
#     col1.write("""Longueur entre 200 et 500, produit par des chercheurs en janvier 2022.
#                Condensé de toutes les séquences viroids recensées ce jour.""")
#     col2.metric("tableau virus 👇", "10000 lignes")
#     col2.write("échantillons de longueur 500, dataset produit par le doctorant. Echantillons de séquences ARN de virus.")

#     # partie du dataframe
#     st.markdown("""### Echantillon random de 10 séquences du Dataset 🔎 """)
#     df = get_dataframe_data()
#     st.write(df.sample(n = 10))



### CODE FOR PAGE MODELE ####
if item_select == 'Transformation':
    st.markdown("""## 🖥️ Transformation des séquences en images pour l'algorithme""")
    st.markdown("""👉 On prend une séquence ARN :""")
    # ramenons qq formules ici pour la suite
    random_result = None
    big_list_kmers = None
    word_map = None
    matrix =  None

    # images
    image = Image.open('images-streamlit/arrow.png')

    # bouton "generate sequence" au centre
    col1, col2, col3 = st.columns(3)
    if col2.button('Generate sequence'):
            df = get_dataframe_data()
            random_result = df.sample(1).reset_index()
            st.write(random_result[['Sequence', 'IsViroid']][0])

    # image flèche 1
    col1, col2, col3 = st.columns(3)
    col2.image(image, caption='', width=200)

    # liste de mots
    col1, col2, col3 = st.columns(3)
    col1.markdown(""" 👉 On la coupe en mots de 4 lettres et on obtient la liste de tous les mots : """)
    if random_result is not None:
        big_list_kmers = enc.list_kmers(random_result['Sequence'][0], 4)
        col2.write(", ".join(big_list_kmers[0:10])+ "...")

    # image flèche 2
    col1, col2, col3 = st.columns(3)
    col2.image(image, caption='', width=200)

    # button matrice
    st.markdown(""" 👉 On met tous les mots possibles dans une matrice (GCAT)  :""")

    # on fait une matrice des mots
    if random_result is not None:
        word_map = enc.initialize_word_map(k=4, map_type="dumb")
        st.write(word_map)
    # image flèche
        col1, col2, col3 = st.columns(3)
        col2.image(image, caption='', width=200)

    # on remplace les mots par leur fréquence d'apparition dans la séquence
    st.markdown(""" 👉 On remplace ces mots par leur fréquence d'apparition :""")
    if random_result is not None:
        word_map = enc.initialize_word_map(k=4, map_type="dumb")
        matrix = enc.matrix_frequencies(random_result['Sequence'][0], 4, word_map)
        st.write(matrix)
        # image flèche
        col1, col2, col3 = st.columns(3)
        col2.image(image, caption='', width=200)
        # on fait l'image
        st.markdown(""" 👉 Voilà une représentation graphique de cette matrice :""")
        fig, ax = plt.subplots()
        plt.axis('off')
        im = ax.imshow(matrix, cmap='viridis')
        fig.colorbar(im, ax=ax)
        st.pyplot(fig)


### CODE FOR PAGE DEMO ####
else:
    st.markdown("""## Demo""")
    # streamlit => API
    seq = st.text_input('Upload une sequence ARN', 'GCAT')
    url = 'https://viroid-docker-image-pvgqbd4luq-ew.a.run.app/viroid'
    if seq is not None:
        params = {
            "seq": seq}
        response = requests.get(url, params=params)
        if response.status_code ==200:
            col1, col2 = st.columns(2)
            viroid_resp = round(float(response.json())*100, 2)
            col1.write("""Probability of being a viroid:""")
            col2.write(f'{viroid_resp} %')
            if float(viroid_resp/100) < 0.5:
                col1.markdown("""### Result :""")
                col2.markdown("""### That's not a viroid !""")
            else:
                col1.markdown("""### Result :""")
                col2.markdown("""### That's a viroid !""")
        else:
            st.write("error de Loup: ", response.status_code)

        # seq covid 19
        st.text("")
        st.text("")
        st.text("")
        st.text("")
        st.text("")
        st.text("")
        st.text("")
        col1, col2= st.columns(2)
        col1.write(""" *En bonus, échantillon de la séquence du COVID :*""")
        col1.write("""TGAGTACAGACACTGGTGTTGAACATGTTACCTTCTTCATCTACAATAAAATTGTTGATGAGCCTGAAGAACATGTCCAAATTCACACAATCGACGGTTCATCCGGAGTTGTTAATCCAGTAATGGAACCAATTTATGATGAACCGACGACGACTACTAGCGTGCCTTTGTAAGCACAAGCTGATGAGTACGAACTTATGTACTCATTCGTTTCGGAAGAGACAGGTACGTTAATAGTTAATAGCGTACTTCTTTTTCTTGCTTTCGTGGTATTCTTGCTAGTTACACTAGCCATCCTTACTGCGCTTCGATTGTGTGCGTACTGCTGCAATATTGTTAACGTGAGTCTTGTAAAACCTTCTTTTTACGTTTACTCTCGTGTTAAAAATCTGAATTCTTCTAGAGTTCCTGATCTTCTGGTCTAAACGAACTAAATATTATATTAGTTTTTCTGTTTGGAACTTTAATTTTAGCCATGGCAGATTCCAACGGTACTATT""")
