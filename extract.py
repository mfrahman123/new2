#Import Statements
import pandas as pd
import numpy as np
import requests
import gzip
import shutil
import time

def get_important_data():
    #get data from url and store content
        try:
            # Extraction of data for Ensembl, Symbol and Family
            url = "http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF"
            req = requests.get(url)
            url_content = req.content
            tsv_file1 = open('download1.tsv', 'wb')
            tsv_file1.write(url_content)
            tsv_file1.close()

            #Extraction of data for full name and chromosomal location
            # old_url = 'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=md_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit'
            new_url = 'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_pub_chrom_map&col=md_prot_id&col=md_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit'
            req2 = requests.get(new_url)
            url_content2 = req2.content
            tsv_file2 = open('download2.tsv', 'wb')
            tsv_file2.write(url_content2)
            tsv_file2.close()

            return 'Data was extracted successfully!'
        except Exception as e:

            return e



def get_htf_data(inp):
    try:
        '''Extract data regarding gene symbol, Ensembl ID and gene family for Human Transcription Factor Data from Human TF DB.'''

        #create pandas dataframe from tsv file
        htf_df = pd.read_csv('download1.tsv',sep='\t')
        htf_df.rename(columns={'Entrez ID':'Entrez_ID'},inplace=True)

        #Remove unnecessary columns
        htf_df.drop(columns=['Species','Protein','Entrez_ID'],inplace=True)

        #provide index name
        htf_df.index.name = 'ID'

        #locate user inputted TF
        idx_loc = htf_df.loc[htf_df.Symbol == inp]

        symbol_data = idx_loc.Symbol.values[0]
        ensembl_data = idx_loc.Ensembl.values[0]
        family_data = idx_loc.Family.values[0]

        # data_dict = {'Symbol':symbol_data, 'Ensembl': ensembl_data , 'Family':family_data}

        return symbol_data, ensembl_data, family_data

    except Exception as e:

        return e


def get_htf_target_data(tf_gene_symbol):

    '''Gets data for which genes the transcription factors target with some additional info.

    Key Arguments:

    tf_gene_symbol --- Gene Symbol of valid human transcription factor '''

    try:

        #Extraction of data for target genes.
        target_url = 'http://bioinfo.life.hust.edu.cn/hTFtarget/static/hTFtarget/tmp_files/targets/' + tf_gene_symbol + '.target.txt.gz' 
        r2 = requests.get(target_url, stream=True) 
        with open('new.txt.gz', 'wb') as f:
            f.write(r2.content)
            #decompress gz file and convert to txt file
        with gzip.open('new.txt.gz', 'rb') as f_in:
            with open('new.txt', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        #create pandas dataframe from txt file
        target_df = pd.read_csv('new.txt', sep = '\t')

        #provide index name and appropriate column names
        target_df.index.name = 'ID'
        target_df.rename(columns={'TF_name':'Symbol', 'target_id':'Target_Ensembl',
                                  'target_name':'Target_Name', 'target_synonyms':'Target_Synonyms'},inplace=True)

        a = np.ndarray.tolist(target_df.Target_Ensembl.values)
        b = np.ndarray.tolist(target_df.Target_Name.values)
        c = np.ndarray.tolist(target_df.Target_Synonyms.values)

        i = 0
        dict1 = {}

        for x in list(zip(b,a,c)) :

            i += 1
            str_i = str(i)

            dict2 = {'Name' : x[0], 'Ensembl' : x[1] }

            dict1['Target' + str_i] = dict2



            if x[2] == '-':

                dict1['Target' + str_i]['Synonyms_for_target'] = 'None'

            else:

                dict1['Target' + str_i]['Synonyms_for_target'] = x[2].replace(',',' / ')



        return dict1

    except Exception as e:

        return 'No Targets Found'


def get_tf_name_location(inp):

    '''Get data for full name and location of transcription factor'''

    try:

        #create dataframe from tsv file.
        location_df = pd.read_csv('download2.tsv', sep='\t')

        #remove irrelevant data
        location_df.drop(columns=['HGNC ID'],inplace=True)

        #provide index name
        location_df.index.name = 'ID'

        #format column names
        location_df.rename(columns={'Approved symbol':'Symbol', 'Approved name': 'Full_name',
                            'Ensembl ID(supplied by Ensembl)' : 'Ensembl', 'UniProt ID(supplied by UniProt)': 'Uniprot'},inplace=True)

        data = location_df.loc[location_df.Symbol == inp]
        chr_location = data.Chromosome.values[0]
        full_name = data.Full_name.values[0]
        uniprot = data.Uniprot.values[0]

        #Retrieve rows and extract relevant data for subcellular location and function from uniprot database
        subcellular_df = pd.read_csv('uniprot.txt',sep='\t')
        subcellular_df.drop(columns=['Status'],inplace=True)
        subcellular_df.rename(columns={'Entry name':'Entry_name', 'Subcellular location [CC]':'SC', 'Function [CC]':'F'}, inplace=True)
        row_data = subcellular_df.loc[subcellular_df.Entry == uniprot]
        subcell_loc = row_data.SC.values[0]
        func_data = row_data.F.values[0]

        return chr_location, full_name, uniprot, subcell_loc, func_data

    except Exception as e:

        return e


def all_data_to_dict(inp):
    a = get_htf_data(inp)
    b = get_tf_name_location(inp)



    c = {'Symbol':a[0], 'Ensembl': a[1], 'Gene_family': a[2], 'Chr_location': b[0], 'Full_name': b[1], 'Uniprot':b[2],
        'Subcellular_loc': b[3], 'Function':b[4]}

    return c

