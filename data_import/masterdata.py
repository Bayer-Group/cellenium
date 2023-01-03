import logging
import os
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from zipfile import ZipFile

import pandas as pd
import requests
import sqlalchemy
import tqdm
from pybiomart import Dataset

logging.basicConfig(level=logging.DEBUG)


# utils
def download(url, fn):
    resp = requests.get(url, stream=True)
    total = int(resp.headers.get('content-length', 0))
    with open(fn, 'wb') as file, tqdm.tqdm(
            desc=fn,
            total=total,
            unit='iB',
            unit_scale=True,
            unit_divisor=1024,
    ) as bar:
        for data in resp.iter_content(chunk_size=1024):
            size = file.write(data)
            bar.update(size)


def walk_down_tree(children):
    result = [nd.id for nd in children]
    for nd in children:
        if len(nd.children):
            result.extend(walk_down_tree(nd.children))
    return result


def get_all_ncit_codes_below(df, ont_code):
    nodes = dict([(id, Node(id)) for id in df.ont_code.drop_duplicates().tolist()])
    hierarchy = df[['ont_code', 'parent_ont_code']].drop_duplicates()
    all_relations = hierarchy.to_dict(orient='records')
    for nd in all_relations:
        try:
            nodes[nd.get('parent_ont_code')].children.append(nodes[nd.get('ont_code')])
        except KeyError:
            pass  # print("parent_ont_code {} does not exist".format(nd.get('parent_ont_code')))
    cur_node = nodes[ont_code]
    all_nodes_below = walk_down_tree(cur_node.children)

    return list(set(all_nodes_below))


@dataclass
class Node:
    id: str
    children: list[str] = field(default_factory=list)


def parse_mesh_ascii_diseases(fn):
    records = list(filter(None, open(fn).read().split('*NEWRECORD')))
    collect = []
    for record in records:
        record = record.replace('PRINT ENTRY = ', 'ENTRY = ')
        rows = [(_[0], _[2]) for _ in [_.partition(' = ') for _ in list(filter(None, record.split('\n')))]]
        record_dict = defaultdict(list)
        for key, value in rows:
            record_dict[key].append(value)
        collect.append(record_dict)
    df = pd.DataFrame(collect).explode(['UI']).explode(['MN'])
    df = df.dropna(subset=['MN'])
    df = df.loc[df.MN.str.startswith('C'), :]  # only those are in fact a disease
    df['ontid'] = 1
    return df


def get_gene_mappings(dataset_name, attributes, tax_id):
    dataset = Dataset(name=dataset_name,
                      host='http://www.ensembl.org',
                      use_cache=True)
    df = dataset.query(attributes=attributes)
    df = df.dropna().rename(columns={'Gene description': 'display_name',
                                     'HGNC symbol': 'hgnc_symbols',
                                     'Gene name': 'hgnc_symbols',
                                     'Gene stable ID': 'ensembl_gene_id',
                                     'NCBI gene (formerly Entrezgene) ID': 'entrez_gene_ids'})
    df.display_name = df.display_name.str.split('\[Source', expand=True)[0].str.strip()
    df.entrez_gene_ids = df.entrez_gene_ids.astype(int).astype(str)
    df = df.groupby(['ensembl_gene_id', 'display_name']).agg(set).reset_index()
    df.hgnc_symbols = df.hgnc_symbols.apply(list)
    df.entrez_gene_ids = df.entrez_gene_ids.apply(list)
    df['tax_id'] = tax_id
    df['display_symbol'] = df.hgnc_symbols.str.join(';')
    df['omics_type'] = 'gene'
    return df


def parse_ncit(fn):
    df = pd.read_csv(fn)
    df = df.loc[pd.isnull(df.Concept_Status), :]
    df = df.query('Obsolete == False')
    df['label'] = df['Preferred Label']
    df['ont_code'] = df['Class ID'].str.split('#', expand=True)[1]
    df['Parents'] = df.Parents.str.split('|')
    df = df.explode(['Parents'])
    df['parent_ont_code'] = df.Parents.str.split('#', expand=True)[1]
    df['synonym'] = df.Synonyms.str.split('|')

    nodes_below = get_all_ncit_codes_below(df, 'C12219')
    df = df.loc[(df.ont_code.isin(['C12219', *nodes_below]) & df.parent_ont_code.isin(
        ['Thing', 'C12219', *nodes_below])), :]
    return df


def parse_ncbi_taxonomy(filesdir):
    tax = pd.read_csv(filesdir / 'nodes.dmp', sep='\t\|\t', header=None, engine='python')
    names = pd.read_csv(filesdir / 'names.dmp', sep='\t\|\t', header=None, engine='python')
    names[3] = names[3].str.split('\t', expand=True)[0]
    pref = names.loc[names[3] == 'scientific name', [0, 1]].reset_index(drop=True).rename(
        columns={0: 'ont_code', 1: 'label'})
    pref.ont_code = pref.ont_code.astype(str)

    synonyms = names.loc[names[3] != 'scientific name', [0, 1]].drop_duplicates().reset_index(drop=True)
    synonyms = synonyms.groupby(0).apply(
        lambda x: pd.Series({'synonym': x[1].tolist()})).reset_index().rename(columns={0: 'ont_code'})
    synonyms.ont_code = synonyms.ont_code.astype(str)

    tax = tax[[0, 1]].rename(columns={0: 'ont_code', 1: 'parent_ont_code'})
    tax = tax.applymap(str)

    tax_joined = tax.merge(pref.merge(synonyms, on='ont_code'))
    return tax_joined


def get_antibody_mappings_from_biolegend():
    tables = pd.read_html('https://www.biolegend.com/en-us/totalseq/barcode-lookup', encoding="utf-8")
    ab = pd.concat([df[['Description', 'Ensembl Gene Id']] for df in tables if
                    len(set(['Description', 'Ensembl Gene Id']).intersection(df.columns)) == 2], axis=0).rename(
        columns={'Description': 'antibody_symbol', 'Ensembl Gene Id': 'ensembl_gene_id'})
    ab['antibody_name'] = ab['antibody_symbol']
    ab['display_name'] = ab['antibody_symbol']
    ab['display_symbol'] = ab['antibody_symbol']
    ab.ensembl_gene_id = ab.ensembl_gene_id.str.split(',')
    ab = ab.explode('ensembl_gene_id')
    ab.ensembl_gene_id = ab.ensembl_gene_id.str.strip()
    ab.ensembl_gene_id = ab.ensembl_gene_id.str.replace(':', '')
    ab = ab.dropna()
    ab.ensembl_gene_id = ab.ensembl_gene_id.str.replace(r'^NSG', 'ENSG')
    ab = ab.drop_duplicates().reset_index(drop=True)
    ab['omics_type'] = 'CITE-seq'
    ab.loc[ab.ensembl_gene_id.str.startswith('ENSMUSG'), 'tax_id'] = 10090
    ab.loc[ab.ensembl_gene_id.str.startswith('ENSRNOG'), 'tax_id'] = 10116
    return ab


class Dataimport(object):
    def __init__(self, db_config):
        self.db_config = db_config
        self.conn = None
        self.engine = None
        self.cursor = None
        self.meshfn = '../scratch/d2023.txt'
        self.ncitfn = '../scratch/ncit.csv.gz'
        self.ncbitaxonomyfn = '../scratch/taxdmp.zip'

    def __enter__(self):
        self.engine = sqlalchemy.create_engine(
            "postgresql+psycopg2://{user}:{password}@{host}:{port}/{dbname}".format(**self.db_config))
        return self

    def __exit__(self, exception_type, exception_value, trace):
        try:
            self.engine.dispose()
        except AttributeError:
            logging.error('Cannot close DB connection')
            return True  # exception handled successfully

    def add_ontology(self, name, ontid):
        df = pd.DataFrame(pd.DataFrame({'ontid': [ontid], 'name': [name]}))
        df.to_sql('ontology', con=self.engine, index=False, if_exists='append')

    def import_ncbi_taxonomy(self):
        logging.info('importing NCBI taxonomy')
        # download
        url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip'
        if not os.path.exists(self.ncbitaxonomyfn):
            download(url, self.ncbitaxonomyfn)
            with ZipFile(self.ncbitaxonomyfn, 'r') as zipObject:
                for fn in ['names.dmp', 'nodes.dmp']:
                    zipObject.extract(fn, Path(self.ncbitaxonomyfn).parent)

        # make entry into ontology table
        try:
            self.add_ontology(name='NCBI_taxonomy', ontid=3)
        except sqlalchemy.exc.IntegrityError:
            logging.warning('NCBI taxonomy already imported into ontology table')

        # parse and prepare the datatable
        logging.info('parse the NCBI taxonomy files')
        df = parse_ncbi_taxonomy(Path(self.ncbitaxonomyfn).parent)

        # fill concept table
        logging.info('import the NCBI taxonomy concepts')
        concept = df[['ont_code', 'label']].drop_duplicates()
        concept['ontid'] = 3
        concept.to_sql('concept', if_exists='append', index=False, con=self.engine)
        # get ids
        cids = pd.read_sql('SELECT cid, ont_code from concept where ontid = 3', con=self.engine)

        # construct the synonym table
        logging.info('import the NCBI taxonomy concept_synonyms')
        synonyms = df[['ont_code', 'synonym']].explode('synonym').drop_duplicates()
        synonyms = synonyms.merge(cids, on='ont_code')
        synonyms[['cid', 'synonym']].to_sql('concept_synonym', if_exists='append', con=self.engine, index=False)

        # construct the hierarchy table
        logging.info('import the NCIT concept_hierarchy')
        hierarchy = df[['ont_code', 'parent_ont_code']].drop_duplicates()
        hierarchy = hierarchy.merge(cids, left_on='parent_ont_code', right_on='ont_code').rename(
            columns={'cid': 'parent_cid'}).drop(['ont_code_y'], axis=1).rename(columns={'ont_code_x': 'ont_code'})

        hierarchy = hierarchy.merge(cids, on='ont_code')

        hierarchy = hierarchy.dropna(subset=['parent_cid'])
        hierarchy.parent_cid = hierarchy.parent_cid.astype(int)
        hierarchy = hierarchy[['cid', 'parent_cid']].drop_duplicates()
        hierarchy.to_sql('concept_hierarchy', if_exists='append', index=False,
                         con=self.engine)

    def import_ncit(self):
        logging.info('importing NCIT')
        # download
        url = 'https://data.bioontology.org/ontologies/NCIT/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=csv'
        if not os.path.exists(self.ncitfn):
            download(url, self.ncitfn)
        # make entry into ontology table
        try:
            self.add_ontology(name='NCIT', ontid=2)
        except sqlalchemy.exc.IntegrityError:
            logging.warning('NCIT already imported into ontology table')

        # parse and prepare the datatable
        df = parse_ncit(self.ncitfn)

        # fill concept table
        logging.info('import the NCIT concepts')
        concept = df[['ont_code', 'label']].drop_duplicates()
        concept['ontid'] = 2
        concept.to_sql('concept', if_exists='append', index=False, con=self.engine)
        # get ids
        cids = pd.read_sql('SELECT cid, ont_code from concept where ontid = 2', con=self.engine)

        # construct the synonym table
        logging.info('import the NCIT concept_synonyms')
        synonyms = df[['ont_code', 'synonym']].explode('synonym').drop_duplicates()
        synonyms = synonyms.merge(cids, on='ont_code')
        synonyms[['cid', 'synonym']].to_sql('concept_synonym', if_exists='append', con=self.engine, index=False)

        # construct the hierarchy table
        logging.info('import the NCIT concept_hierarchy')
        hierarchy = df[['ont_code', 'parent_ont_code']].drop_duplicates()
        hierarchy = hierarchy.merge(cids, left_on='parent_ont_code', right_on='ont_code').rename(
            columns={'cid': 'parent_cid'}).drop(['ont_code_y'], axis=1).rename(columns={'ont_code_x': 'ont_code'})

        hierarchy = hierarchy.merge(cids, on='ont_code')

        hierarchy = hierarchy.dropna(subset=['parent_cid'])
        hierarchy.parent_cid = hierarchy.parent_cid.astype(int)
        hierarchy = hierarchy[['cid', 'parent_cid']].drop_duplicates()
        hierarchy.to_sql('concept_hierarchy', if_exists='append', index=False,
                         con=self.engine)

    def import_mesh(self):
        logging.info('importing MeSH')
        # download
        url = 'https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/d2023.bin'
        if not os.path.exists(self.meshfn):
            download(url, self.meshfn)

        # make entry into ontology table
        try:
            self.add_ontology(name='MeSH', ontid=1)
        except sqlalchemy.exc.IntegrityError:
            logging.warning('MeSH already imported into ontology table')

        # parse ASCII format
        df = parse_mesh_ascii_diseases(self.meshfn)

        # fill concept table
        # self.engine.execute('DELETE FROM concept WHERE ontid = 1')
        logging.info('import the MeSH concepts')
        concept = df.explode(['MH'])[['ontid', 'UI', 'MH']].rename(
            columns={'MH': 'label', 'UI': 'ont_code'}).drop_duplicates().reset_index(drop=True)
        concept = concept[['ontid', 'ont_code', 'label']]
        concept.to_sql('concept', if_exists='append', index=False, con=self.engine)

        # get ids
        cids = pd.read_sql('SELECT cid, ont_code from concept where ontid = 1', con=self.engine)

        # construct the synonym table
        logging.info('import the MeSH concept_synonyms')
        synonyms = df[['UI', 'ENTRY']].rename(columns={'UI': 'ont_code', 'ENTRY': 'synonym'}).merge(cids, on='ont_code')
        synonyms = synonyms.explode('synonym').drop_duplicates()
        synonyms.synonym = synonyms.synonym.str.split('|', expand=True)[0]
        synonyms[['cid', 'synonym']].to_sql('concept_synonym', if_exists='append', con=self.engine, index=False)

        # construct the hierarchy table
        logging.info('import the MeSH concept_hierarchy')
        hierarchy = df[['UI', 'MN']].explode('MN').rename(columns={'UI': 'ont_code', 'ENTRY': 'synonym'}) \
            .merge(cids, on='ont_code').drop_duplicates()
        hierarchy['MN_PARENT'] = hierarchy['MN'].str.split('.').str[:-1].str.join('.')

        hierarchy = hierarchy.merge(
            hierarchy.rename(columns={'cid': 'parent_cid'})[['parent_cid', 'MN']].drop_duplicates(),
            left_on='MN_PARENT', right_on='MN', how='left').rename(columns={})
        hierarchy = hierarchy.dropna(subset=['parent_cid'])
        hierarchy.parent_cid = hierarchy.parent_cid.astype(int)
        hierarchy = hierarchy[['cid', 'parent_cid']].drop_duplicates()
        hierarchy.to_sql('concept_hierarchy', if_exists='append', index=False,
                         con=self.engine)

    def import_gene_mappings(self):
        # impart gene mappings
        logging.info('importing gene mappings')
        human = get_gene_mappings('hsapiens_gene_ensembl',
                                  ['hgnc_symbol', 'ensembl_gene_id', 'entrezgene_id', 'description'],
                                  9606)
        mus = get_gene_mappings('mmusculus_gene_ensembl',
                                ['ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'],
                                10090)
        rat = get_gene_mappings('rnorvegicus_gene_ensembl',
                                ['ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'],
                                10116)
        genes = pd.concat([human, mus, rat]).reset_index(drop=True)
        genes[['omics_type', 'tax_id', 'display_symbol', 'display_name']].to_sql('omics_base',
                                                                                 if_exists='append', index=False,
                                                                                 con=self.engine)
        omics_ids = pd.read_sql("SELECT omics_id from omics_base where omics_type = 'gene' order by omics_id",
                                con=self.engine)
        genes['gene_id'] = omics_ids.omics_id
        genes[['gene_id', 'ensembl_gene_id', 'entrez_gene_ids', 'hgnc_symbols']].to_sql('omics_gene',
                                                                                        if_exists='append', index=False,
                                                                                        con=self.engine)

        # TODO CITE-Seq antibodies
        # logging.info('importing CITE-Seq antibodies')
        # ab = get_antibody_mappings_from_biolegend()
        # ab.to_sql('omics', if_exists='append', index=False,
        #           con=self.engine)
        # ...

        # TODO jasper

    def import_masterdata(self):
        self.import_mesh()
        self.import_ncit()
        self.import_ncbi_taxonomy()
        self.import_gene_mappings()


if __name__ == '__main__':
    db_config = {
        "host": "localhost",
        "port": 5001,
        "user": "postgres",
        "password": "postgres",
        "dbname": "postgres"}
    with Dataimport(db_config) as db:
        db.import_masterdata()