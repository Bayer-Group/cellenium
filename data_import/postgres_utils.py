import os

import boto3
from psycopg2 import extras
from psycopg2.extensions import register_adapter, AsIs
import numpy as np
from sqlalchemy import create_engine
import json


def addapt_numpy_float64(numpy_float64):
    return AsIs(numpy_float64)


def addapt_numpy_int64(numpy_int64):
    return AsIs(numpy_int64)


def addapt_numpy_float32(numpy_float32):
    return AsIs(numpy_float32)


def addapt_numpy_int32(numpy_int32):
    return AsIs(numpy_int32)


def addapt_numpy_uint32(numpy_uint32):
    return AsIs(numpy_uint32)


def addapt_numpy_array(numpy_array):
    return AsIs(tuple(numpy_array))


register_adapter(np.float64, addapt_numpy_float64)
register_adapter(np.int64, addapt_numpy_int64)
register_adapter(np.float32, addapt_numpy_float32)
register_adapter(np.int32, addapt_numpy_int32)
register_adapter(np.uint32, addapt_numpy_uint32)
register_adapter(np.ndarray, addapt_numpy_array)


class NumpyEncoder(json.JSONEncoder):
    """Special json encoder for numpy types"""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return obj.item()
        return json.JSONEncoder.default(self, obj)


def list_to_pgarray(l):
    return "{" + ",".join(l) + "}"


def import_df(df, schema_and_table: str, engine):
    """
    Append a dataframe into a postgres table - works with reordered/missing columns
    """

    def _quote_special_char_column(c):
        if "." in c:
            return f'"{c}"'
        return c

    # Create a list of tuples from the dataframe values
    tuples = [tuple(x) for x in df.to_numpy()]
    # Comma-separated dataframe columns
    cols = ",".join([_quote_special_char_column(str(c)) for c in df.columns])
    # SQL query to execute
    query = "INSERT INTO %s(%s) VALUES %%s" % (schema_and_table, cols)
    with engine.connect() as connection:
        cursor = connection.connection.cursor()
        extras.execute_values(connection.connection.cursor(), query, tuples)
        cursor.close()
        connection.connection.commit()


def get_aws_db_engine():
    session = boto3.session.Session()
    client = session.client(
        service_name="secretsmanager",
        region_name=os.environ.get("AWS_REGION", "eu-central-1"),
    )

    secret_values = json.loads(
        client.get_secret_value(SecretId=os.environ.get("AWS_DB_SECRET"))[
            "SecretString"
        ]
    )
    db = secret_values["DB"]
    db_host = secret_values["IP"]
    db_port = secret_values["PORT"]
    db_user = secret_values["USER"]
    db_password = secret_values["PASSWORD"]

    return create_engine(
        url=f"postgresql://{db_user}:{db_password}@{db_host}:{db_port}/{db}"
    )


def get_local_db_engine(url: str):
    if url is None:
        url = f"postgresql://postgres:postgres@localhost:5001/postgres"
    return create_engine(url=url)
