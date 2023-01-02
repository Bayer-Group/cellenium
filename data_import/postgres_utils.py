from psycopg2 import extras
from psycopg2.extensions import register_adapter, AsIs
import numpy as np
from sqlalchemy import create_engine


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


def import_df(df, schema_and_table: str):
    """
    Append a dataframe into a postgres table - works with reordered/missing columns
    """

    def _quote_special_char_column(c):
        if '.' in c:
            return f'"{c}"'
        return c

    # Create a list of tuples from the dataframe values
    tuples = [tuple(x) for x in df.to_numpy()]
    # Comma-separated dataframe columns
    cols = ','.join([_quote_special_char_column(str(c)) for c in df.columns])
    # SQL query to execute
    query = "INSERT INTO %s(%s) VALUES %%s" % (schema_and_table, cols)
    with engine.connect() as connection:
        cursor = connection.connection.cursor()
        extras.execute_values(connection.connection.cursor(), query, tuples)
        cursor.close()
        connection.connection.commit()


url = f"postgresql://postgres:postgres@localhost:5001/postgres"
engine = create_engine(url)
