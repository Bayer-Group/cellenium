import logging
import os

from postgres_utils import (
    get_aws_db_engine,
    get_local_db_engine,
)

# true if this code is running in AWS Batch instead of locally (e.g. make target). In AWS Batch mode, the user has
# already registered a study in the web UI and the import code here is triggered by the S3 file upload.
IS_AWS_DEPLOYMENT = os.environ.get("AWS") is not None

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s.%(msecs)03d %(process)d %(levelname)s %(name)s:%(lineno)d %(message)s",
    datefmt="%Y%m%d-%H%M%S",
    handlers=[logging.FileHandler("./cellenium_cli.log"), logging.StreamHandler()] if IS_AWS_DEPLOYMENT else [logging.StreamHandler()],
)


def create_db_engine():
    return get_aws_db_engine() if IS_AWS_DEPLOYMENT else get_local_db_engine(os.environ.get("ENGINE_URL"))


engine = create_db_engine()
