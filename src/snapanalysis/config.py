import os
import traceback
from contextlib import contextmanager

import datetime

import sys
import yaml
import logging
import logging.config
import logging.handlers
import warnings
import tables

try:
    RAW_DATA_DIRECTORY = os.environ.get('SNAPANALYSIS_RAW_DATA', 'data/raw')
    EXTERNAL_DATA_DIRECTORY = os.environ.get('SNAPANALYSIS_EXTERNAL_DATA', 'data/external')
    INTERIM_DATA_DIRECTORY = os.environ.get('SNAPANALYSIS_INTERIM_DATA', 'data/interim')
    OUTPUT_DIRECTORY = os.environ.get('SNAPANALYSIS_OUTPUT', 'output/')
except KeyError:
    raise Exception('Please ensure SNAPANALYSIS_RAW_DATA, SNAPANALYSIS_EXTERNAL_DATA, '
                    'SNAPANALYSIS_INTERIM_DATA and SNAPANALYSIS_OUTPUT environment variables are set')


if not os.path.isdir(RAW_DATA_DIRECTORY):
    raise FileNotFoundError('RAW_DATA_DIRECTORY {} does not exist'.format(RAW_DATA_DIRECTORY))

if not os.path.isdir(EXTERNAL_DATA_DIRECTORY):
    raise FileNotFoundError('EXTERNAL_DATA_DIRECTORY {} does not exist'.format(RAW_DATA_DIRECTORY))

if not os.path.isdir(INTERIM_DATA_DIRECTORY):
    raise FileNotFoundError('INTERIM_DATA_DIRECTORY {} does not exist'.format(RAW_DATA_DIRECTORY))

if not os.path.isdir(OUTPUT_DIRECTORY):
    raise FileNotFoundError('OUTPUT_DIRECTORY {} does not exist'.format(RAW_DATA_DIRECTORY))

def get_logger(name):
    return logging.getLogger('snapanalysis.' + name)

def setup_logging(
        default_path='logging.yaml',
        default_level=logging.INFO,
        env_key='SNAPANALYSIS_LOG_CONFIG'):
    """Setup logging configuration
    Code based on https://fangpenlin.com/posts/2012/08/26/good-logging-practice-in-python/
    """
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):
        with open(path, 'rt') as f:
            config = yaml.safe_load(f.read())
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)

    # Disable these natural name warnings and Performance warnings
    warnings.simplefilter('ignore', tables.NaturalNameWarning)
    warnings.simplefilter('ignore', tables.PerformanceWarning)

    logger = get_logger(__name__)


@contextmanager
def timed_segment(message, logger):

    start_time = datetime.datetime.now()
    exception, exc_type, exc_value, exc_traceback = None, None, None, None
    try:
        yield
    except Exception as e:
        exception = e
        exc_type, exc_value, exc_traceback = sys.exc_info()
    finally:
        end_time = datetime.datetime.now()
        diff = (end_time-start_time)
        total_seconds = diff.total_seconds()

        if exception is None:
            status = 'success'
            log_f = logger.info
        else:
            status = 'failure'
            log_f = logger.error

        log_f('Timed segment: {}. Status: {}. Took {:.2f}s'.format(message,
                                                                   status,
                                                                   total_seconds))

        # Re-raise the exception
        if exception is not None:
            raise exception


def ensure_directory_exists(filename, filename_isdir=False):
    if not filename_isdir:
        dir_ = os.path.dirname(filename)
    else:
        dir_ = filename
    if not os.path.isdir(dir_):
        os.makedirs(dir_)


_LOGGING_HAS_BEEN_SETUP = False
if not _LOGGING_HAS_BEEN_SETUP:
    setup_logging()
    _LOGGING_HAS_BEEN_SETUP = True




