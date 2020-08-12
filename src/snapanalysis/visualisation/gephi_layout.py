import os
import sh
import tempfile
import networkx as nx
import shutil

from snapanalysis.config import get_logger

_GEPHI_FORCEATLAS_JAR_LOCATION = '/jars/gephi-forceatlas2.jar'

_gephi_forceatlas2 = sh.java.bake("-XX:+UnlockExperimentalVMOptions",
                                 "-XX:+UseCGroupMemoryLimitForHeap",
                                 "-jar",
                                 _GEPHI_FORCEATLAS_JAR_LOCATION)


def gephi_forceatlas2_layout(network, desired_basename,
                             outdir,
                             **kwargs):
    """
    Generates forceatlas2 layout using gephi.

    Uses `gephi-toolkit-forceatlas2-standalone` behind the scenes.
    See https://github.com/lukauskas/gephi-toolkit-forceatlas2-standalone

    Produces a gexf file and a pdf.

    :param network: network to layout
    :param desired_basename: basename to give to the gexf and pdf files
    :param outdir: output directory to store pdf and gexf
    :param kwargs: kwargs passed to gephi-toolkit-forceatlas2-standalone, see README
    :return:
    """
    logger = get_logger('gephi_forceatlas2_layout')

    tmp_dir = tempfile.mkdtemp()
    try:
        tmp_file = os.path.join(tmp_dir, desired_basename + '.gexf')
        nx.write_gexf(network, tmp_file)

        logger.info('Gephi processing {}'.format(desired_basename))
        _gephi_forceatlas2(tmp_file, outdir=outdir, **kwargs)
        logger.info('Gephi processing for {} done'.format(desired_basename))
    finally:
        try:
            shutil.rmtree(tmp_dir)
            del tmp_dir
        except FileNotFoundError:
            pass

