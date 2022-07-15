import logging

import cmdstanpy


def setup_loggers(logfile):
    birdman_logger = logging.getLogger("birdman")
    birdman_logger.setLevel(logging.INFO)
    fh = logging.FileHandler(logfile, mode="w")
    sh = logging.StreamHandler()
    formatter = logging.Formatter(
        "[%(asctime)s - %(name)s - %(levelname)s] ::  %(message)s"
    )
    fh.setFormatter(formatter)
    fh.setLevel(logging.INFO)
    sh.setFormatter(formatter)
    sh.setLevel(logging.DEBUG)
    birdman_logger.addHandler(fh)
    birdman_logger.addHandler(sh)

    cmdstanpy_logger = cmdstanpy.utils.get_logger()
    cmdstanpy_logger.setLevel(logging.DEBUG)
    cmdstanpy_logger.addHandler(fh)
    for h in cmdstanpy_logger.handlers:
        h.setFormatter(formatter)

    return birdman_logger
