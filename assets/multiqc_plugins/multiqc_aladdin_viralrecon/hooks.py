#!/usr/bin/env python
"""
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_aladdin_rnaseq_version = get_distribution("multiqc_aladdin_viralrecon").version

# Add default config options that can be overriden by user config
def plugin_before_config():
    
    # Use the aladdin template by default
    config.template = 'aladdin'
    
# Add additional config options
def plugin_execution_start():
    """ Code to execute after the config files and
    command line flags have been parsed.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', False) is True:
        log.info("Plugging disable")
        return None

    log.info("Running viralrecon MultiQC Plugins v{}".format(config.multiqc_aladdin_viralrecon_version))

    # Prepend some additional filename cleaning
    config.fn_clean_exts[0:0] = [
        '_1.fastq.gz',
        '_2.fastq.gz',
    ]