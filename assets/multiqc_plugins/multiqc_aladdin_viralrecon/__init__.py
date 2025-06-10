#! /usr/bin/env python

from pkg_resources import get_distribution
from multiqc.utils import config

__version__ = get_distribution("multiqc_aladdin_viralrecon").version
config.multiqc_aladdin_viralrecon_version = __version__