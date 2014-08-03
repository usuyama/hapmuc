# -*- coding: utf8 -*-

import sys
import logging

from colorlog import ColoredFormatter

formatter = ColoredFormatter(
        u"%(log_color)s[%(asctime)s] %(levelname)-8s %(message)s",
        datefmt=u"%m-%d %H:%M:%S",
        reset=True,
        log_colors={
                u'DEBUG':    u'cyan',
                u'INFO':     u'white',
                u'WARNING':  u'yellow',
                u'ERROR':    u'red',
                u'CRITICAL': u'red',
        }
)

log = logging.getLogger(u"search_variants.py")
log.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(formatter)
log.addHandler(handler)


class CustomError(Exception): pass
class TooFewVariantReadsError(Exception): pass
class LowDepthError(Exception): pass
class HighDepthError(Exception): pass
class TriallelicSiteError(Exception): pass
class TooManyNormalVariantReadsError(Exception): pass
class AlleleFreqOutOfRangeError(Exception): pass
class StrandFreqOutOfRangeError(Exception): pass
class LowBaseQualityError(Exception): pass
class IndelCoverError(Exception): pass
