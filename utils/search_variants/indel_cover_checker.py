# -*- coding: utf8 -*-

import settings
from common import log, IndelCoverError
import pileup_unit

class IndelCoverChecker:
    chromosome = None
    cover_position = None

    @classmethod
    def update(cls, current_chromosome, current_position, tumor_profiles, normal_profiles):
        # update cover-range
        total_counts = sum([tumor_profiles.get(base_key, 0) + normal_profiles.get(base_key, 0) for base_key in u"ATGCN"])
        indel_counts =  sum(tumor_profiles.values() + normal_profiles.values()) - total_counts
        if float(indel_counts) / total_counts >= settings.freq_threshold_for_indel_cover_checker:
            for key in set(tumor_profiles.keys() + normal_profiles.keys()):
                if key not in u"ATGCN":
                    length = len(key) - 1 # remove '+' or '-' symbole
                    if ((cls.chromosome is None)
                        or (cls.chromosome != current_chromosome)
                        or (cls.chromosome == current_chromosome and cls.cover_position < current_position + length)):
                        # log.debug("update indel cover: {0}:{1}".format(cls.chromosome, cls.cover_position))
                        cls.chromosome = current_chromosome
                        cls.cover_position = current_position + length

    @classmethod
    def check(cls, current_chromosome, current_position):
        if ((cls.chromosome is not None)
            and (cls.chromosome == current_chromosome)
            and (current_position <= cls.cover_position)):
            raise IndelCoverError(u"filter by indel-cover (<= {0}:{1})".format(cls.chromosome, cls.cover_position))
