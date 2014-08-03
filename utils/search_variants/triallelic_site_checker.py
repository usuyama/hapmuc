# -*- coding: utf8 -*-

import settings
from common import log, TriallelicSiteError


def check(reference, current_chromosome, current_position, tumor_profiles, normal_profiles):
    total_counts = sum(tumor_profiles.values()) + sum(normal_profiles.values())
    count = 0
    for base in u"ATGC":
        read_counts = sum([tumor_profiles.get(base, 0) + normal_profiles.get(base, 0)])
        if (float(read_counts) / total_counts) >= settings.freq_threshold_for_triallelic_site_checker:
            count += 1
    if count > 2:
        raise TriallelicSiteError(u"triallelic sites filter")
