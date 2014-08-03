# -*- coding: utf8 -*-

from __future__ import division
import pileup_unit
import settings
import fisher
import math

import common
from common import log, CustomError


def __my_formatter__(x):
    if x is None:
        return "-"
    elif type(x) == float:
        return str(round(x, 3))
    elif type(x) == int:
        return str(x)
    else:
        return x


class Variant:
    def __init__(self):
        self.chromosome = None
        self.start_pos = None
        self.end_pos = None
        self.ref = u"-"
        self.obs = u"-"
        self.tumor_ref_count = 0
        self.tumor_obs_count = 0
        self.normal_ref_count = 0
        self.normal_obs_count = 0
        self.tumor_freq = 0.0
        self.tumor_strand = 0.0
        self.normal_freq = 0.0
        self.normal_strand = 0.0
        self.tumor_ref_avg_base_quality = 0.0
        self.tumor_obs_avg_base_quality = 0.0
        self.normal_ref_avg_base_quality = 0.0
        self.normal_obs_avg_base_quality = 0.0
        self.indel_cover_check = "ok"
        self.triallelic_site_check = "ok"

    def to_str(self):
        out = [self.chromosome,
               self.start_pos,
               self.end_pos,
               self.ref,
               self.obs,
               self.tumor_ref_count,
               self.tumor_obs_count,
               self.normal_ref_count,
               self.normal_obs_count,
               self.tumor_freq,
               self.tumor_strand,
               self.normal_freq,
               self.normal_strand,
               self.tumor_ref_avg_base_quality,
               self.tumor_obs_avg_base_quality,
               self.normal_ref_avg_base_quality,
               self.normal_obs_avg_base_quality,
               self.indel_cover_check,
               self.triallelic_site_check]
        return "\t".join([__my_formatter__(x) for x in out])

    def is_snv(self):
        return not (self.ref == u"-" or self.obs == u"-")

    def set_basic_info(self, variant_key, chromosome, position, reference_base):
        self.chromosome = chromosome
        if variant_key[0] == u"+":
            self.ref = u"-"
            self.obs = variant_key[1:]
            self.start_pos = position
            self.end_pos = position
        elif variant_key[0] == u"-":
            self.ref = variant_key[1:]
            self.obs = u"-"
            self.start_pos = position
            self.end_pos = position + len(self.ref)
        else:
            self.ref = reference_base
            self.obs = variant_key
            self.start_pos = position - 1
            self.end_pos = position



class HeterozygousGermlineVariant(Variant):
    def __init__(self):
        Variant.__init__(self)

    @classmethod
    def from_pileup_units(cls, raw_tumor_ref_units, raw_tumor_obs_units, raw_normal_ref_units, raw_normal_obs_units):
        v = HeterozygousGermlineVariant()
        params = settings.hetero_germline_variant
        variant_key = raw_tumor_obs_units[0].key()

        tumor_ref_units = pileup_unit.filter_by_base_quality(raw_tumor_ref_units, settings.base_quality_threshold)
        tumor_obs_units = pileup_unit.filter_by_base_quality(raw_tumor_obs_units, settings.base_quality_threshold)
        normal_ref_units = pileup_unit.filter_by_base_quality(raw_normal_ref_units, settings.base_quality_threshold)
        normal_obs_units = pileup_unit.filter_by_base_quality(raw_normal_obs_units, settings.base_quality_threshold)

        v.tumor_ref_count = pileup_unit.get_depth_without_indel(tumor_ref_units)
        v.normal_ref_count = pileup_unit.get_depth_without_indel(normal_ref_units)
        v.tumor_obs_count = len(tumor_obs_units)
        v.normal_obs_count = len(normal_obs_units)

        if not (v.tumor_obs_count >= settings.min_variant_supporting_reads):
            raise common.TooFewVariantReadsError

        if not (v.normal_obs_count >= settings.min_variant_supporting_reads):
            raise common.TooFewVariantReadsError

        tumor_depth = pileup_unit.get_depth_without_indel(tumor_ref_units + tumor_obs_units)
        if tumor_depth < settings.min_depth:
            raise common.LowDepthError

        v.tumor_freq = float(v.tumor_obs_count) / tumor_depth
        if not (params.min_freq <= v.tumor_freq <= params.max_freq):
            raise common.AlleleFreqOutOfRangeError

        v.tumor_strand = pileup_unit.calc_strand_freq(tumor_obs_units)
        if not (params.min_strand_freq <= v.tumor_strand <= params.max_strand_freq):
            raise common.StrandFreqOutOfRangeError

        normal_depth = pileup_unit.get_depth_without_indel(normal_ref_units + normal_obs_units)
        if normal_depth < settings.min_depth:
            raise common.LowDepthError

        v.normal_freq = float(v.normal_obs_count) / normal_depth
        if not (params.min_freq <= v.normal_freq <= params.max_freq):
            raise common.AlleleFreqOutOfRangeError

        v.normal_strand = pileup_unit.calc_strand_freq(normal_obs_units)
        if not (params.min_strand_freq <= v.normal_strand <= params.max_strand_freq):
            raise common.StrandFreqOutOfRangeError

        v.tumor_ref_avg_base_quality = pileup_unit.calc_avg_base_quality(raw_tumor_ref_units)
        v.tumor_obs_avg_base_quality = pileup_unit.calc_avg_base_quality(raw_tumor_obs_units)
        v.normal_ref_avg_base_quality = pileup_unit.calc_avg_base_quality(raw_normal_ref_units)
        v.normal_obs_avg_base_quality = pileup_unit.calc_avg_base_quality(raw_normal_obs_units)

        for avg_bq in [v.tumor_ref_avg_base_quality,
                       v.tumor_obs_avg_base_quality,
                       v.normal_ref_avg_base_quality,
                       v.normal_obs_avg_base_quality]:
            if avg_bq is not None and not (params.min_avg_base_quality <= avg_bq):
                raise common.LowBaseQualityError(u"low_avg_base_quality({0} < {1})".format(avg_bq, params.min_avg_base_quality))

        return v

    def __str__(self):
        return Variant.to_str(self)


class SomaticVariant(Variant):
    def __init__(self):
        Variant.__init__(self)
        self.fisher_score = None

    def set_fisher_score(self):
        # https://pypi.python.org/pypi/fisher/
        [[a, b], [c, d]] = [[self.tumor_ref_count, self.tumor_obs_count],
                            [self.normal_ref_count, self.normal_obs_count]]
        p = fisher.pvalue(a, b, c, d)
        self.fisher_score = - math.log10(p.two_tail)

    @classmethod
    def from_pileup_units(cls, raw_tumor_ref_units, raw_tumor_obs_units, raw_normal_ref_units, raw_normal_obs_units):
        v = SomaticVariant()
        params = settings.somatic_variant

        tumor_ref_units = pileup_unit.filter_by_base_quality(raw_tumor_ref_units, settings.base_quality_threshold)
        tumor_obs_units = pileup_unit.filter_by_base_quality(raw_tumor_obs_units, settings.base_quality_threshold)
        normal_ref_units = pileup_unit.filter_by_base_quality(raw_normal_ref_units, settings.base_quality_threshold)
        normal_obs_units = pileup_unit.filter_by_base_quality(raw_normal_obs_units, settings.base_quality_threshold)

        v.tumor_ref_count = pileup_unit.get_depth_without_indel(tumor_ref_units)
        v.normal_ref_count = pileup_unit.get_depth_without_indel(normal_ref_units)
        v.tumor_obs_count = len(tumor_obs_units)
        v.normal_obs_count = len(normal_obs_units)

        if not (v.tumor_obs_count >= settings.min_variant_supporting_reads):
            raise common.TooFewVariantReadsError

        if not (v.normal_obs_count <= params.max_variant_supporting_reads_normal):
            raise common.TooManyNormalVariantReadsError

        tumor_depth = pileup_unit.get_depth_without_indel(tumor_ref_units + tumor_obs_units)
        if tumor_depth < settings.min_depth:
            raise common.LowDepthError


        v.tumor_freq = float(v.tumor_obs_count) / tumor_depth
        if not (params.min_tumor_freq <= v.tumor_freq):
            if params.sufficient_num_variant_reads is not None\
                and v.tumor_obs_count < params.sufficient_num_variant_reads:
                raise common.AlleleFreqOutOfRangeError

        v.tumor_strand = pileup_unit.calc_strand_freq(tumor_obs_units)
        if not (params.min_strand_freq <= v.tumor_strand <= params.max_strand_freq):
            raise common.StrandFreqOutOfRangeError

        normal_depth = pileup_unit.get_depth_without_indel(normal_ref_units + normal_obs_units)
        if normal_depth < settings.min_depth:
            raise common.LowDepthError

        v.normal_freq = float(v.normal_obs_count) / normal_depth
        if not (v.normal_freq <= params.max_normal_freq):
            raise common.AlleleFreqOutOfRangeError

        v.normal_strand = pileup_unit.calc_strand_freq(normal_obs_units)

        v.tumor_ref_avg_base_quality = pileup_unit.calc_avg_base_quality(raw_tumor_ref_units)
        v.tumor_obs_avg_base_quality = pileup_unit.calc_avg_base_quality(raw_tumor_obs_units)
        v.normal_ref_avg_base_quality = pileup_unit.calc_avg_base_quality(raw_normal_ref_units)
        v.normal_obs_avg_base_quality = pileup_unit.calc_avg_base_quality(raw_normal_obs_units)

        for avg_bq in [v.tumor_ref_avg_base_quality,
                       v.tumor_obs_avg_base_quality,
                       v.normal_ref_avg_base_quality]:
            if avg_bq is not None and not (params.min_avg_base_quality <= avg_bq):
                raise common.LowBaseQualityError(u"low_avg_base_quality({0} < {1}".format(avg_bq, params.min_avg_base_quality))

        return v

    def __str__(self):
        out = Variant.to_str(self)
        out += "\t" + __my_formatter__(self.fisher_score)
        return out
