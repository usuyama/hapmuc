# -*- coding: utf8 -*-

u"""
represents a base or indel with quality
"""

from __future__ import division
import settings
from common import log


class UnitBase(object):
    def __str__(self):
        return u"{0}: {1} {2} {3}".format(self.TYPE, self.seq, self.strand, self.quality)

    def __repr__(self):
        return self.__str__()

    def adjust(self):
        self.seq = self.seq.upper()
        if self.quality is not None and type(self.quality) is unicode:
            self.quality = ord(self.quality) - settings.base_quality_offset


class Base(UnitBase):
    TYPE = u"Base"

    def __init__(self, obs, strand, quality):
        self.seq = obs
        self.strand = strand
        self.quality = quality
        self.adjust()

    def key(self):
        return self.seq


class Ins(UnitBase):
    TYPE = u"Ins"

    def __init__(self, obs, strand, quality):
        self.seq = obs
        self.strand = strand
        self.quality = quality
        self.adjust()

    def key(self):
        return u"+{0}".format(self.seq)


class Del(UnitBase):
    TYPE = u"Del"

    def __init__(self, ref, strand, quality):
        self.seq = ref
        self.strand = strand
        self.quality = quality
        self.adjust()

    def key(self):
        return u"-{0}".format(self.seq)


def get_profiles(pileup_units):
    u"""
    input: array of pileup_units.{Base, Ins, Del}
    output: counts of bases and indels
    """
    profiles = dict()
    for pu in pileup_units:
        if pu.key() in profiles:
            profiles[pu.key()] += 1
        else:
            profiles[pu.key()] = 1
    return profiles


def filter_by_base_quality(pileup_units, base_quality_threshold):
    return [pu for pu in pileup_units if pu.quality is None or pu.quality >= base_quality_threshold]


def calc_avg_base_quality(pileup_units):
    total_count = len([x for x in pileup_units if x.quality is not None])
    if total_count == 0:
        return None
    sum = 0.0
    for pu in pileup_units:
        if pu.quality is not None:
            sum += pu.quality
    return sum / total_count


def calc_strand_freq(pileup_units):
    if len(pileup_units) == 0:
        return None
    num_forward_reads = len([pu for pu in pileup_units if pu.strand == 1])
    freq = float(num_forward_reads) / len(pileup_units)
    return freq

def get_depth_without_indel(pileup_units):
    return len([x for x in pileup_units if x.key()[0] not in u"+-"])
