# -*- coding: utf8 -*-

import pileup
import pileup_unit
import settings
import re
from io import open
import common
from common import log, CustomError, TooFewVariantReadsError, LowDepthError,\
    HighDepthError, TriallelicSiteError, TooManyNormalVariantReadsError, IndelCoverError,\
    AlleleFreqOutOfRangeError, StrandFreqOutOfRangeError, LowBaseQualityError
from variant import HeterozygousGermlineVariant, SomaticVariant
import triallelic_site_checker
from indel_cover_checker import IndelCoverChecker

REGEX_COUNT_W = re.compile("\w")


def get_variants_from_matched_lines(tumor_line, normal_line):
    u"""mathced tumor line and normal line"""

    hetero_germline_variants = []
    somatic_variants = []

    if (tumor_line.depth < settings.min_depth) or (normal_line.depth < settings.min_depth):
        raise LowDepthError

    if (tumor_line.depth > settings.max_depth) or (normal_line.depth > settings.max_depth):
        raise HighDepthError

    if tumor_line.ref == u'N':
        raise CustomError(u"reference_is_N")

    if len(REGEX_COUNT_W.findall(tumor_line.bases)) < settings.min_variant_supporting_reads:
        raise TooFewVariantReadsError

    tumor_pileup_units = tumor_line.get_bases_with_qualities()
    normal_pileup_units = normal_line.get_bases_with_qualities()

    tumor_profiles = pileup_unit.get_profiles(tumor_pileup_units)
    normal_profiles = pileup_unit.get_profiles(normal_pileup_units)

    for variant_key in tumor_profiles.keys():
        if variant_key == tumor_line.ref:
            # skip for the reference base
            continue
        try:
            tumor_count = tumor_profiles[variant_key]
            if tumor_count < settings.min_variant_supporting_reads:
                raise TooFewVariantReadsError

            tumor_ref_units = [x for x in tumor_pileup_units if x.key() != variant_key]
            tumor_obs_units = [x for x in tumor_pileup_units if x.key() == variant_key]
            normal_ref_units = [x for x in normal_pileup_units if x.key() != variant_key]
            normal_obs_units = [x for x in normal_pileup_units if x.key() == variant_key]

            IndelCoverChecker.update(tumor_line.chromosome, tumor_line.position, tumor_profiles, normal_profiles)

            try:
                normal_count = normal_profiles.get(variant_key, 0)
                if normal_count < settings.min_variant_supporting_reads:
                    raise TooFewVariantReadsError
                v = HeterozygousGermlineVariant.from_pileup_units(tumor_ref_units, tumor_obs_units, normal_ref_units, normal_obs_units)
                v.set_basic_info(variant_key, tumor_line.chromosome, tumor_line.position, tumor_line.ref)
                if v.is_snv():
                    try:
                        triallelic_site_checker.check(tumor_line.ref, tumor_line.chromosome, tumor_line.position,
                                                      tumor_profiles, normal_profiles)
                    except TriallelicSiteError:
                        v.triallelic_site_check = "triallelic"
                    try:
                        IndelCoverChecker.check(tumor_line.chromosome, tumor_line.position)
                    except IndelCoverError:
                        v.indel_cover_check = "indel-cover"
                hetero_germline_variants.append(v)
            except AlleleFreqOutOfRangeError: pass
            except StrandFreqOutOfRangeError: pass
            except TooFewVariantReadsError: pass
            except LowDepthError: pass
            except LowBaseQualityError as e:
                log.debug(u"HeteroGermline: {0}, tumor: {1}, normal: {2}".format(e, tumor_line, normal_line))
            except CustomError, e:
                log.warning(u"HeteroGermline CustomError: {0}, tumor: {1}, normal: {2}".format(e, tumor_line, normal_line))

            try:
                v = SomaticVariant.from_pileup_units(tumor_ref_units, tumor_obs_units, normal_ref_units, normal_obs_units)
                v.set_basic_info(variant_key, tumor_line.chromosome, tumor_line.position, tumor_line.ref)
                v.set_fisher_score()
                if v.is_snv():
                    try:
                        triallelic_site_checker.check(tumor_line.ref, tumor_line.chromosome, tumor_line.position,
                                                      tumor_profiles, normal_profiles)
                    except TriallelicSiteError:
                        v.triallelic_site_check = "triallelic"
                    try:
                        IndelCoverChecker.check(tumor_line.chromosome, tumor_line.position)
                    except IndelCoverError:
                        v.indel_cover_check = "indel-cover"
                somatic_variants.append(v)
            except AlleleFreqOutOfRangeError: pass
            except StrandFreqOutOfRangeError: pass
            except TooManyNormalVariantReadsError: pass
            except TooFewVariantReadsError: pass
            except LowDepthError: pass
            except LowBaseQualityError as e:
                log.debug(u"Somatic: {0}, tumor: {1}, normal: {2}".format(e, tumor_line, normal_line))
            except CustomError as e:
                log.warning(u"Somatic CustomError: {0}, tumor: {1}, normal: {2}".format(e, tumor_line, normal_line))

        except TooFewVariantReadsError: pass
        except CustomError as e:
            log.warning(u"CustomError: {0}, tumor: {1}, normal: {2}".format(e, tumor_line, normal_line))

    return [hetero_germline_variants, somatic_variants]


def search_variants(tumor_pileup_filename, normal_pileup_filename,
                    cand_somatic_variant_file, cand_hetero_germline_variant_file):
    tumor_f = open(tumor_pileup_filename, u'r')
    normal_f = open(normal_pileup_filename, u'r')

    normal_l = pileup.PileupLine(normal_f.readline())
    current_chromosome = normal_l.chromosome

    tumor_l = pileup.PileupLine(tumor_f.readline())
    if current_chromosome != tumor_l.chromosome:
        raise CustomError(u"different_chromosome_at_the_first_line")

    line_count = 0
    while True:

        # check how many lines have been processed
        line_count += 1
        if line_count % settings.debug_number_of_lines == 0:
            log.debug(u"""processing...
            \ttumor: {0}
            \tnormal: {1}""".format(tumor_l, normal_l))
            line_count = 0

        try:
            if tumor_l.chromosome != normal_l.chromosome:
                if normal_l.chromosome == current_chromosome:
                    normal_l = pileup.PileupLine(normal_f.readline())
                    continue
                else:
                    tumor_l = pileup.PileupLine(tumor_f.readline())
                    continue
            if tumor_l.position < normal_l.position:
                tumor_l = pileup.PileupLine(tumor_f.readline())
                continue
            elif tumor_l.position > normal_l.position:
                normal_l = pileup.PileupLine(normal_f.readline())
                continue
        except IOError as e:
            log.debug(u"reach the bottom of the file. {0}".format(e))
            break

        try:
            hetero_germline_results, somatic_results = get_variants_from_matched_lines(tumor_l, normal_l)

            [cand_somatic_variant_file.write(u"{0}\n".format(v)) for v in somatic_results]
            [cand_hetero_germline_variant_file.write(u"{0}\n".format(v)) for v in hetero_germline_results]

        except TooFewVariantReadsError: pass
        except LowDepthError: pass
        except HighDepthError: pass
        except CustomError as e:
            log.debug(u"CustomError: {0}, tumor: {1}, normal: {2}".format(e, tumor_l, normal_l))

        try:
            current_chromosome = normal_l.chromosome
            normal_l = pileup.PileupLine(normal_f.readline())
            tumor_l = pileup.PileupLine(tumor_f.readline())
        except IOError, e:
            log.debug(u"reach the bottom of the file. {0}".format(e))
            break

    tumor_f.close()
    normal_f.close()


def set_settings_by_opts(opts):
    common.log.debug(opts)
    settings.min_depth = opts.min_depth
    settings.max_depth = opts.max_depth
    settings.somatic_variant.min_tumor_freq = opts.min_tumor_freq
    settings.somatic_variant.sufficient_num_variant_reads = opts.sufficient_num_variant_reads
    settings.somatic_variant.max_normal_freq = opts.max_normal_freq
    settings.somatic_variant.max_variant_supporting_reads_normal = opts.max_variant_supporting_reads_normal


if __name__ == u"__main__":
    import sys
    (opts, args) = settings.parser.parse_args()
    set_settings_by_opts(opts)
    usage = u"""usage:
    python search_variants.py {tumor.pileup} {normal.pileup} {output_file_prefix} <options>
    \toutputs: {output_file_prefix}somatic_candidates
    \t\t{output_file_prefix}hetero_germline_variants"""
    if len(args) != 3:
        log.info(usage)
    else:
        tumor_pileup_filename = args[0]
        normal_pileup_filename = args[1]
        out_prefix = args[2]
        somatic_candidates_filename = u"{0}somatic_candidates".format(out_prefix)
        hetero_germline_candidates_filename = u"{0}hetero_germline_candidates".format(out_prefix)

        log.info(u"""
        inputs:
        \ttumor pileup file: {0}
        \tnormal pileup file: {1}
        outputs:
        \tcandidate somatic mutations: {2}
        \t candidate heterozygous germline variants: {3}""".format(tumor_pileup_filename, normal_pileup_filename,
                                                                   somatic_candidates_filename,
                                                                   hetero_germline_candidates_filename))

        log.info(u"\nsettings:\n" + settings.to_str())

        somatic_candidates_file = open(somatic_candidates_filename, u"w")
        hetero_germline_candidates_file = open(hetero_germline_candidates_filename, u"w")
        search_variants(tumor_pileup_filename, normal_pileup_filename,
                        somatic_candidates_file, hetero_germline_candidates_file)

        log.debug(u"done.")
